"""HPC cluster config & SGE job queue management.

Ported from auto_opt (`auto_opt/cluster.py`), dropping the AmberTools-specific
`get_amber_tool()` (csp does not use Amber — see CLAUDE.md). Everything else
here (queue config, qstat polling, SGE job script generation) is generic to
running Gaussian16 jobs on the same SGE cluster.

Reads user environment config from ~/.csp.yaml (falls back to built-in
defaults if absent):

  scheduler: sge
  queues:
    - name: gr1.q
      nproc: 40
      pe: OpenMP
    - name: gr2.q
      nproc: 52
      pe: OpenMP
  max_concurrent_jobs: 6
  poll_interval: 30
  nproc_reserve: 2
"""
from __future__ import annotations

import os
import re
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

_DEFAULTS: dict = {
    'scheduler': 'sge',
    'queues': [
        {'name': 'gr1.q', 'nproc': 40, 'pe': 'OpenMP'},
        {'name': 'gr2.q', 'nproc': 52, 'pe': 'OpenMP'},
    ],
    'max_concurrent_jobs': 6,
    'poll_interval': 30,
    'nproc_reserve': 2,
}

_cached: dict | None = None


def load_env(config_path: str | Path | None = None) -> dict:
    """Return environment config.

    Priority: config_path argument > ~/.csp.yaml > built-in defaults.
    Cached after the first call.
    """
    global _cached
    if _cached is not None:
        return _cached

    cfg = dict(_DEFAULTS)
    search = [Path.home() / '.csp.yaml']
    if config_path:
        search = [Path(config_path)] + search

    for p in search:
        if p.exists():
            with open(p, encoding='utf-8') as f:
                user = yaml.safe_load(f) or {}
            cfg.update(user)
            break

    _cached = cfg
    return cfg


def _queue_spec() -> Dict[str, dict]:
    return {q['name']: q for q in load_env()['queues']}


def _queue_names() -> List[str]:
    return [q['name'] for q in load_env()['queues']]


def get_my_job_count() -> int:
    """Number of this user's SGE jobs currently running or queued."""
    user = os.environ.get('USER', '')
    try:
        r = subprocess.run(['qstat', '-u', user],
                           capture_output=True, text=True, timeout=15)
        return sum(1 for ln in r.stdout.splitlines() if re.match(r'^\s*\d+', ln))
    except Exception:
        return 0


def get_free_queue_instances() -> Dict[str, List[str]]:
    """Parse `qstat -f` and return queue instances with used==0 slots."""
    names = _queue_names()
    free: Dict[str, List[str]] = {n: [] for n in names}
    q_pat = '|'.join(re.escape(n) for n in names)
    pat = re.compile(rf'^(({q_pat})@\S+)\s+\S+\s+\d+/(\d+)/\d+', re.MULTILINE)
    try:
        r = subprocess.run(['qstat', '-f', '-u', '*'],
                           capture_output=True, text=True, timeout=15)
        for m in pat.finditer(r.stdout):
            qi, qname, used = m.group(1), m.group(2), int(m.group(3))
            if used == 0:
                free[qname].append(qi)
    except Exception:
        pass
    return free


def pick_free_instance(
    free_by_queue: Dict[str, List[str]],
    prefer_queue: Optional[str] = None,
) -> Optional[Tuple[str, str]]:
    """Pick one free queue instance, returning (qname, queue_instance)."""
    names = _queue_names()
    order = names if not prefer_queue else \
            [prefer_queue] + [q for q in names if q != prefer_queue]
    for qname in order:
        if free_by_queue.get(qname):
            return qname, free_by_queue[qname][0]
    return None


def wait_for_free_node(
    prefer_queue: Optional[str] = None,
    is_test: bool = False,
    test_index: int = 0,
) -> Tuple[str, str, int]:
    """Block until a free node is available, returning (qname, queue_instance, nproc_actual).

    nproc_actual = nproc - nproc_reserve (pass this as the job's process count).
    With is_test=True, skips qstat and returns a deterministic dummy value.
    """
    cfg = load_env()
    spec = _queue_spec()
    names = _queue_names()
    max_conc = cfg['max_concurrent_jobs']
    poll = cfg['poll_interval']
    reserve = cfg['nproc_reserve']

    if is_test:
        qname = names[test_index % len(names)]
        nproc = spec[qname]['nproc']
        return qname, f"{qname}@node{test_index + 1:02d}", nproc - reserve

    while True:
        n_jobs = get_my_job_count()
        if n_jobs < max_conc:
            free = get_free_queue_instances()
            pick = pick_free_instance(free, prefer_queue)
            if pick:
                qname, qi = pick
                nproc = spec[qname]['nproc']
                return qname, qi, nproc - reserve
            print(f"  No free node. Rechecking in {poll}s...")
        else:
            print(f"  {n_jobs} jobs running (limit {max_conc}). Rechecking in {poll}s...")
        time.sleep(poll)


def make_job_script(
    cmd: str,
    queue: str,
    queue_instance: Optional[str] = None,
    job_name: Optional[str] = None,
    stdout: Optional[str] = None,
    stderr: Optional[str] = None,
) -> str:
    """Generate an SGE job script string that runs `cmd`."""
    spec = _queue_spec()
    q = spec[queue]
    nproc = q['nproc']
    pe = q.get('pe', 'OpenMP')
    q_directive = queue_instance if queue_instance else queue

    lines = ['#!/bin/sh\n', '#$ -S /bin/sh\n', '#$ -cwd\n', '#$ -V\n']
    if job_name:
        lines.append(f'#$ -N {job_name}\n')
    lines.append(f'#$ -q {q_directive}\n')
    lines.append(f'#$ -pe {pe} {nproc}\n')
    if stdout:
        lines.append(f'#$ -o {stdout}\n')
    if stderr:
        lines.append(f'#$ -e {stderr}\n')
    lines += ['\n', 'hostname\n', '\n', cmd + '\n', '\n']
    return ''.join(lines)
