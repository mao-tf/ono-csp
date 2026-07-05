#!/usr/bin/env python3
"""
Step 1a (post-processing): Extract candidate (a, b) pairs from vdW contact sweep.

Reads vdW_r_contact_{monomer}.csv produced by sweep.py and finds, for each alpha:
  - a-stack: theta_c at upper end of the valid (TorF=True) region  → a is minimized
  - b-stack: theta_c at lower end of the valid region               → b is minimized
  - local_min (optional): local minima of S = a*b within valid region

Output columns: alpha, theta_c, a, b, S, structure_type

Usage:
  python -m csp.vdw.extract_init \\
      --vdw-csv runs/anthracene/vdW_r_contact_anthracene.csv \\
      --out     runs/anthracene/step1_init_params.csv

  # Include local minima of S=a*b:
  python -m csp.vdw.extract_init --vdw-csv ... --out ... --minima

  # Only a-stack:
  python -m csp.vdw.extract_init --vdw-csv ... --out ... --select a-stack
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _local_minima(values: np.ndarray) -> list[int]:
    """Indices of interior local minima (endpoints excluded)."""
    idx = []
    n = len(values)
    for i in range(1, n - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            if values[i] < values[i - 1] or values[i] < values[i + 1]:
                idx.append(i)
    return idx


def _true_runs(mask: np.ndarray) -> list[list[int]]:
    """Return contiguous runs of True in a boolean array (as lists of indices)."""
    runs: list[list[int]] = []
    i = 0
    while i < len(mask):
        if not mask[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(mask) and mask[j + 1]:
            j += 1
        runs.append(list(range(i, j + 1)))
        i = j + 1
    return runs


# ---------------------------------------------------------------------------
# Main extraction
# ---------------------------------------------------------------------------

def extract_init(
    vdw_csv: str | Path,
    out_csv: str | Path,
    include_minima: bool = False,
    select: list[str] | None = None,
    round_digits: int = 2,
) -> pd.DataFrame:
    """Extract candidate structures from vdW contact sweep CSV.

    Parameters
    ----------
    vdw_csv : path to vdW_r_contact_{monomer}.csv (from sweep.py)
    out_csv : output path for step1_init_params.csv
    include_minima : if True, also include local S=a*b minima
    select : list of structure types to keep ('a-stack', 'b-stack', 'local_min')
             None or ['all'] → keep everything
    round_digits : decimal places for a, b values
    """
    select_set = set(select) if select else {'all'}
    accept_all = 'all' in select_set

    df = pd.read_csv(vdw_csv)
    required = {'alpha', 'theta_c', 'R_clps', 'TorF'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {vdw_csv}: {missing}")

    # Compute geometry
    theta_rad = np.deg2rad(df['theta_c'].to_numpy(float))
    R = df['R_clps'].to_numpy(float)
    df['a'] = 2.0 * R * np.cos(theta_rad)
    df['b'] = 2.0 * R * np.sin(theta_rad)
    df['S'] = df['a'] * df['b']

    # TorF may be stored as bool or string
    if df['TorF'].dtype != bool:
        df['TorF'] = df['TorF'].astype(str).str.lower().isin({'true', '1', 't', 'y', 'yes'})

    rows = []
    for alpha, grp in df.groupby('alpha', sort=True):
        grp = grp.sort_values('theta_c').reset_index(drop=True)
        mask = grp['TorF'].to_numpy(bool)
        if not mask.any():
            continue

        for run in _true_runs(mask):
            candidates: dict[int, str] = {}
            candidates[run[0]]  = 'b-stack'   # smallest theta_c → b is minimized
            candidates[run[-1]] = 'a-stack'   # largest theta_c  → a is minimized

            if include_minima:
                S_vals = grp.loc[run, 'S'].to_numpy(float)
                for k in _local_minima(S_vals):
                    idx = run[k]
                    if idx not in candidates:
                        candidates[idx] = 'local_min'

            for i in sorted(candidates):
                stype = candidates[i]
                if not accept_all and stype not in select_set:
                    continue
                a = round(float(grp.loc[i, 'a']), round_digits)
                b = round(float(grp.loc[i, 'b']), round_digits)
                S = round(a * b, round_digits)
                tc = float(grp.loc[i, 'theta_c'])
                rows.append({
                    'alpha': float(alpha),
                    'theta_c': tc,
                    'a': a,
                    'b': b,
                    'S': S,
                    'structure_type': stype,
                    'status': 'NotYet',
                })

    if rows:
        out_df = (
            pd.DataFrame(rows)
            .drop_duplicates(subset=['alpha', 'a', 'b', 'structure_type'])
            .sort_values(['alpha', 'a', 'b'])
            .reset_index(drop=True)
        )
    else:
        print("Warning: no candidates found with current settings.")
        out_df = pd.DataFrame(columns=['alpha', 'theta_c', 'a', 'b', 'S',
                                        'structure_type', 'status'])

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}  (n={len(out_df)})")
    return out_df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Step 1a: extract candidate (a,b) pairs from vdW contact sweep"
    )
    ap.add_argument('--vdw-csv', required=True,
                    help='vdW_r_contact_{monomer}.csv from sweep.py')
    ap.add_argument('--out', required=True,
                    help='Output CSV path (step1_init_params.csv)')
    ap.add_argument('--minima', action='store_true',
                    help='Also include local minima of S=a×b')
    ap.add_argument('--select', nargs='+',
                    default=['all'],
                    choices=['all', 'a-stack', 'b-stack', 'local_min'],
                    help='Structure types to keep (default: all)')
    ap.add_argument('--round', type=int, default=2, dest='round_digits',
                    help='Decimal places for a, b (default: 2)')

    args = ap.parse_args()
    extract_init(
        vdw_csv=args.vdw_csv,
        out_csv=args.out,
        include_minima=args.minima,
        select=args.select,
        round_digits=args.round_digits,
    )


if __name__ == '__main__':
    main()
