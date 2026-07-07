"""Reconstruct the paper's Fig. 2b-style branch classification from a
Step 1 DFT-D hill-climb history (step1.csv).

step1.csv only records every (theta, a, b) point Ohno's hill-climb visited
plus its energy -- not which of the 2-3 seeds (per theta) each point
belongs to, nor which point a given seed's climb converged to. Both are
recoverable without any new DFT calculation:

- The seeds and their kind (a-stack / b-stack / local min) are exactly
  `csp.vdw.contact.step1a_scan`'s `df_init` (same vdW rough-scan Ohno's own
  `init_process` uses, verified to produce identical (a, b, theta) values).
- Replaying the same greedy 3x3-neighbor descent step1.py's
  `get_opt_params_dict` uses, but against step1.csv's already-computed
  (theta, a, b) -> E lookup instead of launching new Gaussian jobs,
  deterministically finds the same converged point the real hill-climb did.

Polyacenes are also exactly symmetric under theta -> 90-theta with a and b
swapped (verified against the vdW contact model: vdw_R(theta) for 'a' equals
vdw_R(90-theta) for 'b'), since it's the same physical structure with the a/b
labels exchanged. Reflecting the computed branches through this symmetry
extends a theta in [5, 45] scan to the full [5, 85] range shown in the
paper's Fig. 2b, with no extra computation.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd

from csp.structure.molecule import Molecule
from csp.vdw.contact import step1a_scan

KIND_LABEL = {"b_contact": "b-stack", "a_contact": "a-stack", "local_min": "local min"}
KIND_COLOR = {"b_contact": "#1f77b4", "a_contact": "#d62728", "local_min": "#2ca02c"}
_KIND_SWAP = {"b_contact": "a_contact", "a_contact": "b_contact", "local_min": "local_min"}


def _replay_hill_climb(
    lookup: Mapping[tuple, float], theta: float, a0: float, b0: float,
) -> tuple[float, float, Optional[float]]:
    """Greedy 3x3-neighbor descent identical to step1.py's
    get_opt_params_dict, but reading energies from `lookup` instead of
    submitting new Gaussian jobs. Terminates at whatever step1.csv shows
    the real hill-climb converged to (every point it could step to was,
    by construction, already computed there).
    """
    theta = round(theta, 1)
    a, b = round(a0, 1), round(b0, 1)
    while True:
        candidates = []
        for da in (-0.1, 0.0, 0.1):
            for db in (-0.1, 0.0, 0.1):
                key = (theta, round(a + da, 1), round(b + db, 1))
                if key in lookup:
                    candidates.append((lookup[key], key[1], key[2]))
        if not candidates:
            return a, b, lookup.get((theta, a, b))
        best_E, best_a, best_b = min(candidates)
        if (best_a, best_b) == (a, b):
            return a, b, best_E
        a, b = best_a, best_b


def classify_and_fold_step1_results(
    mol: Molecule,
    alphas: Iterable[float],
    df_results: pd.DataFrame,
    radii_overrides: Optional[Mapping[str, float]] = None,
) -> pd.DataFrame:
    """Branch-classified, theta<->90-theta folded Step 1 DFT-D results.

    `df_results` is step1.csv (columns a, b, theta, E, ...). Returns one row
    per (seed, fold) with columns theta, a, b, E, kind ('b_contact' /
    'a_contact' / 'local_min'), folded (bool) -- ready for a Fig. 2b-style
    plot colored by `kind`.
    """
    _, df_init = step1a_scan(mol, alphas, radii_overrides=radii_overrides)

    lookup = {
        (round(r.theta, 1), round(r.a, 1), round(r.b, 1)): r.E
        for r in df_results.itertuples()
    }

    rows = []
    converged = {}  # theta -> {(a_f, b_f) already claimed by a-stack/b-stack}
    # Process a/b-stack endpoints before local_min so a local_min that
    # converges onto the exact same point (a common outcome once the
    # hill-climb actually runs -- the "seed" only picks the starting point,
    # not a separate final structure) can be recognized as redundant and
    # dropped, instead of duplicating an existing branch and creating a
    # jump discontinuity where it happens to differ.
    order = {"b_contact": 0, "a_contact": 0, "local_min": 1}
    for seed in sorted(df_init.itertuples(), key=lambda s: order.get(s.kind, 1)):
        a_f, b_f, E_f = _replay_hill_climb(lookup, seed.theta, seed.a, seed.b)
        if E_f is None:
            continue
        theta_r = round(seed.theta, 1)
        point = (a_f, b_f)
        if seed.kind == "local_min" and point in converged.get(theta_r, set()):
            continue
        converged.setdefault(theta_r, set()).add(point)
        rows.append({"theta": theta_r, "a": a_f, "b": b_f, "E": E_f,
                      "kind": seed.kind, "folded": False})

    df_branches = pd.DataFrame(rows)
    if len(df_branches) == 0:
        return df_branches

    df_folded = df_branches.copy()
    df_folded["theta"] = 90.0 - df_folded["theta"]
    df_folded["a"], df_folded["b"] = df_branches["b"].values, df_branches["a"].values
    df_folded["kind"] = df_folded["kind"].map(_KIND_SWAP)
    df_folded["folded"] = True

    return pd.concat([df_branches, df_folded], ignore_index=True).sort_values(
        ["kind", "theta"]
    ).reset_index(drop=True)
