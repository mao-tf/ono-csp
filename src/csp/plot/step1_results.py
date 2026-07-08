"""Reconstruct the paper's Fig. 2b-style branch classification from a
Step 1 DFT-D hill-climb history (step1.csv).

step1.csv only records every (theta, a, b) point visited plus its energy
-- not which structural branch (a-stack / b-stack / local min) each point
belongs to, nor which point each branch converged to. Both are recovered
directly from the data itself, with no dependency on re-running any vdW
model or knowing which molecule/monomer produced the CSV:

- At each theta, pivot the visited (a, b) into a 2D energy grid and find
  its local minima (`scipy.ndimage.minimum_filter`, the same approach
  used for the other 2D maps in this app) -- these *are* the converged
  branches, directly, since the DFT hill-climb only ever stops at a local
  minimum of its own sampled grid.
- Label them by sorting on `a`: the smallest-a minimum is 'a_contact',
  the largest-a minimum is 'b_contact' (matching the vdW pre-scan's own
  naming for the two beta-sweep endpoints), and anything in between is
  'local_min' -- this reproduces step1a_scan's kind labels without
  calling it.

Polyacenes are also exactly symmetric under theta -> 90-theta with a and b
swapped (same physical structure, axes relabeled). Reflecting the
computed branches through this symmetry extends a theta in [5, 45] scan
to the full [5, 85] range shown in the paper's Fig. 2b, with no extra
computation. `kind` is kept fixed across this fold (not re-derived from
the swapped a/b) precisely because it's a *physical*-identity label
(HB/PS, via KIND_LABEL) -- the smallest-a/largest-a criterion that
assigned it in the first place is geometric and would flip under the a/b
swap, but the branch itself doesn't change identity just because its axes
got relabeled. For the same reason, KIND_LABEL only says "HB"/"PS", not
"HB (b-stack)"/"PS (a-stack)" -- the a-stack/b-stack (smallest-a/
largest-a) description is only true for the unfolded half (theta <= 45);
past the fold it's backwards, so it's dropped from the display label
entirely rather than shown-and-wrong for half the range.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import minimum_filter

KIND_LABEL = {"b_contact": "HB", "a_contact": "PS", "local_min": "local min"}
KIND_COLOR = {"b_contact": "#1f77b4", "a_contact": "#d62728", "local_min": "#2ca02c"}


def _local_minima_at_theta(df_theta: pd.DataFrame) -> pd.DataFrame:
    """Local minima of E over the visited (a, b) grid at one fixed theta,
    labeled by position (smallest a -> a_contact, largest a -> b_contact,
    else local_min).
    """
    pivot = df_theta.pivot_table(index="b", columns="a", values="E", aggfunc="min")
    grid = pivot.values
    is_min = (grid == minimum_filter(grid, size=3)) & np.isfinite(grid)
    bi, ai = np.where(is_min)
    a_vals = pivot.columns.to_numpy()[ai]
    b_vals = pivot.index.to_numpy()[bi]
    e_vals = grid[bi, ai]

    order = np.argsort(a_vals)
    a_vals, b_vals, e_vals = a_vals[order], b_vals[order], e_vals[order]

    n = len(a_vals)
    if n == 0:
        return pd.DataFrame(columns=["a", "b", "E", "kind"])
    kinds = ["local_min"] * n
    kinds[0] = "a_contact"
    kinds[-1] = "b_contact"  # if n==1 this overwrites index 0, giving a single b_contact point
    return pd.DataFrame({"a": a_vals, "b": b_vals, "E": e_vals, "kind": kinds})


def classify_and_fold_step1_results(df_results: pd.DataFrame) -> pd.DataFrame:
    """Branch-classified, theta<->90-theta folded Step 1 DFT-D results.

    `df_results` is step1.csv (columns a, b, theta, E, ...). Returns one row
    per (branch, fold) with columns theta, a, b, E, kind ('b_contact' /
    'a_contact' / 'local_min'), folded (bool) -- ready for a Fig. 2b-style
    plot colored by `kind`.
    """
    rows = []
    for theta, df_theta in df_results.groupby("theta"):
        branches = _local_minima_at_theta(df_theta)
        for r in branches.itertuples():
            rows.append({"theta": round(float(theta), 1), "a": float(r.a), "b": float(r.b),
                         "E": float(r.E), "kind": r.kind, "folded": False})

    df_branches = pd.DataFrame(rows)
    if len(df_branches) == 0:
        return df_branches

    df_folded = df_branches.copy()
    df_folded["theta"] = 90.0 - df_folded["theta"]
    df_folded["a"], df_folded["b"] = df_branches["b"].values, df_branches["a"].values
    # `kind` is NOT re-derived from the swapped (a, b) here: folding a branch
    # to theta -> 90-theta with a/b swapped is the *same physical structure*
    # under axis relabeling (HB stays HB, PS stays PS), even though the
    # a_contact/b_contact criterion itself ("smallest-a" vs "largest-a") is
    # purely geometric and would flip if reapplied to the swapped values.
    # Keeping `kind` fixed is what makes it a stable physical-identity label
    # (and what KIND_LABEL's HB/PS mapping assumes) across the full fold.
    df_folded["folded"] = True

    return pd.concat([df_branches, df_folded], ignore_index=True).sort_values(
        ["kind", "theta"]
    ).reset_index(drop=True)
