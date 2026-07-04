"""Rigid-sphere VdW contact distance and the Step 1a coarse (a, b) scan.

`vdw_contact_distance` is the vectorized equivalent of Ono's `vdw_R`
(legacy/ono_scripts/stepwise_optimization/vdw.py) and auto_opt's `vdw_R`
(both compute the same quantity with O(N^2) Python loops).

`step1a_scan` ports Ono's `init_process`/`get_init_para_csv`
(legacy/ono_scripts/stepwise_optimization/step1.py) exactly:

For each herringbone half-angle alpha (= A3 = the `theta` column of the
legacy CSVs):
  1. a_clps / b_clps: contact distance of two parallel (+alpha) molecules
     pushed along the a-axis (0°) / b-axis (90°) — the smallest possible
     lattice constants from slipped-parallel contact.
  2. Sweep the T-contact direction theta_ab over 0..90°: push a -alpha
     molecule along theta_ab until VdW contact at distance R, giving the
     lattice (a, b) = 2R(cos, sin) that puts the T neighbor at (a/2, b/2).
  3. Keep points with a >= a_clps and b >= b_clps (SP neighbors must not
     overlap); record S = a*b.
  4. Initial candidates for the DFT step: local minima of S along the sweep
     (scipy.signal.argrelmin, order=5) plus both endpoints (b-contact-limited
     and a-contact-limited arrangements).
"""
from __future__ import annotations

import math
from typing import Callable, Iterable, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import signal

from csp.structure.molecule import Molecule
from csp.structure.intralayer import to_layer_frame, place_in_layer


def vdw_contact_distance(
    coords_1: np.ndarray, radii_1: np.ndarray,
    coords_2: np.ndarray, radii_2: np.ndarray,
    theta_deg: float,
) -> float:
    """Push molecule 2 along the xy-direction `theta_deg` until VdW contact.

    Returns the minimum push distance (Å) such that no atom pair of the two
    rigid-sphere molecules overlaps; 0.0 if they are already clear.
    """
    ct = math.cos(math.radians(theta_deg))
    st = math.sin(math.radians(theta_deg))
    eR = np.array([ct, st, 0.0], float)

    D = coords_2[None, :, :] - coords_1[:, None, :]
    R12b = D @ eR
    R12a2 = (D * D).sum(axis=2) - R12b ** 2

    rad_sum_sq = (radii_1[:, None] + radii_2[None, :]) ** 2
    mask = R12a2 < rad_sum_sq
    if not np.any(mask):
        return 0.0

    push = -R12b[mask] + np.sqrt(rad_sum_sq[mask] - R12a2[mask])
    return float(max(np.max(push), 0.0))


def _effective_radii(mol: Molecule, overrides: Optional[Mapping[str, float]]) -> np.ndarray:
    if not overrides:
        return mol.radii
    table = {str(k).strip().upper(): float(v) for k, v in overrides.items()}
    return np.array(
        [table.get(s.strip().upper(), r) for s, r in zip(mol.symbols, mol.radii)],
        dtype=float,
    )


def step1a_scan(
    mol: Molecule,
    alphas: Iterable[float],
    radii_overrides: Optional[Mapping[str, float]] = None,
    theta_ab_step: float = 1.0,
    argrelmin_order: int = 5,
    progress_callback: Optional[Callable[[float], None]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Step 1a coarse scan (see module docstring).

    Returns (df_curves, df_init):
    - df_curves: every VdW-feasible sweep point; columns
      alpha, theta_ab, a, b, S (a, b rounded to 0.1 Å as in the legacy code;
      S computed from the unrounded values).
    - df_init: initial candidates for the DFT hill-climb; columns
      a, b, theta, S, kind ('local_min' / 'b_contact' / 'a_contact'),
      status='NotYet' — matching the legacy step1_init_params.csv format
      (theta = alpha).
    """
    mol_l = to_layer_frame(mol)
    radii = _effective_radii(mol_l, radii_overrides)
    alphas = list(alphas)
    theta_abs = np.arange(0.0, 90.0 + 1e-9, theta_ab_step)

    curve_rows = []
    init_rows = []
    for i_alpha, alpha in enumerate(alphas):
        c_i = place_in_layer(mol_l.coords, 0, 0, 0, 0.0, +alpha)
        c_t = place_in_layer(mol_l.coords, 0, 0, 0, 0.0, -alpha)
        a_clps = vdw_contact_distance(c_i, radii, c_i, radii, 0.0)
        b_clps = vdw_contact_distance(c_i, radii, c_i, radii, 90.0)

        kept = []  # (theta_ab, a_rounded, b_rounded, S)
        for theta_ab in theta_abs:
            R = vdw_contact_distance(c_i, radii, c_t, radii, theta_ab)
            a = 2 * R * math.cos(math.radians(theta_ab))
            b = 2 * R * math.sin(math.radians(theta_ab))
            if (a_clps > a) or (b_clps > b):
                continue
            kept.append((theta_ab, round(a, 1), round(b, 1), a * b))

        for theta_ab, a1, b1, S in kept:
            curve_rows.append((alpha, theta_ab, a1, b1, S))

        if kept:
            S_arr = np.array([k[3] for k in kept])
            for idx in signal.argrelmin(S_arr, order=argrelmin_order)[0]:
                init_rows.append((kept[idx][1], kept[idx][2], alpha, kept[idx][3],
                                  'local_min', 'NotYet'))
            init_rows.append((kept[0][1], kept[0][2], alpha, kept[0][3],
                              'b_contact', 'NotYet'))
            init_rows.append((kept[-1][1], kept[-1][2], alpha, kept[-1][3],
                              'a_contact', 'NotYet'))

        if progress_callback:
            progress_callback((i_alpha + 1) / len(alphas))

    df_curves = pd.DataFrame(curve_rows, columns=['alpha', 'theta_ab', 'a', 'b', 'S'])
    df_init = pd.DataFrame(init_rows, columns=['a', 'b', 'theta', 'S', 'kind', 'status'])
    return df_curves, df_init
