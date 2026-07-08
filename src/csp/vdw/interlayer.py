"""Interlayer vdW pre-scan (Step 3 vdW): V(x,y) map before DFT refinement.

Vectorized, faithful port of Ono's step3_para_vdw.py `get_c_vec_vdw`
(legacy/ono_scripts/stepwise_optimization/step3_para_vdw.py). That function
loops over every atom pair in pure Python for each (Ra, Rb) grid point,
which is far too slow for a GUI-interactive scan at its native 0.1 Å
resolution; this module vectorizes the inner atom-pair loop with numpy
while keeping the exact same physics, so results match to numerical
precision (checked against the original for a small test case).

Physical picture (paper §2.3, Fig. 6a-d upper panels): given a fixed
intralayer arrangement (a, b, theta, Rt, Rp — see spec.md 大野コード対応表 for
how Rt/Rp encode the R/G/N-form direction), slide the upper layer by
(Ra, Rb) in the ab-plane and find the vdW-limited interlayer distance
z(Ra, Rb) that just avoids atomic overlap between the two layers, then
V(Ra, Rb) = a*b*z(Ra, Rb). The paper's stable candidates are the local
minima of z (equivalently V) — the "valley" regions in Fig. 6b-d.
"""
from __future__ import annotations

from typing import Mapping, Optional

import numpy as np
import pandas as pd

from csp.structure.molecule import Molecule
from csp.structure.intralayer import to_layer_frame, place_in_layer


def _effective_radii(mol: Molecule, overrides: Optional[Mapping[str, float]]) -> np.ndarray:
    if not overrides:
        return mol.radii
    table = {str(k).strip().upper(): float(v) for k, v in overrides.items()}
    return np.array(
        [table.get(s.strip().upper(), r) for s, r in zip(mol.symbols, mol.radii)],
        dtype=float,
    )


def _build_layer_cluster(
    coords_p: np.ndarray, coords_t: np.ndarray, radii: np.ndarray,
    a_vec: np.ndarray, b_vec: np.ndarray, pattern: int,
) -> tuple[np.ndarray, np.ndarray]:
    """9-molecule intralayer neighbor cluster, per Ono's arr_list1/arr_list2.

    `pattern` selects which of the two glide-related patterns (1 or 2); the
    'p' slots use `coords_p`, the 't' slots use `coords_t`.
    """
    zero = np.zeros(3)
    t1 = (a_vec + b_vec) / 2
    t2 = (a_vec - b_vec) / 2
    if pattern == 1:
        placements = [(zero, 'p'), (b_vec, 'p'), (-b_vec, 'p'), (a_vec, 'p'), (-a_vec, 'p'),
                      (t1, 't'), (-t1, 't'), (t2, 't'), (-t2, 't')]
    else:
        placements = [(zero, 't'), (b_vec, 't'), (-b_vec, 't'), (a_vec, 't'), (-a_vec, 't'),
                      (t1, 'p'), (-t1, 'p'), (t2, 'p'), (-t2, 'p')]
    coords_out, radii_out = [], []
    for offset, kind in placements:
        base = coords_p if kind == 'p' else coords_t
        coords_out.append(base + offset)
        radii_out.append(radii)
    return np.vstack(coords_out), np.concatenate(radii_out)


def _z_max_for_shift(
    over_coords: np.ndarray, over_radii: np.ndarray,
    under_coords: np.ndarray, under_radii: np.ndarray,
) -> float:
    """Minimum z such that `over_coords` (already xy/z-shifted) clears every
    atom in `under_coords` — vectorized equivalent of Ono's inner double loop.
    """
    d_xy2 = ((over_coords[:, None, :2] - under_coords[None, :, :2]) ** 2).sum(-1)
    rad_sum2 = (over_radii[:, None] + under_radii[None, :]) ** 2
    z_sq = rad_sum2 - d_xy2
    mask = z_sq >= 0
    if not np.any(mask):
        return 0.0
    needed = np.sqrt(np.maximum(z_sq, 0.0)) + (under_coords[None, :, 2] - over_coords[:, None, 2])
    return float(max(0.0, needed[mask].max()))


def interlayer_vdw_scan(
    mol: Molecule,
    a: float, b: float, theta: float, Rt: float, Rp: float,
    A2: float = 0.0,
    radii_overrides: Optional[Mapping[str, float]] = None,
    ra_step: float = 0.1, rb_step: float = 0.1,
) -> pd.DataFrame:
    """Step 3 vdW pre-scan: V(Ra, Rb) map at fixed (a, b, theta, Rt, Rp).

    `A2` is the twist torsion (Step 2 twist / Tab 3 twist's A2, 0 by default
    for the untwisted para case) applied to every molecule in the two
    intralayer clusters, so this same scan also covers the twisted (Type
    III) packing at fixed Rp=0 (glide symmetric, matching step3_twist.py's
    parameter set of a, b, theta, Rt, A2 with no separate Rp).

    Returns a DataFrame with columns Ra, Rb, z, V — one row per grid point,
    Ra in [-a/2, a/2] and Rb in [-b/2, b/2] at the given step (0.1 Å by
    default, matching the legacy script).
    """
    mol_l = to_layer_frame(mol)
    radii = _effective_radii(mol_l, radii_overrides)

    coords_p = place_in_layer(mol_l.coords, 0, 0, 0, A2, theta)
    coords_t = place_in_layer(mol_l.coords, 0, 0, 0, A2, -theta)

    a_vec = np.array([a, 0.0, 2 * Rt - Rp])
    b_vec = np.array([0.0, b, Rp])

    under1, rad1 = _build_layer_cluster(coords_p, coords_t, radii, a_vec, b_vec, pattern=1)
    under2, rad2 = _build_layer_cluster(coords_p, coords_t, radii, a_vec, b_vec, pattern=2)

    # Mirror the legacy grid construction exactly (round the half-width to
    # 0.1 Å first, then snap every grid value to 0.1 Å) so results match
    # step3_para_vdw.py's get_c_vec_vdw to numerical precision.
    a_half = round(a / 2, 1)
    b_half = round(b / 2, 1)
    n_ra = int(round(2 * a_half / ra_step)) + 1
    n_rb = int(round(2 * b_half / rb_step)) + 1
    ra_list = np.round(np.linspace(-a_half, a_half, n_ra), 1)
    rb_list = np.round(np.linspace(-b_half, b_half, n_rb), 1)

    rows = []
    for Ra in ra_list:
        z_shift_a = (2 * Rt - Rp) * Ra / a
        over1_xy = coords_p.copy()
        over1_xy[:, 0] += Ra
        over2_xy = coords_t.copy()
        over2_xy[:, 0] += Ra
        for Rb in rb_list:
            z_shift = z_shift_a + Rp * Rb / b
            over1 = over1_xy.copy()
            over1[:, 1] += Rb
            over1[:, 2] += z_shift
            over2 = over2_xy.copy()
            over2[:, 1] += Rb
            over2[:, 2] += z_shift

            z_max1 = _z_max_for_shift(over1, radii, under1, rad1)
            z_max2 = _z_max_for_shift(over2, radii, under2, rad2)
            z_max = max(z_max1, z_max2)
            # `z_max` is only the *additional* vdW-clearing gap on top of the
            # z_shift baseline already applied to over1/over2 above (matching
            # Ono's original get_c_vec_vdw, which returns this same
            # increment — verified bit-for-bit against it). `cz` is the
            # actual overlayer z-coordinate (z_shift + z_max): use *this* one
            # as the c-vector's z component downstream (3D preview, exported
            # init CSV) — using `z` alone silently drops the z_shift term and
            # places the overlayer too low whenever Rt or Rp is nonzero,
            # causing real atomic overlap.
            rows.append((float(Ra), float(Rb), z_max, a * b * z_max, z_shift + z_max))

    return pd.DataFrame(rows, columns=['Ra', 'Rb', 'z', 'V', 'cz'])


def bilayer_preview(
    mol: Molecule,
    a: float, b: float, theta: float, Rt: float, Rp: float,
    cx: float, cy: float, cz: float,
    A2: float = 0.0,
    radii_overrides: Optional[Mapping[str, float]] = None,
) -> tuple[list, np.ndarray]:
    """9-molecule underlayer cluster (pattern 1) + 1 overlayer molecule at
    (cx, cy, cz), for a 3D preview of one point from `interlayer_vdw_scan`
    (analogous to the paper's Fig. 6e "overlayer vs underlayer" view).

    `A2`: see `interlayer_vdw_scan` -- the twist torsion, 0 for the
    untwisted para case.
    """
    mol_l = to_layer_frame(mol)
    radii = _effective_radii(mol_l, radii_overrides)
    coords_p = place_in_layer(mol_l.coords, 0, 0, 0, A2, theta)
    coords_t = place_in_layer(mol_l.coords, 0, 0, 0, A2, -theta)

    a_vec = np.array([a, 0.0, 2 * Rt - Rp])
    b_vec = np.array([0.0, b, Rp])
    under1, _ = _build_layer_cluster(coords_p, coords_t, radii, a_vec, b_vec, pattern=1)

    over = coords_p + np.array([cx, cy, cz])
    coords = np.vstack([under1, over])
    symbols = list(mol_l.symbols) * 10
    return symbols, coords
