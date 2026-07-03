"""Rigid-sphere VdW contact distance and Step 1a coarse a,b minimization.

`vdw_contact_distance` is ported directly from auto_opt (`utils.vdw_R`) — it
is a generic "how far apart do two rigid-sphere atom sets need to be along a
given direction to just touch" routine and does not encode any
molecule/packing-specific geometry, so it is safe to reuse as-is.

`min_ab_for_alpha` (Step 1a: for a given herringbone angle alpha, find the
smallest unit cell a x b satisfying VdW contact for the 8 T-shaped/SP-shaped
neighbors — spec.md §Step 1a) is NOT yet implemented here. auto_opt's
equivalent sweep (`vdw/sweep_phi.py`) assumes a 9-molecule a-stack/b-stack
cluster built for upright (BTBT/DNTT-like) molecules, which is not the same
neighbor geometry as csp's flat-lying glide layer (T-type x4 + SP-type x4
around one central molecule, per spec.md). Rather than guess that geometry,
this is left as a stub pending Ono's legacy scripts (see legacy/ono_scripts/)
or direct confirmation of the neighbor construction against the paper.
"""
from __future__ import annotations

import math
from typing import List

import numpy as np


def vdw_contact_distance(axyz_1: List[List], axyz_2: List[List], theta_deg: float) -> float:
    """Rigid-sphere contact distance between two atom sets along a direction.

    `axyz_1`, `axyz_2`: lists of [x, y, z, symbol] (see `molecule.read_xyz`).
    `theta_deg`: direction (in the xy-plane) along which molecule 2 is pushed
    away from molecule 1.

    Returns the minimum push distance (Å) along that direction such that no
    atom pair overlaps (sum of VdW radii). Returns 0.0 if already clear of
    all contacts.
    """
    from csp.structure.molecule import vdw_radius

    R1 = np.asarray([[x, y, z] for x, y, z, _ in axyz_1], float)
    R2 = np.asarray([[x, y, z] for x, y, z, _ in axyz_2], float)
    r1 = np.asarray([vdw_radius(a[3]) for a in axyz_1], float)
    r2 = np.asarray([vdw_radius(a[3]) for a in axyz_2], float)

    ct = math.cos(math.radians(theta_deg))
    st = math.sin(math.radians(theta_deg))
    eR = np.array([ct, st, 0.0], float)

    D = R2[None, :, :] - R1[:, None, :]
    R12b = D @ eR
    R12a2 = (D * D).sum(axis=2) - R12b ** 2

    rad_sum = r1[:, None] + r2[None, :]
    rad_sum_sq = rad_sum ** 2
    mask = R12a2 < rad_sum_sq

    if not np.any(mask):
        return 0.0

    sq = rad_sum_sq[mask] - R12a2[mask]
    two_r_need = np.maximum(-R12b[mask] + np.sqrt(sq), 0.0)
    return float(np.max(two_r_need))


def min_ab_for_alpha(*args, **kwargs):
    """Step 1a coarse a,b scan for a given alpha — not yet implemented.

    See module docstring: needs the T-type (x4) / SP-type (x4) neighbor
    construction from spec.md §Step 1a, which should come from Ono's
    existing scripts rather than be guessed here.
    """
    raise NotImplementedError(
        "min_ab_for_alpha: pending Step 1a neighbor geometry — "
        "see legacy/ono_scripts/ once Ono's existing code is integrated, "
        "or spec.md §Step 1a for the algorithm description."
    )
