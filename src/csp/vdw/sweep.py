#!/usr/bin/env python3
"""
Step 1a: vdW contact sweep for the R-form intralayer structure.

For each herringbone angle alpha, sweeps the T-contact direction theta_c (0-90°)
and records the minimum contact distance R_clps. Outputs a CSV that
extract_init.py uses to find candidate (a, b) pairs.

Output CSV columns: alpha, theta_c, R_clps, TorF
  - R_clps: minimum R such that the T-shaped dimer at direction theta_c is just touching
  - TorF:   True if SP contacts (a-dir, b-dir) also satisfy vdW constraints

Usage:
  python -m csp.vdw.sweep \\
      --monomer data/molecules/anthracene.xyz \\
      --out-dir runs/anthracene/ \\
      --alpha-min 5 --alpha-max 45 --alpha-step 5 \\
      --theta-step 1
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from csp.utils import Rod, read_xyz, vdw_radius


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _rotate(mol: list[list], Rz: np.ndarray) -> list[list]:
    """Rotate molecule (list of [x,y,z,el]) by matrix Rz."""
    result = []
    for x, y, z, el in mol:
        r = Rz @ np.array([x, y, z])
        result.append([r[0], r[1], r[2], el])
    return result


def _t_shaped_pair(mol_rot: list[list]) -> tuple[list[list], list[list]]:
    """Create T-shaped (herringbone) dimer pair from a rotated monomer.

    mol1 = mol_rot
    mol2 = glide-reflected mol_rot: x -> -x  (z-shift is 0 for R-form)
    vdw_R then shifts mol2 along theta_c direction to find contact distance.
    """
    mol1 = mol_rot
    mol2 = [[-x, y, z, el] for x, y, z, el in mol_rot]
    return mol1, mol2


def _sp_pair(mol_rot: list[list]) -> tuple[list[list], list[list]]:
    """Slipped-parallel dimer: mol2 starts at same position as mol1.
    vdw_R then shifts mol2 along theta_deg to find SP contact distance.
    """
    return mol_rot, list(mol_rot)


def _vdw_R(mol1: list[list], mol2: list[list], theta_deg: float) -> float:
    """Minimum shift along (cos(theta), sin(theta), 0) for mol2 to just touch mol1.

    Starting from the current relative position of mol2 vs mol1,
    returns R such that placing mol2 at mol2_pos + R*(cos,sin,0) avoids all vdW clashes.
    """
    R1 = np.array([[a[0], a[1], a[2]] for a in mol1], float)
    R2 = np.array([[a[0], a[1], a[2]] for a in mol2], float)
    r1 = np.array([vdw_radius(a[3]) for a in mol1], float)
    r2 = np.array([vdw_radius(a[3]) for a in mol2], float)

    ct = math.cos(math.radians(theta_deg))
    st = math.sin(math.radians(theta_deg))
    eR = np.array([ct, st, 0.0])

    # D[i,j] = R2[j] - R1[i]
    D = R2[None, :, :] - R1[:, None, :]          # (n1, n2, 3)
    R12b = D @ eR                                  # projection along eR
    D2 = (D * D).sum(axis=2)                       # |D|^2
    R12a2 = D2 - R12b ** 2                         # perpendicular component^2

    rad_sum = r1[:, None] + r2[None, :]            # (n1, n2)
    disc = rad_sum ** 2 - R12a2
    disc = np.maximum(disc, 0.0)
    shift_needed = -R12b + np.sqrt(disc)
    shift_needed = np.maximum(shift_needed, 0.0)

    return float(np.max(shift_needed))


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def sweep(
    monomer_path: str | Path,
    out_dir: str | Path,
    alpha_min: float = 5.0,
    alpha_max: float = 45.0,
    alpha_step: float = 5.0,
    alpha_list: list[float] | None = None,
    theta_step: float = 1.0,
    eps: float = 1e-3,
) -> Path:
    """Sweep T-contact direction theta_c for each herringbone angle alpha.

    For each (alpha, theta_c):
      1. Build T-shaped pair (glide-reflected)
      2. Compute R_clps = minimum contact shift along theta_c
      3. a = 2*R_clps*cos(theta_c), b = 2*R_clps*sin(theta_c)
      4. Check SP contacts: a >= R_SP_a and b >= R_SP_b  (TorF)

    Returns path to the output CSV.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    monomer_name = Path(monomer_path).stem

    if alpha_list:
        alphas = sorted(set(alpha_list))
    else:
        alphas = [round(a, 6) for a in np.arange(alpha_min, alpha_max + 1e-9, alpha_step)]

    theta_cs = [round(t, 6) for t in np.arange(0.0, 90.0 + 1e-9, theta_step)]

    mol0 = read_xyz(monomer_path)
    ez = np.array([0.0, 0.0, 1.0])

    print(f"Molecule : {monomer_name} ({len(mol0)} atoms)")
    print(f"Alpha    : {alphas[0]}° – {alphas[-1]}°  ({len(alphas)} pts)")
    print(f"Theta_c  : 0° – 90°  (step {theta_step}°, {len(theta_cs)} pts)")
    print(f"Total    : {len(alphas) * len(theta_cs)} points")

    rows = []
    for alpha in alphas:
        Rz = Rod(ez, alpha)
        mol_rot = _rotate(mol0, Rz)

        # SP contact distances (molecule in a-direction and b-direction)
        sp1, sp2 = _sp_pair(mol_rot)
        R_sp_a = _vdw_R(sp1, sp2, 0.0)    # shift along x (a-direction)
        R_sp_b = _vdw_R(sp1, sp2, 90.0)   # shift along y (b-direction)

        # T-shaped pair
        t1, t2 = _t_shaped_pair(mol_rot)

        for theta_c in theta_cs:
            R_clps = _vdw_R(t1, t2, theta_c)
            a = 2.0 * R_clps * math.cos(math.radians(theta_c))
            b = 2.0 * R_clps * math.sin(math.radians(theta_c))
            ok = (a >= R_sp_a - eps) and (b >= R_sp_b - eps)
            rows.append([alpha, theta_c, R_clps, ok])

    df = pd.DataFrame(rows, columns=['alpha', 'theta_c', 'R_clps', 'TorF'])
    df = df.sort_values(['alpha', 'theta_c']).reset_index(drop=True)

    out_csv = out / f"vdW_r_contact_{monomer_name}.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}  (n={len(df)})")
    return out_csv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Step 1a: vdW contact sweep (R-form intralayer)"
    )
    ap.add_argument('--monomer', required=True, help='Monomer XYZ or CSV path')
    ap.add_argument('--out-dir', required=True, help='Output directory')

    ap.add_argument('--alpha-min',  type=float, default=5.0)
    ap.add_argument('--alpha-max',  type=float, default=45.0)
    ap.add_argument('--alpha-step', type=float, default=5.0)
    ap.add_argument('--alpha-list', type=float, nargs='+',
                    help='Explicit alpha values (overrides min/max/step)')
    ap.add_argument('--theta-step', type=float, default=1.0,
                    help='theta_c step in degrees (default: 1)')
    ap.add_argument('--eps', type=float, default=1e-3,
                    help='Tolerance for SP contact check (Å)')

    args = ap.parse_args()
    sweep(
        monomer_path=args.monomer,
        out_dir=args.out_dir,
        alpha_min=args.alpha_min,
        alpha_max=args.alpha_max,
        alpha_step=args.alpha_step,
        alpha_list=args.alpha_list,
        theta_step=args.theta_step,
        eps=args.eps,
    )


if __name__ == '__main__':
    main()
