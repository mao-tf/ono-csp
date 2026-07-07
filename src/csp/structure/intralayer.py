"""Glide-symmetric herringbone layer construction.

Ported from Ohno's legacy code
(legacy/ono_scripts/stepwise_optimization/make_step1.py: get_monomer_xyzR,
make_xyzfile, make_gjf_xyz).

Frame and variable conventions (Ohno's, kept as-is so results and CSVs stay
comparable with the legacy scripts):

- Monomer frame: molecular long axis along z; at A3=0 the molecular plane is
  the yz-plane (plane normal along x). NOTE: which in-plane axis the plane
  normal points to at A3=0 is a convention; swapping it maps theta -> 90-theta
  with the a/b labels exchanged. To be verified against Ohno's own monomer
  CSVs when they arrive.
- A2: rotation about -x ("torsion"; 0 throughout Step 1, scanned in Step 2
  twist).
- A3: rotation about z = half of the herringbone dihedral angle. This is the
  `theta` column of the step CSVs (paper's alpha, ~25° at the optimum). The
  T-shaped sublattice uses -A3.
- Intralayer lattice: a along x, b along y. T-shaped neighbors sit at
  (±a/2, ±b/2) with -A3; slipped-parallel (SP) neighbors at (±a, 0) and
  (0, ±b) with +A3.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd

from csp.structure.molecule import Molecule, Rod

_EX = np.array([1.0, 0.0, 0.0])
_EZ = np.array([0.0, 0.0, 1.0])

# PCA frame (long=x, in-plane short=y, normal=z; see molecule.load_molecule)
# -> Ohno layer frame (long=z, in-plane short=y, normal=x).
# Proper rotation (det=+1): (X, Y, Z) -> (Z, Y, -X).
_LAYER_FRAME_ROT = np.array([
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
    [-1.0, 0.0, 0.0],
])


def to_layer_frame(mol: Molecule) -> Molecule:
    """Return a copy of `mol` rotated from the PCA frame into the layer frame."""
    return Molecule(
        name=mol.name,
        symbols=list(mol.symbols),
        coords=mol.coords @ _LAYER_FRAME_ROT.T,
        radii=mol.radii.copy(),
    )


def place_in_layer(
    coords: np.ndarray,
    Ta: float, Tb: float, Tc: float,
    A2: float, A3: float,
) -> np.ndarray:
    """Ohno's get_monomer_xyzR rotation+translation on layer-frame coords."""
    xyz = coords @ Rod(-_EX, A2).T
    xyz = xyz @ Rod(_EZ, A3).T
    return xyz + np.array([Ta, Tb, Tc])


# (Ta, Tb, Tc, A3 sign) for the 9-molecule cluster: center, 4 SP, 4 T.
def _cluster9_configs(a: float, b: float) -> List[Tuple[float, float, float, int]]:
    return [
        (0.0,   0.0,  0.0, +1),
        (0.0,   b,    0.0, +1),
        (0.0,  -b,    0.0, +1),
        (a,     0.0,  0.0, +1),
        (-a,    0.0,  0.0, +1),
        (a/2,   b/2,  0.0, -1),
        (a/2,  -b/2,  0.0, -1),
        (-a/2, -b/2,  0.0, -1),
        (-a/2,  b/2,  0.0, -1),
    ]


def cluster9(
    mol: Molecule, a: float, b: float, theta: float, A2: float = 0.0,
    *, in_layer_frame: bool = False,
) -> Tuple[List[str], np.ndarray]:
    """9-molecule intralayer cluster (Ohno's make_xyzfile), for visualization.

    `theta` is A3 (half of the herringbone dihedral). Returns (symbols,
    coords) covering all 9 molecules in order center, SP x4, T x4.
    """
    mol_l = mol if in_layer_frame else to_layer_frame(mol)
    parts = [
        place_in_layer(mol_l.coords, ta, tb, tc, A2, sign * theta)
        for ta, tb, tc, sign in _cluster9_configs(a, b)
    ]
    symbols = list(mol_l.symbols) * 9
    return symbols, np.vstack(parts)


def dimer(
    mol: Molecule, kind: str, a: float, b: float, theta: float, A2: float = 0.0,
    z: float = 0.0, *, in_layer_frame: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Central molecule + one neighbor, as used by the Step 1/2 energy sums.

    kind: 't' -> T-shaped at (a/2, b/2) with -theta;
          'p1' -> SP along b at (0, b); 'p2' -> SP along a at (a, 0).
    `z` shifts the neighbor along the molecular long axis (Step 2's scan
    variable; 0 reproduces the plain Step 1 dimer).
    E_intra(8) = 4*E_t + 2*E_p1 + 2*E_p2 (legacy step1.py).
    """
    mol_l = mol if in_layer_frame else to_layer_frame(mol)
    c_i = place_in_layer(mol_l.coords, 0, 0, 0, A2, theta)
    if kind == 't':
        c_j = place_in_layer(mol_l.coords, a/2, b/2, z, A2, -theta)
    elif kind == 'p1':
        c_j = place_in_layer(mol_l.coords, 0, b, z, A2, theta)
    elif kind == 'p2':
        c_j = place_in_layer(mol_l.coords, a, 0, z, A2, theta)
    else:
        raise ValueError(f"kind must be 't', 'p1' or 'p2': {kind!r}")
    return c_i, c_j


def cluster6_inclined(
    mol: Molecule, a: float, b: float, theta: float, zt1: float, zp: float,
    A2: float = 0.0, *, in_layer_frame: bool = False,
) -> Tuple[List[str], np.ndarray]:
    """6-neighbor + center cluster under a uniform long-axis inclination
    (paper's Fig. 5a/5c), for the Fig. 5(b)-style (theta_incl, phi_incl) map.

    A uniform inclination is the plane z(x, y) = kx*x + ky*y; substituting
    the 4 T-shaped and 2 slipped-parallel (along b) neighbor positions from
    the 9-molecule cluster gives each neighbor's z purely in terms of two
    free parameters `zt1` (T-neighbor at (a/2, b/2)) and `zp` (SP-neighbor
    at (0, b)) -- verified algebraically to reproduce the plane equation
    exactly (spec.md "N形2次元マップの再構成方法"). zt1=zp=0 is the flat
    R-form; zp=0, zt1!=0 is the glide-symmetric G-form direction
    (equivalent to step2_para.py's single-z scan); zt1 != zp/2 breaks glide
    symmetry (N-form direction).
    """
    mol_l = mol if in_layer_frame else to_layer_frame(mol)
    c0 = place_in_layer(mol_l.coords, 0.0, 0.0, 0.0, A2, theta)
    configs = [
        (a/2,  b/2,  zt1,     -theta),
        (-a/2, b/2,  zp-zt1,  -theta),
        (a/2,  -b/2, zt1-zp,  -theta),
        (-a/2, -b/2, -zt1,    -theta),
        (0.0,  b,    zp,       theta),
        (0.0,  -b,   -zp,      theta),
    ]
    parts = [c0] + [
        place_in_layer(mol_l.coords, ta, tb, tc, A2, a3) for ta, tb, tc, a3 in configs
    ]
    symbols = list(mol_l.symbols) * 7
    return symbols, np.vstack(parts)


def monomer_csv(mol: Molecule) -> str:
    """Layer-frame monomer as Ohno's CSV format (columns X, Y, Z, R).

    This is the file the legacy CLI scripts read from
    '~/path/to/monomer/{name}.csv' (get_monomer_xyzR); offer it as a
    download so users can run the DFT steps without hand-building it.
    """
    mol_l = to_layer_frame(mol)
    df = pd.DataFrame(mol_l.coords, columns=["X", "Y", "Z"])
    df["R"] = mol_l.radii
    return df.to_csv(index=False)
