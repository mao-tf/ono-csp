"""Molecule loading and alignment.

`vdw_radius`, `R2atom`, `Rod`, `read_xyz` are ported from auto_opt
(`auto_opt/utils.py`); they are shared with Ono's legacy utils.py, which
uses the same Bondi radii and Rodrigues rotation.

Rigid-body placement into the herringbone layer lives in
`structure.intralayer` (Ono's A2/A3 convention) — this module only loads a
monomer XYZ and aligns it to a reproducible principal-axis frame.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np

_MOLECULE_DIR = Path(__file__).resolve().parents[3] / "data" / "molecules"

# VdW radii (Å) — Bondi values (ported from auto_opt/utils.py)
_VDW: dict[str, float] = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'F': 1.47, 'P': 1.80, 'S': 1.80,
    'CL': 1.75, 'BR': 1.85, 'I': 1.98,
}
_RADIUS_TO_ATOM: dict[float, str] = {v: k for k, v in _VDW.items()}


def vdw_radius(sym: str) -> float:
    """Element symbol -> VdW radius (Å). Unknown elements fall back to carbon."""
    return _VDW.get(sym.strip().upper(), _VDW['C'])


def R2atom(R: float) -> str:
    """VdW radius value -> element symbol. Unknown values fall back to 'C'."""
    return _RADIUS_TO_ATOM.get(round(R, 2), 'C')


def Rod(n: np.ndarray, theta_deg: float) -> np.ndarray:
    """Rodrigues rotation matrix: rotate by theta_deg (degrees) about axis n."""
    nx, ny, nz = n
    c = np.cos(np.radians(theta_deg))
    s = np.sin(np.radians(theta_deg))
    return np.array([
        [c + nx*nx*(1-c),     nx*ny*(1-c) - nz*s,  nx*nz*(1-c) + ny*s],
        [nx*ny*(1-c) + nz*s,  c + ny*ny*(1-c),     ny*nz*(1-c) - nx*s],
        [nx*nz*(1-c) - ny*s,  ny*nz*(1-c) + nx*s,  c + nz*nz*(1-c)   ],
    ])


def read_xyz(path: str | Path) -> List[List]:
    """Read an XYZ file and return [[x, y, z, symbol], ...]."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if len(s) == 4:
                try:
                    x, y, z = float(s[1]), float(s[2]), float(s[3])
                except ValueError:
                    continue
                rows.append([x, y, z, s[0]])
    if not rows:
        raise ValueError(f"No XYZ atom lines found in {path} (expected 'El x y z')")
    return rows


@dataclass
class Molecule:
    """A single monomer in its own principal-axis frame.

    Frame convention (see `load_molecule`): centroid at origin, long axis
    along x, second in-plane axis along y, plane normal along z.
    """
    name: str
    symbols: List[str]
    coords: np.ndarray  # (N, 3) Å, principal-axis frame
    radii: np.ndarray   # (N,) Å, per-atom VdW radius

    @property
    def n_atoms(self) -> int:
        return len(self.symbols)


def load_molecule(name_or_path: str | Path, molecule_dir: str | Path | None = None) -> Molecule:
    """Load a monomer XYZ file and align it to the standard principal-axis frame.

    `name_or_path` is either a preset name (looked up as
    `{molecule_dir or data/molecules}/{name}.xyz`) or a direct path to an
    XYZ file.

    Alignment: subtract the centroid, then PCA (SVD) the coordinates and
    reorder axes by decreasing variance, so the long molecular axis is x,
    the in-plane short axis is y, and the (near-zero, for a planar
    aromatic) out-of-plane axis is z. This matches the frame `place_monomer`
    assumes (see its docstring).
    """
    p = Path(name_or_path)
    if p.suffix.lower() == ".xyz" and p.exists():
        path = p
        name = p.stem
    else:
        name = str(name_or_path)
        _dir = Path(molecule_dir) if molecule_dir else _MOLECULE_DIR
        path = _dir / f"{name}.xyz"

    rows = read_xyz(path)
    symbols = [r[3] for r in rows]
    coords = np.array([[r[0], r[1], r[2]] for r in rows], dtype=float)

    coords = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords, full_matrices=False)
    coords = coords @ vt.T

    radii = np.array([vdw_radius(s) for s in symbols], dtype=float)
    return Molecule(name=name, symbols=symbols, coords=coords, radii=radii)
