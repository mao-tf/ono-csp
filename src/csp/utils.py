"""Common utilities: geometry, I/O, vdW radii."""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np

# Bondi vdW radii (Å)
_VDW: dict[str, float] = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'F': 1.47, 'P': 1.80, 'S': 1.80,
    'CL': 1.75, 'BR': 1.85, 'I': 1.98,
}


def vdw_radius(sym: str) -> float:
    return _VDW.get(sym.strip().upper(), _VDW['C'])


def Rod(n: np.ndarray, theta_deg: float) -> np.ndarray:
    """Rodrigues rotation matrix: rotate by theta_deg around axis n."""
    nx, ny, nz = n
    c = np.cos(np.radians(theta_deg))
    s = np.sin(np.radians(theta_deg))
    return np.array([
        [c + nx*nx*(1-c),     nx*ny*(1-c) - nz*s,  nx*nz*(1-c) + ny*s],
        [nx*ny*(1-c) + nz*s,  c + ny*ny*(1-c),     ny*nz*(1-c) - nx*s],
        [nx*nz*(1-c) - ny*s,  ny*nz*(1-c) + nx*s,  c + nz*nz*(1-c)   ],
    ])


def read_xyz(path: str | Path) -> list[list]:
    """Read XYZ file. Returns [[x, y, z, element], ...].
    Handles standard XYZ (with or without comment line) and also
    the X,Y,Z,R CSV format used by auto_opt.
    """
    rows = []
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == '.csv':
        import pandas as pd
        df = pd.read_csv(p)
        if {'X', 'Y', 'Z', 'R'}.issubset(df.columns):
            from csp.utils import _R_to_sym
            for _, row in df.iterrows():
                sym = _R_to_sym(float(row['R']))
                rows.append([float(row['X']), float(row['Y']), float(row['Z']), sym])
            return rows
    # Standard XYZ: skip lines that don't parse as "El x y z"
    with open(p) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    rows.append([x, y, z, parts[0]])
                except ValueError:
                    continue
    if not rows:
        raise ValueError(f"No atom lines found in {path}")
    return rows


_R_TO_SYM: dict[float, str] = {round(v, 2): k for k, v in _VDW.items()}


def _R_to_sym(R: float) -> str:
    return _R_TO_SYM.get(round(R, 2), 'C')
