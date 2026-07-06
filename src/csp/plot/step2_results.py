"""Reconstruct the paper's Fig. 5(b)-style (theta_incl, phi_incl) map from
step2_para.csv's independent Et(z)/Ep(z) 1D scans.

Confirmed against Ono's own `step2_para.py::plot2d()` (added 2026-07-06,
after this was independently derived and verified algebraically -- see
spec.md "N形2次元マップの再構成方法"): a uniform long-axis inclination is
the plane z(x, y) = kx*x + ky*y, so Eintra(6) at any (theta_incl, phi_incl)
is just a combination of the *already-computed* Et(z)/Ep(z) values at the
z implied by that plane for each of the 6 neighbors -- no new DFT needed.
"""
from __future__ import annotations

import cmath
import math

import numpy as np
import pandas as pd


def build_theta_phi_map(df_step2: pd.DataFrame, a: float, b: float) -> pd.DataFrame:
    """(theta_incl, phi_incl) Eintra(6) map from step2_para.csv (columns z,
    Et, Ep). Returns columns zt, zp, theta_incl, phi_incl, x, y, E, Et1, Et2,
    Ep.

    `x`/`y` are Ono's own plot2d() axes -- theta_incl as a radius and
    phi_incl as the angle, in a polar-to-Cartesian transform
    (x = theta_incl*cos(phi_incl), y = theta_incl*sin(phi_incl)), NOT plain
    (phi_incl, theta_incl) axes. This is what actually reproduces the
    paper's Fig. 5(b) layout (x in [-45, 45], y in [-30, 30]): the flat
    R-form sits at the origin, and the radial distance from it is the
    inclination magnitude in whatever direction phi_incl points.
    """
    z_vals = sorted(round(float(z), 1) for z in df_step2["z"])
    et_lookup = dict(zip((round(float(z), 1) for z in df_step2["z"]), df_step2["Et"]))
    ep_lookup = dict(zip((round(float(z), 1) for z in df_step2["z"]), df_step2["Ep"]))
    z_min, z_max = min(z_vals), max(z_vals)

    rows = []
    for zt in z_vals:
        for zp in z_vals:
            zt2 = round(zt - zp, 1)
            if zt2 < z_min or zt2 > z_max:
                continue
            Et1, Et2, Ep = et_lookup.get(zt), et_lookup.get(zt2), ep_lookup.get(zp)
            if Et1 is None or Et2 is None or Ep is None:
                continue
            E = 2 * (Et1 + Et2 + Ep)
            za, zb = 2 * zt - zp, zp
            Z = 1.0 / math.sqrt(1.0 + (za / a) ** 2 + (zb / b) ** 2)
            theta_incl = math.degrees(math.acos(Z))
            phi_incl = math.degrees(cmath.phase(complex(za / a, zb / b)))
            phi_rad = math.radians(phi_incl)
            rows.append({
                "zt": zt, "zp": zp, "theta_incl": theta_incl, "phi_incl": phi_incl,
                "x": theta_incl * math.cos(phi_rad), "y": theta_incl * math.sin(phi_rad),
                "E": E, "Et1": Et1, "Et2": Et2, "Ep": Ep,
            })
    return pd.DataFrame(rows)
