"""py3Dmol-based 3D molecule/cluster viewing, for embedding in Streamlit.

Ported from auto_opt (`app.py`'s `_render_3d` / `plot/make_cluster_xyz.py`'s
`_xyz_block`), split out as reusable, streamlit-agnostic helpers: this module
only builds XYZ strings and py3Dmol HTML; callers embed the HTML themselves
(e.g. via `st.components.v1.html(html, height=...)`).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np


def xyz_block(symbols: Sequence[str], coords) -> str:
    """Format (symbol, xyz) pairs as XYZ-file atom lines (no header)."""
    return "".join(
        f"{s:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n"
        for s, (x, y, z) in zip(symbols, coords)
    )


def to_xyz_string(symbols: Sequence[str], coords, comment: str = "") -> str:
    """Format (symbol, xyz) pairs as a complete XYZ file string."""
    coords = np.asarray(coords)
    return f"{len(symbols)}\n{comment}\n" + xyz_block(symbols, coords)


def render_3d_html(
    xyz_str: str,
    *,
    style: str = "sticks",
    width: int = 500,
    height: int = 400,
) -> str:
    """Render an XYZ string to a standalone py3Dmol HTML snippet.

    `style`: "sticks" (capped sticks + small spheres) or "spacefill".
    """
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    view.addModel(xyz_str, "xyz")
    if style == "spacefill":
        view.setStyle({"sphere": {"scale": 1.0}})
    else:
        view.setStyle({"stick": {"radius": 0.15}, "sphere": {"radius": 0.3}})
    view.setProjection("orthographic")
    view.zoomTo()
    return view._make_html()
