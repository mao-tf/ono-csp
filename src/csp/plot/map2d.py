"""Plotly interactive 2D energy/map heatmaps, for embedding in Streamlit.

Adapted from auto_opt (`app.py`'s `_render_heatmap`), split so this module
only builds the `go.Figure` — callers wire up `st.plotly_chart(...,
on_select="rerun")` themselves and handle the returned selection event, since
that interaction loop is app-specific.
"""
from __future__ import annotations

from typing import Optional, Sequence

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def build_heatmap_figure(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    val_col: str,
    val_label: str,
    *,
    markers: Optional[Sequence[dict]] = None,
) -> go.Figure:
    """Build a 2D heatmap of `val_col` over (`x_col`, `y_col`).

    `markers`: optional overlay points, each a dict with keys
    `x`, `y`, and optionally `label`, `color` (default "white"),
    `symbol` (default "circle-open"), `size` (default 16) — e.g. to mark
    R-form/G-form/N-form points on a Step 2 inclination map (spec.md Tab 4).
    """
    pivot = df.pivot_table(values=val_col, index=y_col, columns=x_col, aggfunc="min")
    if pivot.shape[0] >= 2 and pivot.shape[1] >= 2:
        fig = px.imshow(
            pivot,
            color_continuous_scale="RdBu_r",
            labels={"color": val_label},
            aspect="auto",
        )
    else:
        fig = go.Figure()

    for m in (markers or []):
        fig.add_trace(go.Scatter(
            x=[m["x"]], y=[m["y"]],
            mode="markers",
            marker=dict(
                symbol=m.get("symbol", "circle-open"),
                size=m.get("size", 16),
                color=m.get("color", "white"),
                line=dict(width=2, color=m.get("color", "white")),
            ),
            name=m.get("label", ""),
            showlegend=bool(m.get("label")),
            hoverinfo="name" if m.get("label") else "skip",
        ))

    fig.update_layout(margin=dict(l=20, r=20, t=30, b=20), clickmode="event+select")
    return fig
