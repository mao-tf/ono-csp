#!/usr/bin/env python3
"""csp — Crystal Structure Prediction UI (Streamlit)

Usage:
  streamlit run app.py

Reproduces the structure-search pipeline of Ono et al. (JACS, submitted).
See spec.md for the full method-to-tab mapping and for how each tab maps to
Ono's legacy scripts (legacy/ono_scripts/).

Execution policy (spec.md "実行方式", meeting 2026-07-04):
- The Step 1a vdW scan runs directly in this GUI (it completes on a laptop).
- The DFT steps run from the command line with the legacy scripts — each tab
  shows the commands and displays results instead.
- Every results section accepts a drag & dropped CSV; until one is dropped it
  falls back to the precomputed sample results bundled in example/pentacene/.

Tab layout (spec.md "大野コード対応表", 2026-07-05): Ono's actual pipeline is
3 steps (step1 -> step2 para/twist -> step3 para/twist -> tcal), not the
5-section paper draft this spec.md started from. Tabs are grouped by step
number, with para/twist as sub-tabs within Step 2 and Step 3 (they are
alternative paths, not sequential sub-stages) — Step 1's vdW pre-scan and
DFT refinement are sequential instead, so they're plain sections rather than
sub-tabs.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from csp.structure.molecule import load_molecule, _VDW  # noqa: E402
from csp.structure.intralayer import cluster9, monomer_csv  # noqa: E402
from csp.vdw.contact import step1a_scan  # noqa: E402
from csp.plot.viewer3d import to_xyz_string, render_3d_html  # noqa: E402
from csp.plot.map2d import build_heatmap_figure  # noqa: E402

ROOT = Path(__file__).resolve().parent
MOLECULE_DIR = ROOT / "data" / "molecules"
EXAMPLE_DIR = ROOT / "example" / "pentacene"
PRESET_MOLECULES = ["naphthalene", "anthracene", "tetracene", "pentacene", "hexacene"]

# Route line used by Ono's legacy input generators (make_step*.py)
GAUSSIAN_ROUTE = "#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2"
LEGACY_DIR = "legacy/ono_scripts/stepwise_optimization"

st.set_page_config(page_title="csp — Crystal Structure Prediction", layout="wide")
st.title("csp — Crystal Structure Prediction")
st.caption(
    "Reproduces Ono et al., \"Origin of Layered Herringbone Packing and "
    "Polymorphism in Polyacenes\" (JACS, submitted). See spec.md for the "
    "full method-to-tab mapping."
)


# ══════════════════════════════════════════════════════════
#  Shared helpers
# ══════════════════════════════════════════════════════════
def load_results_csv(key: str, sample_file: str) -> tuple[Optional[pd.DataFrame], str]:
    """Results data source: user-dropped CSV first, bundled sample otherwise.

    Returns (df, source) where source is "uploaded", "sample", or "" when
    neither is available. Implements the sample-until-dropped display policy
    from spec.md "実行方式".
    """
    uploaded = st.file_uploader(
        "Drop your results CSV here (a precomputed sample is shown until then)",
        type="csv", key=key,
    )
    if uploaded is not None:
        try:
            return pd.read_csv(uploaded), "uploaded"
        except Exception as e:
            st.error(f"Could not read the dropped CSV: {e}")
            return None, ""
    sample_path = EXAMPLE_DIR / sample_file
    if sample_path.exists():
        try:
            return pd.read_csv(sample_path), "sample"
        except Exception as e:
            st.error(f"Could not read bundled sample {sample_path.name}: {e}")
            return None, ""
    st.caption(
        f"Sample results are not bundled yet (expected at "
        f"`example/pentacene/{sample_file}`); drop your own CSV above to plot it."
    )
    return None, ""


def source_badge(source: str) -> None:
    if source == "uploaded":
        st.caption("Showing **your uploaded results**.")
    elif source == "sample":
        st.caption("Showing **precomputed sample results** (example/pentacene/). Drop a CSV above to replace.")
    elif source == "scan":
        st.caption("Showing **the scan you just ran in this GUI**.")


def _numeric_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]


def _default_index(cols: Sequence[str], candidates: Sequence[str], fallback: int = 0) -> int:
    """Index of the first column matching a candidate name (case-insensitive)."""
    lowered = [c.lower() for c in cols]
    for cand in candidates:
        if cand in lowered:
            return lowered.index(cand)
    return min(fallback, len(cols) - 1)


def line_plot_section(
    df: pd.DataFrame,
    key_prefix: str,
    x_candidates: Sequence[str],
    y_candidates: Sequence[str],
    agg_min: bool = False,
) -> None:
    """x-vs-y line plot with column pickers.

    Axes default to the first match from `x_candidates`/`y_candidates` but
    stay user-selectable. With `agg_min`, offers taking the minimum of y over
    all other parameters at each x (the optimization envelope, e.g. best E
    per theta from a hill-climb history CSV).
    """
    numeric = _numeric_columns(df)
    if len(numeric) < 2:
        st.warning("Need at least two numeric columns to plot.")
        st.dataframe(df, width="stretch")
        return
    c1, c2 = st.columns(2)
    x_col = c1.selectbox(
        "X axis", numeric, index=_default_index(numeric, x_candidates, 0),
        key=f"{key_prefix}_x",
    )
    y_col = c2.selectbox(
        "Y axis", numeric, index=_default_index(numeric, y_candidates, 1),
        key=f"{key_prefix}_y",
    )
    df_plot = df
    if agg_min:
        if st.checkbox(
            "Minimum over all other parameters at each X (optimization envelope)",
            value=True, key=f"{key_prefix}_agg",
        ):
            df_plot = df.groupby(x_col, as_index=False)[y_col].min()
    fig = px.line(df_plot.sort_values(x_col), x=x_col, y=y_col, markers=True)
    fig.update_layout(margin=dict(l=20, r=20, t=30, b=20))
    st.plotly_chart(fig, width="stretch")
    with st.expander("Data table"):
        st.dataframe(df, width="stretch")


def heatmap_section(
    df: pd.DataFrame,
    key_prefix: str,
    x_candidates: Sequence[str],
    y_candidates: Sequence[str],
    val_candidates: Sequence[str],
    val_label: str,
) -> None:
    """2D heatmap with column pickers (same column-name caveat as line_plot_section)."""
    numeric = _numeric_columns(df)
    if len(numeric) < 3:
        st.warning("Need at least three numeric columns (x, y, value) for a map.")
        st.dataframe(df, width="stretch")
        return
    c1, c2, c3 = st.columns(3)
    x_col = c1.selectbox(
        "X axis", numeric, index=_default_index(numeric, x_candidates, 0),
        key=f"{key_prefix}_x",
    )
    y_col = c2.selectbox(
        "Y axis", numeric, index=_default_index(numeric, y_candidates, 1),
        key=f"{key_prefix}_y",
    )
    val_col = c3.selectbox(
        "Value", numeric, index=_default_index(numeric, val_candidates, 2),
        key=f"{key_prefix}_v",
    )
    fig = build_heatmap_figure(df, x_col, y_col, val_col, val_label)
    st.plotly_chart(fig, width="stretch")
    with st.expander("Data table"):
        st.dataframe(df, width="stretch")


def render_molecule_3d(symbols, coords, comment: str, key_suffix: str, style_key: str) -> None:
    style = st.radio(
        "Style", ["Capped sticks", "Space fill"], horizontal=True, key=style_key
    )
    xyz_str = to_xyz_string(symbols, coords, comment=comment)
    html = render_3d_html(
        xyz_str,
        style="spacefill" if style == "Space fill" else "sticks",
        width=480, height=420,
    )
    st.iframe(html, height=440)
    st.download_button(
        "Download XYZ", data=xyz_str,
        file_name=f"{comment.split()[0]}_{key_suffix}.xyz", mime="text/plain",
        key=f"dl_xyz_{key_suffix}",
    )


def cli_howto(
    *,
    what: str,
    prepare: str,
    setup: str,
    command: str,
    output: str,
    show_route: bool = True,
) -> None:
    """Standard 4-part "how to run this step with Ono's package" block.

    Every CLI section (Step 1 DFT, Step 2/3 para/twist, Transfer Integrals)
    uses the same shape so a user who has read one knows where to look in
    the others: what the step computes, what input file(s) to prepare by
    hand (and where their values come from), what one-time environment
    setup is needed, then the command and the resulting output file.
    """
    st.info(what)
    st.markdown("**1. Prepare inputs**")
    st.markdown(prepare)
    st.markdown("**2. One-time setup**")
    st.markdown(setup)
    st.markdown("**3. Run**")
    st.code(command, language="bash")
    if show_route:
        st.code(GAUSSIAN_ROUTE, language="text")
    st.markdown("**4. Output**")
    st.markdown(output)


_SETUP_MONOMER_ENV = (
    "- Set `CSP_MONOMER_DIR` to the folder holding the monomer CSV "
    "(X, Y, Z, R columns — download it from **Tab 1**), e.g. "
    "`export CSP_MONOMER_DIR=/path/to/data/monomer` "
    "(defaults to `~/path/to/monomer/` if unset)."
)
_SETUP_SCHEDULER = (
    "- The legacy scripts submit jobs with `pjsub` (Fugaku/PJM). On our SGE "
    "cluster, replace the `pjsub`/batch-script generation with `qsub` "
    "(see `src/csp/dft/job_cluster.py` for an SGE job-script template)."
)


# ══════════════════════════════════════════════════════════
#  Tabs
# ══════════════════════════════════════════════════════════
tab_setup, tab_step1, tab_step2, tab_step3, tab_transfer = st.tabs([
    "1. Molecule Setup",
    "2. Step 1 – Intralayer",
    "3. Step 2 – Long-Axis Shift",
    "4. Step 3 – Interlayer Stacking",
    "5. Transfer Integrals",
])


# ══════════════════════════════════════════════════════════
#  Tab 1: Molecule Setup
# ══════════════════════════════════════════════════════════
with tab_setup:
    col_pick, col_view = st.columns([1, 1])

    with col_pick:
        st.subheader("Molecule")
        available_presets = [
            m for m in PRESET_MOLECULES if (MOLECULE_DIR / f"{m}.xyz").exists()
        ]
        source = st.radio("Source", ["Preset", "Upload custom XYZ"], horizontal=True)

        molecule = None
        if source == "Preset":
            if not available_presets:
                st.warning(f"No preset XYZ files found in {MOLECULE_DIR}")
            else:
                preset_name = st.selectbox("Preset molecule", available_presets)
                molecule = load_molecule(preset_name, molecule_dir=MOLECULE_DIR)
        else:
            uploaded = st.file_uploader("Upload a monomer XYZ file", type=["xyz"])
            if uploaded is not None:
                with tempfile.NamedTemporaryFile(
                    suffix=".xyz", mode="w", delete=False, encoding="utf-8"
                ) as tf:
                    tf.write(uploaded.getvalue().decode("utf-8"))
                    tmp_path = tf.name
                molecule = load_molecule(tmp_path)
                molecule.name = Path(uploaded.name).stem

        if molecule is not None:
            st.session_state["molecule"] = molecule
            st.caption(f"**{molecule.name}** — {molecule.n_atoms} atoms")
            st.download_button(
                "Download monomer CSV (X,Y,Z,R — input for the legacy CLI scripts)",
                data=monomer_csv(molecule),
                file_name=f"{molecule.name}.csv", mime="text/csv",
            )

        st.divider()
        st.subheader("VdW radius table")
        st.caption("Editable; used by the vdW scan in Tab 2.")
        default_table = pd.DataFrame(
            {"atom": list(_VDW.keys()), "vdW radius (Å)": list(_VDW.values())}
        )
        vdw_table = st.data_editor(
            default_table, num_rows="dynamic", key="vdw_radius_table", hide_index=True
        )
        st.session_state["vdw_radii_overrides"] = dict(
            zip(vdw_table["atom"], vdw_table["vdW radius (Å)"])
        )

    with col_view:
        st.subheader("3D viewer")
        if molecule is None:
            st.caption("Select or upload a molecule to preview it here.")
        else:
            render_molecule_3d(
                molecule.symbols, molecule.coords,
                f"{molecule.name} (principal-axis frame)",
                key_suffix="setup", style_key="setup_style",
            )


# ══════════════════════════════════════════════════════════
#  Tab 2: Step 1 – Intralayer (vdW pre-scan, then DFT refinement)
# ══════════════════════════════════════════════════════════
with tab_step1:
    mol2 = st.session_state.get("molecule")

    # ─── Section: vdW pre-scan (runs in this GUI) ─────────────────
    st.subheader("vdW pre-scan — runs in this GUI")
    st.caption(
        "Rigid-sphere vdW contact model (same algorithm as legacy "
        "step1.py --init): for each herringbone half-angle alpha, sweep the "
        "T-contact direction and find the smallest cell S = a×b compatible "
        "with T and SP contacts. Initial DFT candidates are the local minima "
        "of S plus both endpoints."
    )
    c1, c2, c3 = st.columns(3)
    alpha_min = c1.number_input("alpha min (deg)", value=5.0, key="s1vdw_amin")
    alpha_max = c2.number_input("alpha max (deg)", value=45.0, key="s1vdw_amax")
    alpha_step = c3.number_input("alpha step (deg)", value=5.0, min_value=0.5, key="s1vdw_astep")

    if mol2 is None:
        st.info("Select a molecule in Tab 1 first.")
    elif st.button(f"Run vdW Scan ({mol2.name})", key="s1vdw_run"):
        alphas = list(np.arange(alpha_min, alpha_max + 1e-9, alpha_step))
        prog = st.progress(0.0, text="Scanning...")
        df_curves, df_init = step1a_scan(
            mol2, alphas,
            radii_overrides=st.session_state.get("vdw_radii_overrides"),
            progress_callback=lambda f: prog.progress(f),
        )
        prog.empty()
        st.session_state["s1vdw_curves"] = df_curves
        st.session_state["s1vdw_init"] = df_init
        st.session_state["s1vdw_scan_mol"] = mol2.name
        # Cache the S-axis range once per scan so it's byte-identical on every
        # rerun (recomputing it per-render caused the y-axis floor to flicker
        # between 0 and 40 on each click — see spec.md).
        s_max = float(df_curves["S"].max())
        y_lo = 40.0 if s_max > 40.0 else 0.0
        st.session_state["s1vdw_yrange"] = [y_lo, s_max * 1.03]

    st.divider()
    st.subheader("vdW results — alpha vs S = a×b")

    df_curves = st.session_state.get("s1vdw_curves")
    df_init = st.session_state.get("s1vdw_init")
    src2 = ""
    if df_init is not None:
        src2 = "scan"
        if st.button("Clear scan results", key="s1vdw_clear"):
            for k in ("s1vdw_curves", "s1vdw_init", "s1vdw_scan_mol", "s1vdw_yrange"):
                st.session_state.pop(k, None)
            st.rerun()
    else:
        df_up, src2 = load_results_csv("s1vdw_results", "step1_init_params.csv")
        if df_up is not None:
            df_up = df_up.rename(columns={c: c.lower() for c in df_up.columns})
            if {"a", "b", "theta"} <= set(df_up.columns):
                df_init = df_up.rename(columns={"s": "S"})
                if "S" not in df_init.columns:
                    df_init["S"] = df_init["a"] * df_init["b"]
            else:
                st.warning("Expected columns a, b, theta (legacy step1_init_params.csv format).")

    if df_init is not None:
        source_badge(src2)

        _KIND_LABEL = {
            "b_contact": "b-stack",
            "a_contact": "a-stack",
            "local_min": "local min",
        }
        _KIND_COLOR = {
            "b_contact": "#1f77b4",
            "a_contact": "#d62728",
            "local_min": "#2ca02c",
        }

        def _default_current(df: pd.DataFrame) -> Optional[dict]:
            if len(df) == 0:
                return None
            row = df.iloc[0]
            return {
                "alpha": float(row["theta"]), "a": float(row["a"]), "b": float(row["b"]),
                "label": f"default: alpha={row['theta']} a={row['a']} b={row['b']}",
            }

        def _render_preview(df: pd.DataFrame, key_suffix: str) -> None:
            st.markdown("**Structure preview**")
            if mol2 is None:
                st.caption("Select a molecule in Tab 1.")
                return
            current = st.session_state.get("s1vdw_current") or _default_current(df)
            if current is None:
                st.caption("No candidates to preview.")
                return
            st.caption(current["label"])
            syms, coords = cluster9(mol2, current["a"], current["b"], current["alpha"])
            render_molecule_3d(
                syms, coords,
                f"{mol2.name} cluster9 alpha={current['alpha']} a={current['a']} b={current['b']}",
                key_suffix=key_suffix, style_key=f"s1vdw_style_{key_suffix}",
            )

        plot_fig2b, plot_figS1c = st.tabs([
            "Fig. 2(b)-style: S vs α", "Fig. S1(c)-style: S vs β (per α)",
        ])

        # ─── Fig. 2(b)-style: S vs alpha (vdW analogue) ───────────────
        with plot_fig2b:
            col_a, col_3d_a = st.columns([2, 1])
            with col_a:
                st.caption(
                    "vdW analogue of Fig. 2(b): candidate cells at each alpha, "
                    "connected within each category (paper plots the "
                    "DFT-optimized E_intra(8) here instead of S; HB/PS/CH "
                    "there ↔ endpoint/local-min category here). Click a point "
                    "to preview its 9-molecule structure."
                )
                figA = go.Figure()
                if df_curves is not None:
                    env = df_curves[df_curves["valid"]].groupby("alpha", as_index=False)["S"].min()
                    figA.add_trace(go.Scatter(
                        x=env["alpha"], y=env["S"], mode="lines",
                        name="min S (feasible)", line=dict(color="lightgray"),
                        hoverinfo="skip",
                    ))
                for kind, grp in (df_init.groupby("kind") if "kind" in df_init.columns else []):
                    grp = grp.sort_values("theta")
                    color = _KIND_COLOR.get(kind, "gray")
                    figA.add_trace(go.Scatter(
                        x=grp["theta"], y=grp["S"], mode="lines+markers",
                        name=_KIND_LABEL.get(kind, kind),
                        line=dict(color=color), marker=dict(size=9, color=color),
                        customdata=grp.index.to_numpy(),
                        hovertemplate="alpha=%{x}<br>S=%{y:.1f} Å²<extra>" + kind + "</extra>",
                    ))
                if "kind" not in df_init.columns:
                    grp = df_init.sort_values("theta")
                    figA.add_trace(go.Scatter(
                        x=grp["theta"], y=grp["S"], mode="lines+markers",
                        name="candidates", marker=dict(size=9, color="tomato"),
                        line=dict(color="tomato"),
                        customdata=grp.index.to_numpy(),
                    ))
                figA.update_layout(
                    xaxis_title="alpha (deg, half of herringbone dihedral)",
                    yaxis_title="S = a×b (Å²)",
                    margin=dict(l=20, r=20, t=30, b=20),
                )
                event_a = st.plotly_chart(
                    figA, width="stretch", on_select="rerun", key="s1vdw_figA"
                )
                pts = event_a.selection.points if (event_a and event_a.selection) else []
                if pts and pts[0].get("customdata") is not None:
                    idx = int(pts[0]["customdata"])
                    if st.session_state.get("s1vdw_figA_prev") != idx:
                        row = df_init.loc[idx]
                        st.session_state["s1vdw_current"] = {
                            "alpha": float(row["theta"]), "a": float(row["a"]), "b": float(row["b"]),
                            "label": f"Fig.2b click: alpha={row['theta']} a={row['a']} b={row['b']}",
                        }
                        st.session_state["s1vdw_figA_prev"] = idx

                st.download_button(
                    "Download step1_init_params.csv (input for the DFT step below)",
                    data=df_init[["a", "b", "theta"]].assign(status="NotYet").to_csv(index=False),
                    file_name="step1_init_params.csv", mime="text/csv",
                    key="dl_s1vdw_init",
                )
                with st.expander("Candidates table"):
                    st.dataframe(df_init, width="stretch")

                labels = [
                    f"alpha={r['theta']}  a={r['a']}  b={r['b']}  S={r['S']:.1f}"
                    + (f"  [{r['kind']}]" if "kind" in df_init.columns else "")
                    for _, r in df_init.iterrows()
                ]
                if labels:
                    sel = st.selectbox(
                        "Or pick a candidate from the list", range(len(labels)),
                        format_func=lambda i: labels[i], key="s1vdw_sel",
                    )
                    if st.session_state.get("s1vdw_sel_prev") != sel:
                        row = df_init.iloc[sel]
                        st.session_state["s1vdw_current"] = {
                            "alpha": float(row["theta"]), "a": float(row["a"]), "b": float(row["b"]),
                            "label": f"list pick: alpha={row['theta']} a={row['a']} b={row['b']}",
                        }
                        st.session_state["s1vdw_sel_prev"] = sel

            with col_3d_a:
                _render_preview(df_init, key_suffix="s1vdw_A")

        # ─── Fig. S1(c)-style: S vs beta (theta_ab) per alpha ─────────
        with plot_figS1c:
            if df_curves is None:
                st.info(
                    "This plot needs the full sweep curves, which are only "
                    "available right after running the scan in this GUI "
                    "(not from an uploaded step1_init_params.csv)."
                )
            else:
                col_b, col_3d_b = st.columns([2, 1])
                with col_b:
                    st.caption(
                        "vdW analogue of Fig. S1(c): S vs contact-direction angle "
                        "β for a few alpha values. **Solid, thin, full color** = "
                        "SP-neighbor vdW spheres clear (valid); **long-dashed, "
                        "thick, faded** = they overlap (infeasible). Click a "
                        "point to preview its structure."
                    )
                    alphas_avail = sorted(df_curves["alpha"].unique())
                    n_default = min(5, len(alphas_avail))
                    default_idx = np.linspace(0, len(alphas_avail) - 1, n_default).round().astype(int)
                    default_alphas = [alphas_avail[i] for i in sorted(set(default_idx))]
                    sel_alphas = st.multiselect(
                        "alpha values to show", alphas_avail, default=default_alphas,
                        key="s1vdw_figB_alphas",
                    )
                    palette = px.colors.qualitative.Plotly
                    figB = go.Figure()
                    for i, alpha in enumerate(sel_alphas):
                        color = palette[i % len(palette)]
                        sub = df_curves[df_curves["alpha"] == alpha].sort_values("theta_ab")
                        # split into contiguous valid/invalid runs for solid/dashed segments
                        group_id = (sub["valid"] != sub["valid"].shift()).cumsum()
                        first = True
                        for _, seg in sub.groupby(group_id):
                            is_valid = bool(seg["valid"].iloc[0])
                            figB.add_trace(go.Scatter(
                                x=seg["theta_ab"], y=seg["S"], mode="lines+markers",
                                line=dict(
                                    color=color,
                                    dash="solid" if is_valid else "longdash",
                                    width=2.5 if is_valid else 4,
                                ),
                                marker=dict(size=5 if is_valid else 3, color=color),
                                opacity=1.0 if is_valid else 0.55,
                                name=f"alpha={alpha}°", legendgroup=f"a{alpha}",
                                showlegend=first,
                                customdata=[[alpha, ta, a, b] for ta, a, b in
                                            zip(seg["theta_ab"], seg["a"], seg["b"])],
                                hovertemplate="beta=%{x}<br>S=%{y:.1f} Å²<extra>alpha="
                                              + str(alpha) + "</extra>",
                            ))
                            first = False

                    y_range = st.session_state.get("s1vdw_yrange")
                    figB.update_layout(
                        xaxis_title="beta (deg, T-contact direction)",
                        yaxis_title="S = a×b (Å²)",
                        margin=dict(l=20, r=20, t=30, b=20),
                    )
                    if y_range:
                        # fixedrange locks the y-axis so the click-triggered
                        # rerun (via on_select) can never leave it at Plotly's
                        # own autorange — without this the floor flickered
                        # between 0 and 40 on alternating clicks.
                        figB.update_yaxes(range=y_range, autorange=False, fixedrange=True)
                    event_b = st.plotly_chart(
                        figB, width="stretch", on_select="rerun", key="s1vdw_figB"
                    )
                    pts_b = event_b.selection.points if (event_b and event_b.selection) else []
                    if pts_b and pts_b[0].get("customdata") is not None:
                        cd = pts_b[0]["customdata"]
                        ident = tuple(cd)
                        if st.session_state.get("s1vdw_figB_prev") != ident:
                            st.session_state["s1vdw_current"] = {
                                "alpha": float(cd[0]), "a": float(cd[2]), "b": float(cd[3]),
                                "label": f"Fig.S1c click: alpha={cd[0]} beta={cd[1]} a={cd[2]} b={cd[3]}",
                            }
                            st.session_state["s1vdw_figB_prev"] = ident

                    st.download_button(
                        "Download full sweep curves CSV",
                        data=df_curves.to_csv(index=False),
                        file_name="step1_vdw_curves.csv", mime="text/csv",
                        key="dl_s1vdw_curves",
                    )

                with col_3d_b:
                    _render_preview(df_init, key_suffix="s1vdw_B")

    st.divider()

    # ─── Section: DFT-D refinement (CLI) ──────────────────────────
    st.subheader("DFT-D refinement — CLI")
    st.subheader("How to run (CLI)")
    cli_howto(
        what=(
            "Step 1 DFT-D refinement runs on an HPC with Gaussian16 using "
            f"Ono's scripts in `{LEGACY_DIR}/` (`step1.py` + `make_step1.py`). "
            "It hill-climbs (a, b) in 0.1 Å steps at each fixed alpha, "
            "minimizing E_intra(8) = 4·E_t + 2·E_p1 + 2·E_p2 (BSSE-corrected "
            "dimer interaction energies)."
        ),
        prepare=(
            "- `step1_init_params.csv` (columns: `a, b, theta, status`) — "
            "generated automatically by `--init` below, using the same vdW "
            "model as the pre-scan above. You don't need to write this by hand."
        ),
        setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
        command=(
            "python step1.py --init --auto-dir /path/to/workdir "
            "--monomer-name pentacene \\\n    --num-nodes 4 --num-init 2\n"
            "# --init (re)builds step1_init_params.csv; omit it on later runs\n"
            "# to keep hill-climbing from where step1.csv left off."
        ),
        output="`<auto-dir>/step1.csv` — columns `a, b, theta, E, E_p1, E_p2, E_t, status, file_name`.",
    )

    st.divider()
    st.subheader("DFT results — alpha vs E_intra(8)  (step1.csv)")
    df, src = load_results_csv("s1dft_results", "step1.csv")
    if df is not None:
        source_badge(src)
        line_plot_section(
            df, "s1dft_plot",
            x_candidates=["theta", "alpha"],
            y_candidates=["e", "e_intra8"],
            agg_min=True,
        )
        st.caption(
            "E = 4·E_t + 2·E_p1 + 2·E_p2 (kcal/mol). The minimum over the "
            "hill-climb history at each alpha is the optimized E_intra(8); "
            "its global minimum identifies the R-form (alpha ≈ 25°)."
        )


# ══════════════════════════════════════════════════════════
#  Tab 3: Step 2 – Long-Axis Shift  (para / twist sub-tabs)
# ══════════════════════════════════════════════════════════
with tab_step2:
    sub_para, sub_twist = st.tabs(["para", "twist (→ Type III)"])

    with sub_para:
        st.subheader("How to run (CLI)")
        cli_howto(
            what=(
                "Step 2 (para variant) scans the molecular shift z along the "
                "long axis at the T-shaped and slipped-parallel contacts, "
                "with the lattice fixed to the Step 1 optimum "
                f"(`{LEGACY_DIR}/step2_para.py` + `make_step2_para.py`; "
                "41 z-points from 0 to 4 Å in one Gaussian input via "
                "--Link1). See spec.md \"タブ別 詳細ガイド\" for how to read "
                "the plot below."
            ),
            prepare=(
                "- `step2_init_params.csv`: **one row**, columns `a, b, "
                "theta` — take the best (a, b) at your chosen alpha from "
                "**Tab 2**'s DFT results (`step1.csv`) and write it "
                "yourself (no generator script for this one)."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command="python step2_para.py --auto-dir /path/to/workdir --monomer-name pentacene",
            output="`<auto-dir>/step2_para.csv` — columns `z, Et, Ep` (mirrored to -4..4 Å).",
        )

        st.divider()
        st.subheader("Results — Et / Ep vs z  (step2_para.csv)")
        df, src = load_results_csv("s2_results", "step2_para.csv")
        if df is not None:
            source_badge(src)
            cols_lower = {c.lower(): c for c in df.columns}
            if {"z", "et", "ep"} <= set(cols_lower):
                zc, etc, epc = cols_lower["z"], cols_lower["et"], cols_lower["ep"]
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=df[zc], y=df[etc], mode="lines+markers", name="E_t (T-shaped)"))
                fig.add_trace(go.Scatter(x=df[zc], y=df[epc], mode="lines+markers", name="E_p (slipped parallel)"))
                fig.add_trace(go.Scatter(
                    x=df[zc], y=4 * df[etc] + 2 * df[epc],
                    mode="lines", name="4·E_t + 2·E_p", line=dict(dash="dash"),
                ))
                fig.update_layout(
                    xaxis_title="z shift along long axis (Å)",
                    yaxis_title="E (kcal/mol)",
                    margin=dict(l=20, r=20, t=30, b=20),
                )
                st.plotly_chart(fig, width="stretch")
                with st.expander("Data table"):
                    st.dataframe(df, width="stretch")
            else:
                line_plot_section(
                    df, "s2_plot",
                    x_candidates=["z"], y_candidates=["et", "ep", "e"],
                )

    with sub_twist:
        st.subheader("How to run (CLI)")
        cli_howto(
            what=(
                "Step 2 (twist variant, → Type III packing with glide "
                "symmetry): introduces the T-contact long-axis shift Rt and "
                "the torsion A2, then re-optimizes (a, b) by hill-climbing "
                f"at each fixed (theta, Rt, A2) (`{LEGACY_DIR}/step2_twist.py`; "
                "E = 4·E_t + 2·E_p1)."
            ),
            prepare=(
                "- `step2_twist_init_params.csv`: rows of `theta, Rt, A2, "
                "a, b, status`. `theta` and starting `a, b` come from Step 1 "
                "(**Tab 2**'s DFT results); `Rt` (twist shift) and `A2` "
                "(torsion) are the new values you want to try (e.g. a grid "
                "of A2 = 0..20° per spec.md's expected "
                "naphthalene/anthracene optimum). Set `status='NotYet'` for "
                "each row."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command=(
                "python step2_twist.py --auto-dir /path/to/workdir "
                "--monomer-name naphthalene \\\n    --num-nodes 4 --num-init 2"
            ),
            output=(
                "`<auto-dir>/step2_twist.csv` — columns `a, b, theta, Rt, A2, "
                "E, E_p1, E_t, status, file_name`."
            ),
        )
        st.divider()
        st.subheader("Results  (step2_twist.csv)")
        df, src = load_results_csv("s4a_results", "step2_twist.csv")
        if df is not None:
            source_badge(src)
            line_plot_section(
                df, "s4a_plot",
                x_candidates=["a2", "rt"],
                y_candidates=["e"],
                agg_min=True,
            )
            st.caption(
                "Expected: E minimum at a twist of ≈ 13° for naphthalene, "
                "≈ 9° for anthracene; no gain for tetracene and longer "
                "(paper §2.4)."
            )


# ══════════════════════════════════════════════════════════
#  Tab 4: Step 3 – Interlayer Stacking  (para / twist sub-tabs)
# ══════════════════════════════════════════════════════════
with tab_step3:
    sub_para3, sub_twist3 = st.tabs(["para", "twist (→ Type III)"])

    with sub_para3:
        st.subheader("How to run (CLI)")
        cli_howto(
            what=(
                "Step 3 (para variant) optimizes the interlayer c-vector "
                "(cx, cy, cz) by 3×3×3 hill-climbing at fixed intralayer "
                f"parameters (a, b, theta, Rt, Rp) (`{LEGACY_DIR}/step3_para.py`). "
                "Each point computes 10 interlayer dimers; E averages the two "
                "stacking patterns and sums the four T-shaped pairs."
            ),
            prepare=(
                "- `step3_para_init_params.csv`: rows of `a, b, theta, Rt, Rp, "
                "cx, cy, cz` (starting points). `a, b, theta` come from Step 1 "
                "(**Tab 2**'s DFT results); `Rt, Rp` (long-axis shifts at the "
                "T-shaped / slipped-parallel contacts) and the initial "
                "`cx, cy, cz` come from the vdW interlayer-distance map "
                "`step3_para_vdw.py` (not yet wired into this GUI — run it "
                "separately to pick starting values)."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command=(
                "python step3_para.py --auto-dir /path/to/workdir "
                "--monomer-name pentacene \\\n    --num-nodes 4 --num-init 2"
            ),
            output=(
                "`<auto-dir>/step3_para.csv` — columns `cx, cy, cz, a, b, theta, "
                "Rt, Rp, E, E_i01, E_ip1, E_ip2, E_it1..4, E_i02, E_ip3, E_ip4, "
                "status, file_name`."
            ),
        )

        st.divider()
        st.subheader("Results — interlayer energy map  (step3_para.csv)")
        df, src = load_results_csv("s3_results", "step3_para.csv")
        if df is not None:
            source_badge(src)
            heatmap_section(
                df, "s3_map",
                x_candidates=["cx", "x"],
                y_candidates=["cy", "y"],
                val_candidates=["e", "e_inter7", "v"],
                val_label="E_inter (kcal/mol)",
            )
            st.caption(
                "Hill-climb output is sparse (visited points only); a dense "
                "V(x,y) map from step3_para_vdw can be dropped here too."
            )

    with sub_twist3:
        st.subheader("How to run (CLI)")
        cli_howto(
            what=(
                "Step 3 (twist variant): interlayer c-vector optimization "
                f"for the twisted (Type III) packing (`{LEGACY_DIR}/step3_twist.py`)."
            ),
            prepare=(
                "- `step3_twist_init_params.csv`: rows of `a, b, theta, Rt, "
                "A2, cx, cy, cz` (starting points). `a, b, theta, Rt, A2` "
                "come from the Step 2 twist result (**Tab 3**, twist "
                "sub-tab) at the energy minimum; initial `cx, cy, cz` are a "
                "guess (0 or small values) that the hill-climb will refine."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command=(
                "python step3_twist.py --auto-dir /path/to/workdir "
                "--monomer-name naphthalene \\\n    --num-nodes 4 --num-init 2"
            ),
            output="`<auto-dir>/step3_twist.csv` — columns `cx, cy, cz, ..., E`.",
        )
        st.divider()
        st.subheader("Results  (step3_twist.csv)")
        df, src = load_results_csv("s4b_results", "step3_twist.csv")
        if df is not None:
            source_badge(src)
            heatmap_section(
                df, "s4b_map",
                x_candidates=["cx", "x"],
                y_candidates=["cy", "y"],
                val_candidates=["e"],
                val_label="E (kcal/mol)",
            )


# ══════════════════════════════════════════════════════════
#  Tab 5: Transfer Integrals  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_transfer:
    st.subheader("How to run (CLI)")
    cli_howto(
        what=(
            "Transfer integrals are computed with the CSV-batch workflow in "
            "`legacy/ono_scripts/tcal_csv/` (`tcal_csv.py`) — a wrapper "
            "around Prof. Matsui's tcal program "
            "(https://github.com/matsui-lab-yamagata/tcal; `tcal_1.py` is "
            "Ono's slightly modified copy). For each arrangement row it "
            "builds T-shaped and slipped-parallel dimer inputs, runs "
            "Gaussian16, and extracts the HOMO-HOMO transfer integral. MOs "
            "are computed at B3LYP/6-31G* (paper METHOD)."
        ),
        prepare=(
            "- `<auto-dir>/init_params.csv`: rows of `a, b, theta, A2, z` — "
            "the representative arrangements (e.g. the R-form optimum from "
            "**Tab 2**'s DFT results, and any Type II/III/IV structures "
            "from **Tab 4**'s para/twist sub-tabs) you want J for."
        ),
        setup=(
            _SETUP_MONOMER_ENV.replace("Tab 1", "Tab 1 — this feeds `get_monomer_xyzR` in tcal_csv.py too") + "\n"
            + "- `job.sh` (next to `tcal_csv.py`) is a Fugaku/PJM batch script "
            "template; edit it (or `qsub_process()`'s `pjsub job.sh` call) "
            "for our SGE cluster."
        ),
        command=(
            "python tcal_csv.py --init   --auto-dir /path/to/workdir --monomer-name pentacene   # build inputs\n"
            "python tcal_csv.py --qsub   --auto-dir /path/to/workdir --monomer-name pentacene   # submit Gaussian\n"
            "python tcal_csv.py --tcal   --auto-dir /path/to/workdir --monomer-name pentacene   # run tcal_1.py\n"
            "python tcal_csv.py --result --auto-dir /path/to/workdir --monomer-name pentacene   # collect result.txt"
        ),
        output="`<auto-dir>/result.txt` — space-separated HOMO transfer integrals per row (T-shaped, then slipped-parallel).",
        show_route=False,
    )

    st.divider()
    st.subheader("Results — J per arrangement / contact type")
    st.caption(
        "Drop a CSV (e.g. result.txt converted to CSV with columns like "
        "a, b, theta, J_t, J_p)."
    )
    df, src = load_results_csv("s5_results", "transfer_integrals.csv")
    if df is not None:
        source_badge(src)
        numeric = _numeric_columns(df)
        non_numeric = [c for c in df.columns if c not in numeric]
        if numeric and non_numeric:
            c1, c2, c3 = st.columns(3)
            cat_col = c1.selectbox("Category (x)", non_numeric, key="s5_cat")
            val_col = c2.selectbox(
                "J value", numeric,
                index=_default_index(numeric, ["j", "j_mev", "j (mev)"], 0),
                key="s5_val",
            )
            color_col = c3.selectbox(
                "Group by (color)", ["(none)"] + non_numeric, key="s5_color",
            )
            fig = px.bar(
                df, x=cat_col, y=val_col,
                color=None if color_col == "(none)" else color_col,
                barmode="group",
            )
            fig.update_layout(margin=dict(l=20, r=20, t=30, b=20))
            st.plotly_chart(fig, width="stretch")
            with st.expander("Data table"):
                st.dataframe(df, width="stretch")
        else:
            line_plot_section(
                df, "s5_plot",
                x_candidates=["theta", "alpha", "z"],
                y_candidates=["j_t", "j_p", "j"],
            )
