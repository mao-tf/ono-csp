#!/usr/bin/env python3
"""csp — Crystal Structure Prediction UI (Streamlit)

Usage:
  streamlit run app.py

Reproduces the structure-search pipeline of Ono et al. (JACS, submitted).
See spec.md for the full method-to-tab mapping.

Execution policy (spec.md "実行方式", meeting 2026-07-04):
- The Step 1a vdW scan runs directly in this GUI (it completes on a laptop).
- The DFT steps (1b, 2, 3, 4a/4b, 5) are run from the command line, not from
  the GUI — each tab shows how to run them and displays results instead.
- Every results section accepts a drag & dropped CSV; until one is dropped it
  falls back to the precomputed sample results bundled in example/pentacene/.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
import plotly.express as px
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from csp.structure.molecule import load_molecule, _VDW  # noqa: E402
from csp.plot.viewer3d import to_xyz_string, render_3d_html  # noqa: E402
from csp.plot.map2d import build_heatmap_figure  # noqa: E402

ROOT = Path(__file__).resolve().parent
MOLECULE_DIR = ROOT / "data" / "molecules"
EXAMPLE_DIR = ROOT / "example" / "pentacene"
PRESET_MOLECULES = ["naphthalene", "anthracene", "tetracene", "pentacene", "hexacene"]

GAUSSIAN_KEYWORDS = "# B3LYP empiricaldispersion=gd3bj 6-311g** counterpoise=2 nosymm"

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
) -> None:
    """x-vs-y line plot with column pickers.

    Column names in Ono's result files are not fixed yet, so the axes default
    to the first match from `x_candidates`/`y_candidates` but stay
    user-selectable.
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
    fig = px.line(df.sort_values(x_col), x=x_col, y=y_col, markers=True)
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


def cli_intro(what: str, note: str = "") -> None:
    """Standard "this step runs from the command line" introduction block."""
    st.info(
        f"{what} runs from the command line, not from this GUI "
        "(it needs Gaussian16 and typically an HPC cluster). "
        "The packaged CLI commands will be documented here once Ono's "
        "scripts (legacy/ono_scripts/) are integrated."
        + (f" {note}" if note else "")
    )


# ══════════════════════════════════════════════════════════
#  Tabs
# ══════════════════════════════════════════════════════════
tab_setup, tab_step1_vdw, tab_step1_dft, tab_step2, tab_step3, tab_step4, tab_step5 = st.tabs([
    "1. Molecule Setup",
    "2. Step 1 – Intralayer vdW Scan",
    "3. Step 1 – DFT-D Optimization",
    "4. Step 2 – Inclination Map",
    "5. Step 3 – Interlayer Stacking",
    "6. Step 4 – Refinement",
    "7. Transfer Integrals",
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

        st.divider()
        st.subheader("VdW radius table")
        st.caption("Editable; used by the vdW contact scans in later tabs.")
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
            style = st.radio(
                "Style", ["Capped sticks", "Space fill"], horizontal=True, key="setup_style"
            )
            xyz_str = to_xyz_string(
                molecule.symbols, molecule.coords,
                comment=f"{molecule.name} (principal-axis frame)",
            )
            html = render_3d_html(
                xyz_str,
                style="spacefill" if style == "Space fill" else "sticks",
                width=480, height=420,
            )
            st.iframe(html, height=440)
            st.download_button(
                "Download aligned XYZ", data=xyz_str,
                file_name=f"{molecule.name}_aligned.xyz", mime="text/plain",
            )


# ══════════════════════════════════════════════════════════
#  Tab 2: Step 1 – Intralayer vdW Scan  (runs in the GUI)
# ══════════════════════════════════════════════════════════
with tab_step1_vdw:
    st.subheader("Run in this GUI")
    st.caption(
        "The rigid-sphere vdW scan completes on a laptop, so it runs right "
        "here (no HPC needed)."
    )
    c1, c2, c3 = st.columns(3)
    alpha_min = c1.number_input("alpha min (deg)", value=5.0, key="s1vdw_amin")
    alpha_max = c2.number_input("alpha max (deg)", value=85.0, key="s1vdw_amax")
    alpha_step = c3.number_input("alpha step (deg)", value=5.0, min_value=0.1, key="s1vdw_astep")
    vdw_tol = st.number_input("VdW tolerance (Å, overlap threshold)", value=0.0, key="s1vdw_tol")

    if st.button("Run vdW Scan", key="s1vdw_run"):
        st.info(
            "The scan engine is pending integration of Ono's vdW code "
            "(legacy/ono_scripts/). Once integrated, the scan will run here "
            "and its results will appear below."
        )

    st.divider()
    st.subheader("Results — alpha vs S = a×b")
    df, src = load_results_csv("s1vdw_results", "step1_vdw.csv")
    if df is not None:
        source_badge(src)
        line_plot_section(
            df, "s1vdw_plot",
            x_candidates=["alpha", "α"],
            y_candidates=["s", "a*b", "ab", "area"],
        )
        st.caption(
            "Clicking a point to preview the corresponding layer structure "
            "will be added together with the scan engine."
        )


# ══════════════════════════════════════════════════════════
#  Tab 3: Step 1 – DFT-D Optimization  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_step1_dft:
    st.subheader("How to run")
    cli_intro("Step 1b (DFT-D refinement of a, b at each alpha)")
    st.code(GAUSSIAN_KEYWORDS, language="text")

    st.divider()
    st.subheader("Results — alpha vs E_intra(8)")
    df, src = load_results_csv("s1dft_results", "step1_dft.csv")
    if df is not None:
        source_badge(src)
        line_plot_section(
            df, "s1dft_plot",
            x_candidates=["alpha", "α"],
            y_candidates=["e_intra8", "e_intra(8)", "e_intra", "e"],
        )
        st.caption("The energy minimum identifies the R-form (alpha ≈ 25°).")


# ══════════════════════════════════════════════════════════
#  Tab 4: Step 2 – Inclination Map  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_step2:
    st.subheader("How to run")
    cli_intro(
        "Step 2 (long-axis inclination map from the R-form)",
        note="Scan ranges per spec.md: theta_incl 0–40° (1° step), phi_incl 0–360° (5° step).",
    )
    st.code(GAUSSIAN_KEYWORDS, language="text")

    st.divider()
    st.subheader("Results — E_intra(6) map")
    df, src = load_results_csv("s2_results", "step2_map.csv")
    if df is not None:
        source_badge(src)
        heatmap_section(
            df, "s2_map",
            x_candidates=["theta_incl", "theta"],
            y_candidates=["phi_incl", "phi"],
            val_candidates=["e_intra6", "e_intra(6)", "e_intra", "e"],
            val_label="E_intra(6) (kcal/mol)",
        )
        st.caption(
            "G-form: phi_incl = 0°/180° (glide symmetry kept). "
            "N-form: phi_incl ≈ ±48°/±132° (glide symmetry broken). "
            "The polar-style axes from the paper (x = theta·cos(phi), "
            "y = theta·sin(phi)) will be added once the result format is fixed."
        )


# ══════════════════════════════════════════════════════════
#  Tab 5: Step 3 – Interlayer Stacking  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_step3:
    st.subheader("How to run")
    cli_intro("Step 3 (interlayer V(x,y) and E_inter(7) maps)")
    st.code(GAUSSIAN_KEYWORDS, language="text")

    st.divider()
    form = st.radio("Stacking dataset", ["N1 (→ Type II)", "N2 (→ Type IV)"],
                    horizontal=True, key="s3_form")
    sample_file = "step3_N1.csv" if form.startswith("N1") else "step3_N2.csv"

    st.subheader(f"Results — {form}")
    df, src = load_results_csv(f"s3_results_{sample_file}", sample_file)
    if df is not None:
        source_badge(src)
        heatmap_section(
            df, f"s3_map_{sample_file}",
            x_candidates=["x"],
            y_candidates=["y"],
            val_candidates=["e_inter7", "e_inter(7)", "e_inter", "e", "v"],
            val_label="E_inter(7) or V(x,y)",
        )
        st.caption(
            "Switch the Value column between V (vdW volume) and E_inter(7) "
            "if your CSV contains both."
        )


# ══════════════════════════════════════════════════════════
#  Tab 6: Step 4 – Refinement  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_step4:
    sub_twist, sub_incl = st.tabs([
        "6a. Twist (G-form → Type III)",
        "6b. Non-uniform Inclination (N-form → Type IV)",
    ])

    with sub_twist:
        st.subheader("How to run")
        cli_intro(
            "Step 4a (twist refinement of the G-form)",
            note="Expected: naphthalene ≈ 13°, anthracene ≈ 9°; no gain for tetracene and longer.",
        )
        st.divider()
        st.subheader("Results — theta_twist vs E_int(near)")
        df, src = load_results_csv("s4a_results", "step4a_twist.csv")
        if df is not None:
            source_badge(src)
            line_plot_section(
                df, "s4a_plot",
                x_candidates=["theta_twist", "twist"],
                y_candidates=["e_int_near", "e_int(near)", "e_int", "e"],
            )

    with sub_incl:
        st.subheader("How to run")
        cli_intro(
            "Step 4b (non-uniform inclination of the N-form)",
            note="Compares Type II vs Type IV; expected crossover at tetracene (theta'_incl ≈ 2.5°).",
        )
        st.divider()
        st.subheader("Results — theta'_incl vs E_int(near)")
        df, src = load_results_csv("s4b_results", "step4b_incl.csv")
        if df is not None:
            source_badge(src)
            line_plot_section(
                df, "s4b_plot",
                x_candidates=["theta_incl2", "theta_prime", "theta'_incl", "theta"],
                y_candidates=["e_int_near", "e_int(near)", "e_int", "e"],
            )


# ══════════════════════════════════════════════════════════
#  Tab 7: Transfer Integrals  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_step5:
    st.subheader("How to run")
    st.info(
        "Transfer integrals are computed from the command line with the "
        "CSV-batch wrapper around Prof. Matsui's HOMO-HOMO overlap code "
        "(Gaussian16, B3LYP/6-31g*), prepared by Ono — pending integration "
        "via legacy/ono_scripts/. The command will be documented here."
    )
    st.code("B3LYP/6-31g*  (HOMO → transfer integral J)", language="text")

    st.divider()
    st.subheader("Results — J per polymorph / contact type")
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
                x_candidates=["alpha", "theta_incl"],
                y_candidates=["j", "j_mev"],
            )
