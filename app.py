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
from plotly.subplots import make_subplots
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from csp.structure.molecule import load_molecule, _VDW  # noqa: E402
from csp.structure.intralayer import cluster9, cluster6_inclined, dimer, monomer_csv  # noqa: E402
from csp.vdw.contact import step1a_scan  # noqa: E402
from csp.vdw.interlayer import interlayer_vdw_scan, bilayer_preview  # noqa: E402
from csp.plot.viewer3d import to_xyz_string, render_3d_html  # noqa: E402
from csp.plot.map2d import build_heatmap_figure  # noqa: E402
from csp.plot.step1_results import (  # noqa: E402
    classify_and_fold_step1_results, KIND_LABEL as _S1R_KIND_LABEL, KIND_COLOR as _S1R_KIND_COLOR,
)
from csp.plot.step2_results import build_theta_phi_map  # noqa: E402

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
    "Polymorphism in Polyacenes\" (JACS, submitted). See the project "
    "README for a full tab-by-tab guide."
)


# ══════════════════════════════════════════════════════════
#  Shared helpers
# ══════════════════════════════════════════════════════════
def load_results_csv(
    key: str, sample_file: str, label: Optional[str] = None,
) -> tuple[Optional[pd.DataFrame], str]:
    """Results data source: user-dropped CSV first, bundled sample otherwise.

    Returns (df, source) where source is "uploaded", "sample", or "" when
    neither is available. Implements the sample-until-dropped display policy
    from spec.md "実行方式".
    """
    uploaded = st.file_uploader(
        label or "Drop your results CSV here (a precomputed sample is shown until then)",
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
        "Style", ["Space fill", "Capped sticks"], horizontal=True, key=style_key
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


_TREE_STEPWISE = """your-working-dir/
├── stepX_init_params.csv   # starting points -- from a GUI download, or hand-written
├── stepX.csv               # results, created/updated as jobs finish
├── gaussian/               # created automatically
│   ├── *.inp               # Gaussian16 input files
│   ├── *.log               # Gaussian16 output files
│   └── *.r1                # batch job scripts
└── gaussview/              # created automatically
    └── *.xyz               # structures, for visualization"""

_TREE_TCAL = """your-working-dir/
├── init_params.csv                # arrangements to compute (a, b, alpha, A2, z)
├── result.txt                     # summarized transfer integrals, after --result
└── a=.._b=.._alpha=.._A2=.._z=../ # one folder per arrangement
    ├── job.sh, tcal_1.py          # copied in automatically
    ├── test_t/test_p*.gjf         # Gaussian16 inputs (T-shaped / slipped-parallel)
    ├── test_t/test_p*.log         # Gaussian16 outputs
    └── test_t/test_p.txt          # per-arrangement transfer integral"""

# Not in Ono's own readmes (they only mention --monomer-name as a flag) --
# added separately since CSP_MONOMER_DIR/what's actually inside it was a
# real point of confusion (someone saw the 4-column lattice CSV and asked
# where the atom coordinates were).
_TREE_MONOMER = """/path/to/monomer/         <- CSP_MONOMER_DIR points here
├── pentacene.csv         # X, Y, Z, R columns -- download from the picker above
├── naphthalene.csv       # one CSV per molecule name you pass to --monomer-name
└── ..."""


def cli_howto(
    *,
    what: str,
    prepare: str,
    setup: str,
    command: str,
    output: str,
    show_route: bool = True,
    tree: Optional[str] = _TREE_STEPWISE,
    scripts_dir: str = LEGACY_DIR,
    scripts_files: str = "the `stepX.py`/`make_stepX.py` pair for this step, plus `utils.py`",
) -> None:
    """Standard 5-part "how to run this step with Ono's package" block.

    Every CLI section (Step 1 DFT, Step 2/3 para/twist, Transfer Integrals)
    uses the same shape so a user who has read one knows where to look in
    the others: where the scripts actually are inside the downloaded
    folder, what input file(s) to prepare by hand (and where their values
    come from), what one-time environment setup is needed, then the
    command and the resulting output file.

    `tree`: an optional `--auto-dir` working-directory layout to show in a
    collapsed reference section (pass `_TREE_TCAL` for the transfer-integral
    workflow's per-arrangement layout, or `None` to omit).
    """
    st.info(what)
    st.markdown("**1. Locate the scripts**")
    st.markdown(
        f"They live at `{scripts_dir}/` inside the folder you downloaded "
        f"or cloned (e.g. `ono-csp-main/{scripts_dir}/` from a zip "
        f"download, or `ono-csp/{scripts_dir}/` from `git clone`). Either "
        f"`cd` there directly, or copy {scripts_files} to wherever you'll "
        "actually run this (e.g. your HPC login node)."
    )
    st.markdown("**2. Prepare inputs**")
    st.markdown(prepare)
    st.markdown("**3. One-time setup**")
    st.markdown(setup)
    st.markdown("**4. Run**")
    st.code(command, language="bash")
    if show_route:
        st.code(GAUSSIAN_ROUTE, language="text")
    st.markdown("**5. Output**")
    st.markdown(output)
    if tree:
        with st.expander("What the working directory (`--auto-dir`) looks like"):
            st.code(tree, language="text")
            st.caption("...and the separate monomer directory (`CSP_MONOMER_DIR`) referenced above:")
            st.code(_TREE_MONOMER, language="text")


_SETUP_MONOMER_ENV = (
    "- Set `CSP_MONOMER_DIR` to the folder holding the monomer CSV "
    "(X, Y, Z, R columns — download it from **Tab 1**), e.g. "
    "`export CSP_MONOMER_DIR=/path/to/data/monomer` "
    "(defaults to `~/path/to/monomer/` if unset). This matters because "
    "the init_params CSV below only has the *lattice* numbers (a, b, "
    "alpha, ...) — no atom coordinates. The script reads the monomer's "
    "actual X, Y, Z, R from this directory and combines the two to build "
    "the real 3D dimer/cluster geometry that Gaussian16 runs on."
)
_SETUP_SCHEDULER = (
    "- The legacy scripts as uploaded submit jobs with `pjsub` "
    "(Fujitsu/PJM, the scheduler on the cluster they were originally "
    "written for). They're kept as-is on purpose — adapt the "
    "batch-script generation and submission call to whatever scheduler "
    "*your* cluster actually uses (`qsub` for SGE/PBS, `sbatch` for "
    "Slurm, etc.), not necessarily `qsub`. `src/csp/dft/job_cluster.py` "
    "has a worked example for one SGE cluster if useful as a reference."
)


# ══════════════════════════════════════════════════════════
#  Persistent molecule picker — visible on every tab, so switching the
#  active molecule doesn't require flipping back to Tab 1.
# ══════════════════════════════════════════════════════════
available_presets = [
    m for m in PRESET_MOLECULES if (MOLECULE_DIR / f"{m}.xyz").exists()
]
pick_col1, pick_col2, pick_col3 = st.columns([1, 2, 3])
molecule = None
with pick_col1:
    mol_source = st.radio(
        "Source", ["Preset", "Upload custom XYZ"], horizontal=True, key="mol_source_top",
    )
with pick_col2:
    if mol_source == "Preset":
        if not available_presets:
            st.warning(f"No preset XYZ files found in {MOLECULE_DIR}")
        else:
            preset_name = st.selectbox("Molecule", available_presets, key="mol_preset_top")
            molecule = load_molecule(preset_name, molecule_dir=MOLECULE_DIR)
    else:
        uploaded = st.file_uploader(
            "Upload a monomer XYZ file", type=["xyz"], key="mol_upload_top",
        )
        if uploaded is not None:
            with tempfile.NamedTemporaryFile(
                suffix=".xyz", mode="w", delete=False, encoding="utf-8"
            ) as tf:
                tf.write(uploaded.getvalue().decode("utf-8"))
                tmp_path = tf.name
            molecule = load_molecule(tmp_path)
            molecule.name = Path(uploaded.name).stem
with pick_col3:
    if molecule is not None:
        st.session_state["molecule"] = molecule
        st.caption(f"Active molecule: **{molecule.name}** — {molecule.n_atoms} atoms")
    else:
        st.caption("No molecule selected yet.")
st.divider()


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
    st.caption(
        "The molecule picker above (visible on every tab) is the shared "
        "source for Tabs 2-5 below -- no need to come back here just to "
        "switch molecules. This tab holds the less-frequently-changed "
        "settings: exporting the monomer CSV, the vdW radius table, and a "
        "3D sanity-check view."
    )
    col_pick, col_view = st.columns([1, 1])
    molecule = st.session_state.get("molecule")

    with col_pick:
        if molecule is not None:
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
            st.caption("Select or upload a molecule above to preview it here.")
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
        "Rigid-sphere vdW contact model — a fast, vectorized "
        "reimplementation (not a call into legacy step1.py itself; "
        "verified to reproduce its `--init` output numerically) that "
        "runs directly in this GUI: for each herringbone half-angle "
        "alpha, sweep the T-contact direction and find the smallest "
        "cell S = a×b compatible with T and slipped-parallel contacts. "
        "Initial DFT candidates are the local minima of S plus both "
        "endpoints."
    )
    c1, c2, c3 = st.columns(3)
    alpha_min = c1.number_input("alpha min (deg)", value=5.0, key="s1vdw_amin")
    alpha_max = c2.number_input("alpha max (deg)", value=45.0, key="s1vdw_amax")
    alpha_step = c3.number_input("alpha step (deg)", value=5.0, min_value=0.5, key="s1vdw_astep")

    if mol2 is None:
        st.info("Select a molecule using the picker above first.")
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
            # Accept the current "alpha" column name or the legacy "theta"
            # one (Ono's scripts originally called this column theta; both
            # bundled samples and freshly-generated files now say alpha).
            if "alpha" not in df_up.columns and "theta" in df_up.columns:
                df_up = df_up.rename(columns={"theta": "alpha"})
            if {"a", "b", "alpha"} <= set(df_up.columns):
                df_init = df_up.rename(columns={"s": "S"})
                if "S" not in df_init.columns:
                    df_init["S"] = df_init["a"] * df_init["b"]
                if src2 == "sample":
                    # step1_init_params.csv alone can't feed the Fig. S1(c)
                    # sub-tab (needs the full per-alpha beta sweep, not just
                    # the extracted candidates) -- load the matching bundled
                    # curves file directly if one's there (no separate
                    # uploader for it; it's a fixed companion to the sample).
                    curves_path = EXAMPLE_DIR / "step1_vdw_curves.csv"
                    if curves_path.exists():
                        try:
                            df_curves = pd.read_csv(curves_path)
                        except Exception:
                            df_curves = None
            else:
                st.warning("Expected columns a, b, alpha (legacy step1_init_params.csv format).")

    if df_init is not None:
        source_badge(src2)

        _KIND_LABEL = {
            "b_contact": "HB",
            "a_contact": "PS",
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
            row = df.loc[df["S"].idxmin()] if "S" in df.columns else df.iloc[0]
            return {
                "alpha": float(row["alpha"]), "a": float(row["a"]), "b": float(row["b"]),
                "label": f"default (min S): alpha={row['alpha']} a={row['a']} b={row['b']}",
            }

        def _render_preview(df: pd.DataFrame, key_suffix: str) -> None:
            st.markdown("**Structure preview**")
            if mol2 is None:
                st.caption("Select a molecule using the picker above.")
                return
            current = st.session_state.get("s1vdw_current")
            if current is None:
                current = _default_current(df)
                if current is not None:
                    # Same reasoning as the DFT Fig. 2(b) default: persist it
                    # so Tab 3/Tab 4 inherit "the min-S candidate" instead of
                    # an arbitrary placeholder before any click happens here.
                    st.session_state["s1vdw_current"] = current
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
                fold_vdw = st.checkbox(
                    "Reflect alpha -> 90-alpha (a/b swapped) to show the full "
                    "0-90° range", value=True, key="s1vdw_fold",
                )
                df_init_plot = df_init
                env = None
                if df_curves is not None:
                    env = df_curves[df_curves["valid"]].groupby("alpha", as_index=False)["S"].min()
                if fold_vdw and "kind" in df_init.columns:
                    folded = df_init.copy()
                    folded["alpha"] = 90.0 - folded["alpha"]
                    folded["a"], folded["b"] = df_init["b"].values, df_init["a"].values
                    # `kind` stays fixed under the fold -- same physical
                    # branch (HB/PS), just axes relabeled; see
                    # step1_results.classify_and_fold_step1_results's
                    # docstring for why re-deriving it from the swapped a/b
                    # would flip HB/PS identity for half the displayed range.
                    df_init_plot = pd.concat([df_init, folded], ignore_index=True).reset_index(drop=True)
                    if env is not None:
                        env_folded = env.copy()
                        env_folded["alpha"] = 90.0 - env_folded["alpha"]
                        env = pd.concat([env, env_folded], ignore_index=True).sort_values("alpha")
                figA = go.Figure()
                if env is not None:
                    figA.add_trace(go.Scatter(
                        x=env["alpha"], y=env["S"], mode="lines",
                        name="min S (feasible)", line=dict(color="lightgray"),
                        hoverinfo="skip",
                    ))
                for kind, grp in (df_init_plot.groupby("kind") if "kind" in df_init_plot.columns else []):
                    grp = grp.sort_values("alpha")
                    color = _KIND_COLOR.get(kind, "gray")
                    figA.add_trace(go.Scatter(
                        x=grp["alpha"], y=grp["S"], mode="lines+markers",
                        name=_KIND_LABEL.get(kind, kind),
                        line=dict(color=color), marker=dict(size=9, color=color),
                        customdata=grp.index.to_numpy(),
                        hovertemplate="alpha=%{x}<br>S=%{y:.1f} Å²<extra>" + kind + "</extra>",
                    ))
                if "kind" not in df_init.columns:
                    grp = df_init.sort_values("alpha")
                    figA.add_trace(go.Scatter(
                        x=grp["alpha"], y=grp["S"], mode="lines+markers",
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
                        row = df_init_plot.loc[idx]
                        st.session_state["s1vdw_current"] = {
                            "alpha": float(row["alpha"]), "a": float(row["a"]), "b": float(row["b"]),
                            "label": f"Fig.2b click: alpha={row['alpha']} a={row['a']} b={row['b']}",
                        }
                        st.session_state["s1vdw_figA_prev"] = idx

                st.download_button(
                    "Download step1_init_params.csv (input for the DFT step below)",
                    data=df_init[["a", "b", "alpha"]].assign(status="NotYet").to_csv(index=False),
                    file_name="step1_init_params.csv", mime="text/csv",
                    key="dl_s1vdw_init",
                )
                with st.expander("Candidates table"):
                    st.dataframe(df_init, width="stretch")

                labels = [
                    f"alpha={r['alpha']}  a={r['a']}  b={r['b']}  S={r['S']:.1f}"
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
                            "alpha": float(row["alpha"]), "a": float(row["a"]), "b": float(row["b"]),
                            "label": f"list pick: alpha={row['alpha']} a={row['a']} b={row['b']}",
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
                        "slipped-parallel-neighbor vdW spheres clear (valid); **long-dashed, "
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
    st.subheader("How to run (CLI, i.e. commands typed in a terminal)")
    cli_howto(
        what=(
            "Step 1 DFT-D refinement runs on an HPC with Gaussian16 using "
            f"Ono's scripts in `{LEGACY_DIR}/` (`step1.py` + `make_step1.py`). "
            "It hill-climbs (a, b) in 0.1 Å steps at each fixed alpha, "
            "minimizing E_intra(8) = 4·E_t + 2·E_p1 + 2·E_p2 (BSSE-corrected "
            "dimer interaction energies)."
        ),
        prepare=(
            "- `step1_init_params.csv` (columns: `a, b, alpha, status`) — "
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
        output="`<auto-dir>/step1.csv` — columns `a, b, alpha, E, E_p1, E_p2, E_t, status, file_name`.",
        scripts_files="`step1.py`, `make_step1.py`, and `utils.py`",
    )

    st.divider()
    st.subheader("DFT results — alpha vs E_intra(8)  (step1.csv)")
    df, src = load_results_csv("s1dft_results", "step1.csv")
    if df is not None:
        source_badge(src)
        st.caption(
            "E = 4·E_t + 2·E_p1 + 2·E_p2 (kcal/mol). Its global minimum "
            "identifies the R-form (alpha ≈ 25°)."
        )

        st.caption(
            "step1.csv alone doesn't record which structural branch "
            "(HB / PS / local min) each hill-climb point belongs to. This "
            "is recovered directly from the data: at "
            "each alpha, the visited (a, b) points are treated as a 2D "
            "energy grid, and its local minima *are* the converged "
            "branches (a DFT hill-climb only ever stops at one). "
            "Polyacenes are also exactly symmetric under alpha -> "
            "90-alpha with a/b swapped (same physical structure, axes "
            "relabeled), so a scan run only up to alpha=45 is reflected "
            "to cover the full range, matching the paper's Fig. 2(b). "
            "Click a point for its 9-molecule structure."
        )
        if "alpha" not in df.columns and "theta" not in df.columns:
            st.warning("Expected an 'alpha' column in this CSV to classify branches.")
        else:
            df_branches = classify_and_fold_step1_results(df)
            if len(df_branches) == 0:
                st.warning("No branches found to classify in this CSV.")
            else:
                col_fig2b, col_3d_fig2b = st.columns([2, 1])
                with col_fig2b:
                    fig2b = go.Figure()
                    for kind, grp in df_branches.groupby("kind"):
                        grp = grp.sort_values("alpha")
                        color = _S1R_KIND_COLOR.get(kind, "gray")
                        fig2b.add_trace(go.Scatter(
                            x=grp["alpha"], y=grp["E"], mode="lines+markers",
                            name=_S1R_KIND_LABEL.get(kind, kind),
                            line=dict(color=color), marker=dict(size=8, color=color),
                            hovertemplate="alpha=%{x}<br>E=%{y:.2f}<extra>"
                                          + _S1R_KIND_LABEL.get(kind, kind) + "</extra>",
                        ))
                    fig2b.update_layout(
                        xaxis_title="alpha (deg)", yaxis_title="E_intra(8) (kcal/mol)",
                        margin=dict(l=20, r=20, t=30, b=20),
                    )
                    event_2b = st.plotly_chart(
                        fig2b, width="stretch", on_select="rerun", key="s1fig2b_chart"
                    )
                    pts_2b = event_2b.selection.points if (event_2b and event_2b.selection) else []
                    if pts_2b:
                        p0 = pts_2b[0]
                        alpha_sel, E_sel = p0.get("x"), p0.get("y")
                        if alpha_sel is not None and E_sel is not None:
                            match = df_branches[
                                np.isclose(df_branches["alpha"], alpha_sel)
                                & np.isclose(df_branches["E"], E_sel)
                            ]
                            if len(match):
                                r = match.iloc[0]
                                st.session_state["s1fig2b_current"] = {
                                    "alpha": float(r["alpha"]), "a": float(r["a"]), "b": float(r["b"]),
                                    "label": f"clicked [{_S1R_KIND_LABEL.get(r['kind'], r['kind'])}]: "
                                             f"alpha={r['alpha']} a={r['a']} b={r['b']} E={r['E']:.2f}",
                                }
                    with st.expander("Branches table"):
                        st.dataframe(df_branches, width="stretch")
                with col_3d_fig2b:
                    st.markdown("**Structure preview**")
                    if mol2 is None:
                        st.caption("Select a molecule using the picker above for the 3D preview.")
                    else:
                        current_2b = st.session_state.get("s1fig2b_current")
                        if current_2b is None and len(df_branches):
                            # The global minimum is physically identical whether
                            # shown via its original alpha (<=45) or its folded
                            # mirror (alpha -> 90-alpha, a/b swapped) -- prefer
                            # the unfolded one so the displayed alpha matches the
                            # commonly-cited value (e.g. ~25 deg) instead of its
                            # (equally correct but less recognizable) >45 mirror.
                            e_min = df_branches["E"].min()
                            tied = df_branches[np.isclose(df_branches["E"], e_min)]
                            unfolded_tied = tied[~tied["folded"]]
                            best = (unfolded_tied if len(unfolded_tied) else tied).iloc[0]
                            current_2b = {
                                "alpha": float(best["alpha"]), "a": float(best["a"]), "b": float(best["b"]),
                                "label": f"default (min E) [{_S1R_KIND_LABEL.get(best['kind'], best['kind'])}]: "
                                         f"alpha={best['alpha']} a={best['a']} b={best['b']}",
                            }
                            # Persist this so Tab 3/Tab 4's own a/b/alpha defaults
                            # pick up the actual R-form optimum immediately, even
                            # before the user clicks a point here -- otherwise
                            # they'd fall back to an arbitrary placeholder instead
                            # of "whatever this tab is already showing".
                            st.session_state["s1fig2b_current"] = current_2b
                        if current_2b is not None:
                            st.caption(current_2b["label"])
                            syms_2b, coords_2b = cluster9(mol2, current_2b["a"], current_2b["b"], current_2b["alpha"])
                            render_molecule_3d(
                                syms_2b, coords_2b,
                                f"{mol2.name} cluster9 alpha={current_2b['alpha']} "
                                f"a={current_2b['a']} b={current_2b['b']}",
                                key_suffix="s1fig2b", style_key="s1fig2b_style",
                            )


# ══════════════════════════════════════════════════════════
#  Tab 3: Step 2 – Long-Axis Shift  (para / twist sub-tabs)
# ══════════════════════════════════════════════════════════
with tab_step2:
    sub_para, sub_twist = st.tabs(["para", "twist (→ Type III)"])

    with sub_para:
        st.subheader("How to run (CLI, i.e. commands typed in a terminal)")
        cli_howto(
            what=(
                "Step 2 (para variant) scans the molecular shift z along the "
                "long axis at the T-shaped and slipped-parallel contacts, "
                "with the lattice fixed to the Step 1 optimum "
                f"(`{LEGACY_DIR}/step2_para.py` + `make_step2_para.py`; "
                "41 z-points from 0 to 4 Å in one Gaussian input via "
                "--Link1). This directly implements the paper's §2.2 "
                "\"Optimization by Uniform Long-Axis Inclination\": z here "
                "*is* the physical long-axis slide (ΔzT) the paper describes "
                "(Fig. 5a) — the paper's inclination angle theta_incl is "
                "just this same z reparametrized as an angle via "
                "ΔzT = (b/2)·tan(theta_incl) (SI Fig. S6a), with the "
                "G-form minimum expected around ΔzT ≈ 2.4 Å (well within "
                "this scan's 0–4 Å range). Applying the *same* z to both "
                "T-shaped and slipped-parallel contacts keeps glide symmetry "
                "(phi_incl = 0°/90°, the G-form direction); finding the "
                "N-form (broken glide symmetry) needs *different* shifts "
                "per T-contact — see Step 3's Rt/Rp in **Tab 4**."
            ),
            prepare=(
                "- `step2_init_params.csv`: **one row**, columns `a, b, "
                "alpha` — take the best (a, b) at your chosen alpha from "
                "**Tab 2**'s DFT results (`step1.csv`) and write it "
                "yourself (no generator script for this one)."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command="python step2_para.py --auto-dir /path/to/workdir --monomer-name pentacene",
            output="`<auto-dir>/step2_para.csv` — columns `z, Et, Ep` (mirrored to -4..4 Å).",
            scripts_files="`step2_para.py`, `make_step2_para.py`, and `utils.py`",
        )

        st.divider()
        st.subheader("Results — Et / Ep vs z  (step2_para.csv)")
        df, src = load_results_csv("s2_results", "step2_para.csv")
        if df is not None:
            source_badge(src)
            cols_lower = {c.lower(): c for c in df.columns}
            if {"z", "et", "ep"} <= set(cols_lower):
                zc, etc, epc = cols_lower["z"], cols_lower["et"], cols_lower["ep"]
                st.caption(
                    "Et(z) and Ep(z) are independent dimer scans (T-shaped "
                    "neighbor at (a/2, b/2, z); slipped-parallel neighbor at "
                    "(0, b, z)). Note there's no single '4·Et+2·Ep' curve to "
                    "show here in general -- that combination only holds for "
                    "the glide-symmetric G-form (same z on every contact); "
                    "the Fig. 5(b)-style map below combines Et/Ep at "
                    "*independent* z per contact to cover the N-form too."
                )
                _s2_default = (
                    st.session_state.get("s1fig2b_current")
                    or st.session_state.get("s1vdw_current") or {}
                )
                # Fallback defaults (when Tab 2 hasn't been used this session):
                # pentacene's Type II parameters from the paper's SI Table S2
                # (a=7.2, b=5.9, alpha=25 deg, theta_incl=27, phi_incl=48) --
                # a physically real point instead of an arbitrary placeholder,
                # so the sample 3D previews/maps show something meaningful.
                c1, c2, c3 = st.columns(3)
                s2_a = c1.number_input("a (Å)", value=float(_s2_default.get("a", 7.2)), key="s2_a")
                s2_b = c2.number_input("b (Å)", value=float(_s2_default.get("b", 5.9)), key="s2_b")
                s2_theta = c3.number_input(
                    "alpha (deg, from Step 1)",
                    value=float(_s2_default.get("alpha", _s2_default.get("theta", 25.0))),
                    key="s2_theta",
                )

                col_s2plot, col_s2_3d = st.columns([2, 1])
                with col_s2plot:
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(
                        x=df[zc], y=df[etc], mode="lines+markers", name="E_t (T-shaped)",
                    ))
                    fig.add_trace(go.Scatter(
                        x=df[zc], y=df[epc], mode="lines+markers", name="E_p (slipped parallel)",
                    ))
                    fig.update_layout(
                        xaxis_title="z shift along long axis (Å)",
                        yaxis_title="E (kcal/mol)",
                        margin=dict(l=20, r=20, t=30, b=20),
                    )
                    event_s2 = st.plotly_chart(
                        fig, width="stretch", on_select="rerun", key="s2_chart"
                    )
                    pts_s2 = event_s2.selection.points if (event_s2 and event_s2.selection) else []
                    if pts_s2:
                        p0 = pts_s2[0]
                        z_sel = p0.get("x")
                        curve_num = p0.get("curve_number", 0)
                        if z_sel is not None:
                            # step2_para.py's own rule: the slipped-parallel
                            # contact runs along whichever of a/b is the
                            # *shorter* axis (along b when a>b, else along a).
                            kind = "t" if curve_num == 0 else ("p1" if s2_a > s2_b else "p2")
                            st.session_state["s2_current"] = {
                                "kind": kind, "z": float(z_sel),
                                "label": f"clicked {'E_t' if kind == 't' else 'E_p'}: z={z_sel}",
                            }
                    with st.expander("Data table"):
                        st.dataframe(df, width="stretch")
                with col_s2_3d:
                    st.markdown("**Structure preview**")
                    if mol2 is None:
                        st.caption("Select a molecule using the picker above.")
                    else:
                        current_s2 = st.session_state.get("s2_current") or {
                            "kind": "t", "z": 0.0, "label": "default: E_t at z=0",
                        }
                        st.caption(current_s2["label"])
                        c_i, c_j = dimer(
                            mol2, current_s2["kind"], s2_a, s2_b, s2_theta, z=current_s2["z"],
                        )
                        syms_s2 = list(mol2.symbols) * 2
                        render_molecule_3d(
                            syms_s2,
                            np.vstack([c_i, c_j]),
                            f"{mol2.name} dimer {current_s2['kind']} z={current_s2['z']}",
                            key_suffix="s2dimer", style_key="s2dimer_style",
                        )

                st.divider()
                st.subheader("Fig. 5(b)-style: theta_incl / phi_incl map (Eintra(6))")
                st.caption(
                    "Reconstructed by combining Et(z)/Ep(z) at independently "
                    "chosen z per contact (zt for the T-neighbor at (a/2, b/2), "
                    "zp for the slipped-parallel neighbor at (0, b)): "
                    "Eintra(6) = 2·(Et(zt) + Et(zt−zp) + Ep(zp)). No new DFT "
                    "needed. zp=0 is the glide-symmetric G-form direction "
                    "(matches the plain Et/Ep scan above); zt != zp/2 breaks "
                    "glide symmetry (N-form). Click a point for its 6-molecule "
                    "structure."
                )
                df_map = build_theta_phi_map(df.rename(columns={zc: "z", etc: "Et", epc: "Ep"}), s2_a, s2_b)
                if len(df_map) == 0:
                    st.warning("Not enough z-range overlap to build this map.")
                else:
                    col_map, col_map3d = st.columns([2, 1])
                    with col_map:
                        n_side = max(1, int(round(len(df_map) ** 0.5)))
                        marker_px = max(3, 650 / n_side)
                        # R/G/N-form local minima aren't necessarily the
                        # global minimum (the paper shows them as separate
                        # local minima on this same landscape) -- a plain
                        # color scale washes them out against the flat
                        # R-form's much deeper global minimum, so find and
                        # mark them explicitly (same approach as Tab 4's vdW
                        # map: 2D minimum_filter over the (zt, zp) grid).
                        from scipy.ndimage import minimum_filter
                        pivot_map = df_map.pivot(index="zp", columns="zt", values="E")
                        grid = pivot_map.values
                        is_min = (grid == minimum_filter(grid, size=5)) & np.isfinite(grid)
                        min_i, min_j = np.where(is_min)
                        min_zp = pivot_map.index.to_numpy()[min_i]
                        min_zt = pivot_map.columns.to_numpy()[min_j]
                        min_pairs = set(zip(min_zt, min_zp))
                        min_mask = [
                            (zt, zp) in min_pairs for zt, zp in zip(df_map["zt"], df_map["zp"])
                        ]
                        min_rows = df_map[min_mask]
                        fig5b = go.Figure()
                        fig5b.add_trace(go.Scatter(
                            x=df_map["x"], y=df_map["y"], mode="markers",
                            marker=dict(
                                symbol="square", size=marker_px,
                                color=df_map["E"], colorscale="RdBu_r",
                                colorbar=dict(title="E_intra(6)"),
                            ),
                            showlegend=False,
                            hovertemplate="x=%{x:.1f}<br>y=%{y:.1f}<br>E=%{marker.color:.2f}<extra></extra>",
                        ))
                        fig5b.add_trace(go.Scatter(
                            x=min_rows["x"], y=min_rows["y"], mode="markers",
                            name="local min", showlegend=False,
                            marker=dict(symbol="square-open", size=marker_px + 6, color="gold", line=dict(width=2)),
                            hovertemplate="x=%{x:.1f}<br>y=%{y:.1f}<extra>local min</extra>",
                        ))
                        fig5b.update_layout(
                            # Ono's plot2d() axes: x = theta_incl*cos(phi_incl),
                            # y = theta_incl*sin(phi_incl) (a polar->Cartesian
                            # transform, radius=theta_incl, angle=phi_incl) --
                            # matches the paper's Fig. 5(b) layout (x in
                            # [-45,45], y in [-30,30]), unlike plain
                            # (phi_incl, theta_incl) axes.
                            xaxis_title="theta_incl · cos(phi_incl) (deg)",
                            yaxis_title="theta_incl · sin(phi_incl) (deg)",
                            margin=dict(l=20, r=20, t=30, b=20),
                        )
                        fig5b.update_yaxes(scaleanchor="x", scaleratio=1)
                        event_5b = st.plotly_chart(
                            fig5b, width="stretch", on_select="rerun", key="s2fig5b_chart"
                        )
                        pts_5b = event_5b.selection.points if (event_5b and event_5b.selection) else []
                        if pts_5b:
                            # customdata proved unreliable for click
                            # identification elsewhere in this app (Tab 4's
                            # vdW map) -- match the point's plain x/y back to
                            # df_map's own x/y columns instead, which is the
                            # pattern that's actually held up everywhere else.
                            x_sel, y_sel = pts_5b[0].get("x"), pts_5b[0].get("y")
                            if x_sel is not None and y_sel is not None:
                                match5b = df_map[
                                    np.isclose(df_map["x"], x_sel) & np.isclose(df_map["y"], y_sel)
                                ]
                                if len(match5b):
                                    r5b = match5b.iloc[0]
                                    st.session_state["s2fig5b_current"] = {
                                        "zt": float(r5b["zt"]), "zp": float(r5b["zp"]),
                                        "label": f"clicked: zt={r5b['zt']} zp={r5b['zp']}",
                                    }
                        st.caption(f"{len(min_rows)} local minima marked (black squares) out of {len(df_map)} grid points.")
                        with st.expander("Map data table"):
                            st.dataframe(df_map, width="stretch")
                        with st.expander("Local minima table"):
                            st.dataframe(
                                min_rows[["zt", "zp", "theta_incl", "phi_incl", "E"]].sort_values("E"),
                                width="stretch",
                            )
                    with col_map3d:
                        st.markdown("**Structure preview**")
                        if mol2 is None:
                            st.caption("Select a molecule using the picker above.")
                        else:
                            current_5b = st.session_state.get("s2fig5b_current") or {
                                "zt": 0.0, "zp": 0.0, "label": "default: zt=0 zp=0 (flat R-form)",
                            }
                            st.caption(current_5b["label"])
                            syms6, coords6 = cluster6_inclined(
                                mol2, s2_a, s2_b, s2_theta, current_5b["zt"], current_5b["zp"],
                            )
                            render_molecule_3d(
                                syms6, coords6,
                                f"{mol2.name} cluster6 zt={current_5b['zt']} zp={current_5b['zp']}",
                                key_suffix="s2fig5b", style_key="s2fig5b_style",
                            )
            else:
                line_plot_section(
                    df, "s2_plot",
                    x_candidates=["z"], y_candidates=["et", "ep", "e"],
                )

    with sub_twist:
        # Independent fallback defaults (don't rely on sub_para's s2_a/s2_b/
        # s2_theta -- those are only defined when step2_para.csv is loaded
        # there, which isn't guaranteed if this sub-tab is opened first).
        _twist_default = (
            st.session_state.get("s1fig2b_current")
            or st.session_state.get("s1vdw_current") or {}
        )
        _twist_a = float(_twist_default.get("a", 6.6))
        _twist_b = float(_twist_default.get("b", 6.7))
        _twist_theta = float(_twist_default.get("alpha", _twist_default.get("theta", 21.0)))

        st.subheader("How to run (CLI, i.e. commands typed in a terminal)")
        cli_howto(
            what=(
                "Step 2 (twist variant, → Type III packing with glide "
                "symmetry): introduces the T-contact long-axis shift Rt and "
                "the torsion A2, then re-optimizes (a, b) by hill-climbing "
                f"at each fixed (alpha, Rt, A2) (`{LEGACY_DIR}/step2_twist.py`; "
                "E = 4·E_t + 2·E_p1)."
            ),
            prepare=(
                "- `step2_twist_init_params.csv`: rows of `alpha, Rt, A2, "
                "a, b, status`. `alpha` and starting `a, b` come from Step 1 "
                "(**Tab 2**'s DFT results); `Rt` (twist shift) and `A2` "
                "(torsion) are the new values you want to try (e.g. a grid "
                "of A2 = 0..20°, since the paper's §2.4 reports optima "
                "around 13° for naphthalene and 9° for anthracene). Set "
                "`status='NotYet'` for each row."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command=(
                "python step2_twist.py --auto-dir /path/to/workdir "
                "--monomer-name naphthalene \\\n    --num-nodes 4 --num-init 2"
            ),
            output=(
                "`<auto-dir>/step2_twist.csv` — columns `a, b, alpha, Rt, A2, "
                "E, E_p1, E_t, status, file_name`."
            ),
            scripts_files="`step2_twist.py`, `make_step2_twist.py`, and `utils.py`",
        )
        st.divider()
        st.subheader("Fig. 7(c)-style: Eintra(6) vs theta_twist / Rt (ΔzT)")
        df, src = load_results_csv("s4a_results", "step2_twist.csv")
        if df is not None:
            source_badge(src)
            st.caption(
                "Expected: E minimum at a twist of ≈ 13° for naphthalene, "
                "≈ 9° for anthracene; no gain for tetracene and longer "
                "(paper §2.4). At each (A2, Rt), (a, b) may have been "
                "hill-climbed to convergence, or evaluated once at a fixed "
                "starting point (step2_twist.py's get_opt_params_dict with "
                "its neighbor step set to 0, see README) -- either way, the "
                "map below takes the minimum E recorded at each grid point. "
                "Click a point for its T-shaped dimer at that twist/shift."
            )
            if src == "sample":
                st.info(
                    "This sample is Ono's real naphthalene data (a, b fixed "
                    "at the vdW-derived value, single-point energy at every "
                    "A2/Rt grid point — the Fig. 7(c) method above). "
                    "Reproducing the *rest* of the paper's Fig. 7 pipeline "
                    "is more involved and isn't fully wired up in this GUI "
                    "yet: from a map like this, Fig. 7(d) narrows to the "
                    "promising A2 range (paper: 0-12°) and, for each A2, "
                    "picks the Rt where E is lowest; at *those* (A2, Rt) "
                    "pairs only, (a, b) is then re-optimized by hill-climbing "
                    "(a separate, denser CSV — not this sample); finally the "
                    "interlayer energy Eint(near) is computed on top of that "
                    "re-optimized intralayer structure (Tab 4's twist "
                    "sub-tab) and added to Eintra(6) to get the Fig. 7(d) "
                    "total. Tab 4 twist's sample below is that last "
                    "interlayer piece — on its own it is not yet the full "
                    "Fig. 7(d) total energy."
                )
            cols_lower_t = {c.lower(): c for c in df.columns}
            need = {"a2", "rt", "e"}
            if not need <= set(cols_lower_t):
                st.warning(f"Expected columns A2, Rt, E (got {list(df.columns)}); "
                           "showing a generic plot instead.")
                line_plot_section(df, "s4a_plot", x_candidates=["a2", "rt"], y_candidates=["e"])
            else:
                a2c, rtc, ec = cols_lower_t["a2"], cols_lower_t["rt"], cols_lower_t["e"]
                ac = cols_lower_t.get("a")
                bc = cols_lower_t.get("b")
                thc = cols_lower_t.get("alpha", cols_lower_t.get("theta"))
                df_conv = df.loc[df.groupby([a2c, rtc])[ec].idxmin()].reset_index(drop=True)
                col_twist, col_twist3d = st.columns([2, 1])
                with col_twist:
                    pivot_t = df_conv.pivot(index=rtc, columns=a2c, values=ec)
                    grid_t = pivot_t.values
                    if grid_t.shape[0] >= 5 and grid_t.shape[1] >= 5:
                        from scipy.ndimage import minimum_filter
                        is_min_t = (grid_t == minimum_filter(grid_t, size=5)) & np.isfinite(grid_t)
                    else:
                        is_min_t = grid_t == np.nanmin(grid_t)
                    mi, mj = np.where(is_min_t)
                    min_rt = pivot_t.index.to_numpy()[mi]
                    min_a2 = pivot_t.columns.to_numpy()[mj]
                    # go.Heatmap (not scaled scatter markers) tiles the grid
                    # exactly regardless of how uneven the A2/Rt axis spacing
                    # or grid shape is -- fixed-pixel square markers left
                    # visible gaps whenever the grid wasn't square.
                    figT = go.Figure(go.Heatmap(
                        x=pivot_t.columns, y=pivot_t.index, z=grid_t,
                        colorscale="RdBu_r", colorbar=dict(title="E_intra(6)"),
                        hovertemplate="A2=%{x}<br>Rt=%{y}<br>E=%{z:.2f}<extra></extra>",
                    ))
                    figT.add_trace(go.Scatter(
                        x=min_a2, y=min_rt, mode="markers",
                        showlegend=False,
                        marker=dict(symbol="square-open", size=14, color="gold", line=dict(width=2)),
                        hovertemplate="A2=%{x}<br>Rt=%{y}<extra>local min</extra>",
                    ))
                    figT.update_layout(
                        xaxis_title="theta_twist / A2 (deg)", yaxis_title="Rt / ΔzT (Å)",
                        margin=dict(l=20, r=20, t=30, b=20),
                    )
                    event_t = st.plotly_chart(figT, width="stretch", on_select="rerun", key="s4a_chart")
                    pts_t = event_t.selection.points if (event_t and event_t.selection) else []
                    if pts_t:
                        # Match plain x/y back to df_conv's own columns
                        # instead of relying on customdata (unreliable for
                        # click identification elsewhere in this app).
                        x_sel_t, y_sel_t = pts_t[0].get("x"), pts_t[0].get("y")
                        if x_sel_t is not None and y_sel_t is not None:
                            match_t = df_conv[
                                np.isclose(df_conv[a2c], x_sel_t) & np.isclose(df_conv[rtc], y_sel_t)
                            ]
                        else:
                            match_t = df_conv.iloc[0:0]
                        if len(match_t):
                            r = match_t.iloc[0]
                            st.session_state["s4a_current"] = {
                                "a2": float(r[a2c]), "rt": float(r[rtc]),
                                "a": float(r[ac]) if ac else _twist_a, "b": float(r[bc]) if bc else _twist_b,
                                "theta": float(r[thc]) if thc else _twist_theta,
                                "label": f"clicked: A2={r[a2c]} Rt={r[rtc]} E={r[ec]:.2f}",
                            }
                    st.caption(f"{len(min_a2)} local minima marked (black squares) out of {len(df_conv)} grid points.")
                    with st.expander("Converged (a,b) at each (A2,Rt) grid point"):
                        st.dataframe(df_conv, width="stretch")
                with col_twist3d:
                    st.markdown("**Structure preview**")
                    if mol2 is None:
                        st.caption("Select a molecule using the picker above.")
                    else:
                        current_t = st.session_state.get("s4a_current")
                        if current_t is None and len(df_conv):
                            best_t = df_conv.loc[df_conv[ec].idxmin()]
                            current_t = {
                                "a2": float(best_t[a2c]), "rt": float(best_t[rtc]),
                                "a": float(best_t[ac]) if ac else _twist_a, "b": float(best_t[bc]) if bc else _twist_b,
                                "theta": float(best_t[thc]) if thc else _twist_theta,
                                "label": f"default (min E): A2={best_t[a2c]} Rt={best_t[rtc]}",
                            }
                            st.session_state["s4a_current"] = current_t
                        st.caption(current_t["label"])
                        c_i_t, c_j_t = dimer(
                            mol2, "t", current_t["a"], current_t["b"], current_t["theta"],
                            A2=current_t["a2"], z=current_t["rt"],
                        )
                        render_molecule_3d(
                            list(mol2.symbols) * 2, np.vstack([c_i_t, c_j_t]),
                            f"{mol2.name} twist dimer A2={current_t['a2']} Rt={current_t['rt']}",
                            key_suffix="s4a", style_key="s4a_style",
                        )


# ══════════════════════════════════════════════════════════
#  Tab 4: Step 3 – Interlayer Stacking  (para / twist sub-tabs)
# ══════════════════════════════════════════════════════════
with tab_step3:
    sub_para3, sub_twist3 = st.tabs(["para", "twist (→ Type III)"])

    with sub_para3:
        st.subheader("vdW pre-scan — runs in this GUI")
        st.caption(
            "Rigid-sphere interlayer contact model — a fast, vectorized "
            "reimplementation (not a call into legacy step3_para_vdw.py "
            "itself; verified to reproduce its output numerically) that "
            "runs directly in this GUI: at fixed intralayer parameters "
            "(a, b, alpha, Rt, Rp), slide the upper layer by (Ra, Rb) and find "
            "the vdW-limited interlayer distance z — this is the paper's "
            "Fig. 6(b–d) upper-panel V(x,y) = a·b·z map. Click a point (or one "
            "of the marked local minima) to preview the two-layer structure "
            "and download it as a starting point for the DFT step below."
        )
        mol3 = st.session_state.get("molecule")
        s1_current = st.session_state.get("s1vdw_current") or {}
        # Fallback defaults: pentacene's R-form (Type I) parameters from the
        # paper's SI Table S2 (a=7.2, b=6.0, alpha=25 deg), used when Tab 2
        # hasn't been run this session -- a real point instead of a made-up
        # placeholder.
        c1, c2, c3 = st.columns(3)
        s3_a = c1.number_input("a (Å)", value=float(s1_current.get("a", 7.2)), key="s3vdw_a")
        s3_b = c2.number_input("b (Å)", value=float(s1_current.get("b", 6.0)), key="s3vdw_b")
        s3_theta = c3.number_input(
            "alpha (deg, from Step 1)", value=float(s1_current.get("alpha", 25.0)), key="s3vdw_theta"
        )
        c4, c5 = st.columns(2)
        s3_rt = c4.number_input(
            "Rt (Å, T-contact long-axis shift)", value=0.0, step=0.1, key="s3vdw_rt"
        )
        s3_rp = c5.number_input(
            "Rp (Å, 0 = G-form direction, ≠0 = N-form direction)", value=0.0, step=0.1, key="s3vdw_rp"
        )

        if mol3 is None:
            st.info("Select a molecule using the picker above first.")
        else:
            # Auto-run once with the current inputs (Fig. 6(b)'s parameters
            # by default) so the map shows up without an extra click, same
            # as the bundled-sample sections elsewhere in this app -- but
            # only when nothing has been scanned yet this session, so it
            # never overrides a scan the user already ran or is about to
            # change inputs for.
            auto_run = (
                "s3vdw_df" not in st.session_state
                and mol3.name == "pentacene"
            )
            if auto_run or st.button(f"Run interlayer vdW scan ({mol3.name})", key="s3vdw_run"):
                with st.spinner("Scanning..."):
                    df_vdw = interlayer_vdw_scan(
                        mol3, a=s3_a, b=s3_b, theta=s3_theta, Rt=s3_rt, Rp=s3_rp,
                        radii_overrides=st.session_state.get("vdw_radii_overrides"),
                    )
                st.session_state["s3vdw_df"] = df_vdw
                st.session_state["s3vdw_params"] = {
                    "a": s3_a, "b": s3_b, "theta": s3_theta, "Rt": s3_rt, "Rp": s3_rp,
                }
                if auto_run:
                    st.caption(
                        "Showing the Fig. 6(b) default (pentacene R-form) "
                        "-- change the inputs above and click Run for your "
                        "own scan."
                    )

        df_vdw = st.session_state.get("s3vdw_df")
        params3 = st.session_state.get("s3vdw_params")
        if df_vdw is not None and "cz" not in df_vdw.columns:
            # Stale cache from before the cz column existed (session state
            # survives code hot-reloads) — drop it so the UI just asks for a
            # fresh scan instead of crashing on the missing column.
            for k in ("s3vdw_df", "s3vdw_params", "s3vdw_current"):
                st.session_state.pop(k, None)
            df_vdw = params3 = None
            st.info("The cached vdW scan was from an older version of this app — please run it again.")
        if df_vdw is not None and params3 is not None:
            if st.button("Clear vdW scan results", key="s3vdw_clear"):
                for k in ("s3vdw_df", "s3vdw_params", "s3vdw_current"):
                    st.session_state.pop(k, None)
                st.rerun()

            col_map3, col_3d3 = st.columns([2, 1])
            with col_map3:
                value_choice = st.radio(
                    "Color by", ["z (interlayer distance)", "V (unit cell volume)"],
                    horizontal=True, key="s3vdw_valchoice",
                )
                val_col = "z" if value_choice.startswith("z") else "V"
                pivot = df_vdw.pivot(index="Rb", columns="Ra", values=val_col)

                from scipy.ndimage import minimum_filter
                grid = pivot.values
                is_min = (grid == minimum_filter(grid, size=5)) & np.isfinite(grid)
                min_rb_idx, min_ra_idx = np.where(is_min)
                min_ra = pivot.columns.to_numpy()[min_ra_idx]
                min_rb = pivot.index.to_numpy()[min_rb_idx]

                # A plotly Heatmap/imshow trace doesn't reliably fire a
                # single-click "selection" event (that's really a box/lasso
                # select affordance for marker traces), which is why clicking
                # anywhere but the marker overlay never registered. Building
                # the whole map out of square markers instead — the same
                # mechanism already proven to work for Tab 2 — makes every
                # cell clickable.
                n_ra_cells = max(len(pivot.columns), 1)
                n_rb_cells = max(len(pivot.index), 1)
                # Ra and Rb use the same 0.1 Å step, so cells are square in
                # data units — but without locking the axes to a 1:1 scale
                # (below), Plotly renders Ra/Rb at different px-per-Å and a
                # single "square" marker size leaves gaps along one axis
                # (the visible black line). Pick one marker size from the
                # denser axis; scaleanchor makes both axes share that scale.
                marker_px = max(3, 650 / max(n_ra_cells, n_rb_cells))
                fig3 = go.Figure()
                fig3.add_trace(go.Scatter(
                    x=df_vdw["Ra"], y=df_vdw["Rb"], mode="markers", name=val_col,
                    showlegend=False,  # clicking a legend entry toggles trace
                    # visibility in Plotly — for this trace that hides the
                    # entire map, which just looks like it broke. The
                    # colorbar already conveys what this trace is.
                    marker=dict(
                        symbol="square", size=marker_px,
                        color=df_vdw[val_col], colorscale="RdBu_r",
                        colorbar=dict(title=val_col),
                    ),
                    hovertemplate="Ra=%{x}<br>Rb=%{y}<br>" + val_col + "=%{marker.color:.3f}<extra></extra>",
                ))
                fig3.add_trace(go.Scatter(
                    x=min_ra, y=min_rb, mode="markers", name="local min",
                    showlegend=False,  # explained in the caption below instead
                    marker=dict(symbol="square-open", size=marker_px + 6, color="gold", line=dict(width=2)),
                    hovertemplate="Ra=%{x}<br>Rb=%{y}<extra>local min</extra>",
                ))
                fig3.update_layout(
                    xaxis_title="Ra (Å, layer offset along a)",
                    yaxis_title="Rb (Å, layer offset along b)",
                    margin=dict(l=20, r=20, t=30, b=20),
                )
                # Lock Rb's px-per-Å scale to match Ra's, so the square
                # markers (equal size in px, equal 0.1 Å step in data) tile
                # without gaps regardless of the plot's rendered aspect ratio.
                fig3.update_yaxes(scaleanchor="x", scaleratio=1)
                event3 = st.plotly_chart(
                    fig3, width="stretch", on_select="rerun", key="s3vdw_fig"
                )
                pts3 = event3.selection.points if (event3 and event3.selection) else []
                if pts3:
                    # Both the heatmap and the local-min scatter overlay are
                    # built with x=Ra, y=Rb directly, so the point's plain x/y
                    # already is the (Ra, Rb) pair — no need for customdata
                    # (whose shape differs unpredictably between trace types
                    # and caused crashes both as a short list and as a dict).
                    p0 = pts3[0]
                    ra_sel, rb_sel = p0.get("x"), p0.get("y")
                    if ra_sel is not None and rb_sel is not None:
                        ident = (round(float(ra_sel), 1), round(float(rb_sel), 1))
                        if st.session_state.get("s3vdw_fig_prev") != ident:
                            row = df_vdw[
                                np.isclose(df_vdw["Ra"], ident[0]) & np.isclose(df_vdw["Rb"], ident[1])
                            ]
                            if len(row):
                                r = row.iloc[0]
                                st.session_state["s3vdw_current"] = {
                                    "cx": float(r["Ra"]), "cy": float(r["Rb"]), "cz": float(r["cz"]),
                                    "label": f"clicked: cx={r['Ra']} cy={r['Rb']} cz={r['cz']:.2f} V={r['V']:.1f}",
                                }
                                st.session_state["s3vdw_fig_prev"] = ident

                st.caption(f"{len(min_ra)} local minima marked (black squares) out of {len(df_vdw)} grid points.")

            with col_3d3:
                st.markdown("**Structure preview**")
                current3 = st.session_state.get("s3vdw_current")
                if current3 is None and len(df_vdw):
                    best = df_vdw.loc[df_vdw["z"].idxmin()]
                    current3 = {
                        "cx": float(best["Ra"]), "cy": float(best["Rb"]), "cz": float(best["cz"]),
                        "label": f"default (min z): cx={best['Ra']} cy={best['Rb']} cz={best['cz']:.2f}",
                    }
                if mol3 is not None and current3 is not None:
                    st.caption(current3["label"])
                    syms3, coords3 = bilayer_preview(
                        mol3, params3["a"], params3["b"], params3["theta"],
                        params3["Rt"], params3["Rp"],
                        current3["cx"], current3["cy"], current3["cz"],
                        radii_overrides=st.session_state.get("vdw_radii_overrides"),
                    )
                    render_molecule_3d(
                        syms3, coords3,
                        f"{mol3.name} bilayer cx={current3['cx']} cy={current3['cy']} cz={current3['cz']:.2f}",
                        key_suffix="s3vdw", style_key="s3vdw_style",
                    )
                    st.download_button(
                        "Download step3_para_init_params.csv (1 row for the DFT step below)",
                        data=pd.DataFrame([{
                            "a": params3["a"], "b": params3["b"], "theta": params3["theta"],
                            "Rt": params3["Rt"], "Rp": params3["Rp"],
                            "cx": current3["cx"], "cy": current3["cy"], "cz": current3["cz"],
                            "status": "NotYet",
                        }]).to_csv(index=False),
                        file_name="step3_para_init_params.csv", mime="text/csv",
                        key="dl_s3vdw_init",
                    )

        st.divider()
        st.subheader("How to run (CLI, i.e. commands typed in a terminal)")
        cli_howto(
            what=(
                "Step 3 (para variant) optimizes the interlayer c-vector "
                "(cx, cy, cz) by 3×3×3 hill-climbing at fixed intralayer "
                f"parameters (a, b, alpha, Rt, Rp) (`{LEGACY_DIR}/step3_para.py`). "
                "Each point computes 10 interlayer dimers; E averages the two "
                "stacking patterns and sums the four T-shaped pairs. "
                "**Rt/Rp are how this reaches the N-form** that Tab 3's "
                "uniform-z scan can't: the two T-shaped contacts sit at "
                "z = Rt and z = Rt−Rp, and the slipped-parallel contact at "
                "z = Rt−(Rt−Rp) = Rp. Rp = 0 keeps both T-contacts equal "
                "(glide symmetry / G-form direction, phi_incl = 0°/90° in "
                "the paper's §2.2 language); Rp ≠ 0 gives the two T-contacts "
                "different shifts — glide symmetry broken, the N-form "
                "direction (intermediate phi_incl)."
            ),
            prepare=(
                "- `step3_para_init_params.csv`: rows of `a, b, alpha, Rt, Rp, "
                "cx, cy, cz` (starting points). `a, b, alpha` come from Step 1 "
                "(**Tab 2**'s DFT results). `Rt, Rp` are **not** hill-climbed "
                "automatically (only cx, cy, cz are) — you supply the grid of "
                "(Rt, Rp) values to try yourself: Rp = 0 rows explore the "
                "G-form direction, Rp ≠ 0 rows explore the N-form direction. "
                "The initial `cx, cy, cz` guess for each row can come from the "
                "vdW pre-scan above (download the 1-row CSV after clicking a "
                "point) — concatenate rows for more than one starting point."
            ),
            setup=_SETUP_MONOMER_ENV + "\n" + _SETUP_SCHEDULER,
            command=(
                "python step3_para.py --auto-dir /path/to/workdir "
                "--monomer-name pentacene \\\n    --num-nodes 4 --num-init 2"
            ),
            output=(
                "`<auto-dir>/step3_para.csv` — columns `cx, cy, cz, a, b, alpha, "
                "Rt, Rp, E, E_i01, E_ip1, E_ip2, E_it1..4, E_i02, E_ip3, E_ip4, "
                "status, file_name`."
            ),
            scripts_files="`step3_para.py`, `make_step3_para.py`, and `utils.py`",
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
        st.subheader("vdW pre-scan — runs in this GUI")
        st.caption(
            "Same rigid-sphere interlayer contact model as Tab 4 para's vdW "
            "pre-scan, but with the twist torsion A2 applied to both "
            "intralayer clusters (Rp is implicitly 0 -- glide symmetric, "
            "matching step3_twist.py's a, b, alpha, Rt, A2 parameter set, "
            "no separate Rp). Click a point (or a marked local minimum) to "
            "preview the two-layer structure and download it as a starting "
            "point for the DFT step below."
        )
        mol3t = st.session_state.get("molecule")
        # Defaults are naphthalene's Fig. 6(c) G-form parameters (paper SI
        # Table S1, Type III calc row: a=5.8, b=7.5, alpha=65, ΔZT=1.6),
        # A2=0 since Fig. 6(c)'s V(x,y) map itself is computed *before* the
        # twist optimization introduced in Fig. 7/Sect. 2.4. Deliberately NOT
        # falling back to Tab 3 twist's s4a_current here (unlike Tab 4 para's
        # s1_current pattern) -- s4a_current is auto-populated the instant
        # Tab 3 twist's own sample map renders (its min-E point), so a
        # fallback to it would silently override these Fig. 6(c) defaults
        # with Tab 3's own sample values (a=6.0, b=7.2) every time -- numbers
        # easy to mistake for Fig. 6(b)'s pentacene (7.2, 6.0) at a glance.
        c1t, c2t, c3t = st.columns(3)
        s3t_a = c1t.number_input("a (Å)", value=5.8, key="s3twist_a")
        s3t_b = c2t.number_input("b (Å)", value=7.5, key="s3twist_b")
        s3t_theta = c3t.number_input(
            "alpha (deg, from Step 1)", value=65.0, key="s3twist_theta",
        )
        c4t, c5t = st.columns(2)
        s3t_rt = c4t.number_input(
            "Rt (Å, T-contact long-axis shift)", value=1.6, step=0.1, key="s3twist_rt",
        )
        s3t_a2 = c5t.number_input(
            "A2 (deg, twist torsion, from Tab 3 twist)", value=0.0, step=1.0, key="s3twist_a2",
        )

        if mol3t is None:
            st.info("Select a molecule using the picker above first.")
        else:
            # Auto-run once with the current inputs (Fig. 6(c)'s parameters
            # by default) so the map shows up without an extra click, same
            # as the bundled-sample sections elsewhere in this app -- but
            # only when nothing has been scanned yet this session, so it
            # never overrides a scan the user already ran or is about to
            # change inputs for.
            auto_run_t = (
                "s3twist_df" not in st.session_state
                and mol3t.name == "naphthalene"
            )
            if auto_run_t or st.button(f"Run interlayer vdW scan ({mol3t.name})", key="s3twist_run"):
                with st.spinner("Scanning..."):
                    df_vdw_t = interlayer_vdw_scan(
                        mol3t, a=s3t_a, b=s3t_b, theta=s3t_theta, Rt=s3t_rt, Rp=0.0, A2=s3t_a2,
                        radii_overrides=st.session_state.get("vdw_radii_overrides"),
                    )
                st.session_state["s3twist_df"] = df_vdw_t
                st.session_state["s3twist_params"] = {
                    "a": s3t_a, "b": s3t_b, "theta": s3t_theta, "Rt": s3t_rt, "A2": s3t_a2,
                }
                if auto_run_t:
                    st.caption(
                        "Showing the Fig. 6(c) default (naphthalene G-form) "
                        "-- change the inputs above and click Run for your "
                        "own scan."
                    )

        df_vdw_t = st.session_state.get("s3twist_df")
        params3t = st.session_state.get("s3twist_params")
        if df_vdw_t is not None and params3t is not None:
            if st.button("Clear vdW scan results", key="s3twist_clear"):
                for k in ("s3twist_df", "s3twist_params", "s3twist_current"):
                    st.session_state.pop(k, None)
                st.rerun()

            col_mapt, col_3dt = st.columns([2, 1])
            with col_mapt:
                value_choice_t = st.radio(
                    "Color by", ["z (interlayer distance)", "V (unit cell volume)"],
                    horizontal=True, key="s3twist_valchoice",
                )
                val_col_t = "z" if value_choice_t.startswith("z") else "V"
                pivot_vt = df_vdw_t.pivot(index="Rb", columns="Ra", values=val_col_t)

                from scipy.ndimage import minimum_filter
                grid_vt = pivot_vt.values
                is_min_vt = (grid_vt == minimum_filter(grid_vt, size=5)) & np.isfinite(grid_vt)
                min_rb_idx_t, min_ra_idx_t = np.where(is_min_vt)
                min_ra_t = pivot_vt.columns.to_numpy()[min_ra_idx_t]
                min_rb_t = pivot_vt.index.to_numpy()[min_rb_idx_t]

                n_ra_cells_t = max(len(pivot_vt.columns), 1)
                n_rb_cells_t = max(len(pivot_vt.index), 1)
                marker_px_vt = max(3, 650 / max(n_ra_cells_t, n_rb_cells_t))
                figVT = go.Figure()
                figVT.add_trace(go.Scatter(
                    x=df_vdw_t["Ra"], y=df_vdw_t["Rb"], mode="markers", name=val_col_t,
                    showlegend=False,
                    marker=dict(
                        symbol="square", size=marker_px_vt,
                        color=df_vdw_t[val_col_t], colorscale="RdBu_r",
                        colorbar=dict(title=val_col_t),
                    ),
                    hovertemplate="Ra=%{x}<br>Rb=%{y}<br>" + val_col_t + "=%{marker.color:.3f}<extra></extra>",
                ))
                figVT.add_trace(go.Scatter(
                    x=min_ra_t, y=min_rb_t, mode="markers", name="local min",
                    showlegend=False,
                    marker=dict(symbol="square-open", size=marker_px_vt + 6, color="gold", line=dict(width=2)),
                    hovertemplate="Ra=%{x}<br>Rb=%{y}<extra>local min</extra>",
                ))
                figVT.update_layout(
                    xaxis_title="Ra (Å, layer offset along a)",
                    yaxis_title="Rb (Å, layer offset along b)",
                    margin=dict(l=20, r=20, t=30, b=20),
                )
                figVT.update_yaxes(scaleanchor="x", scaleratio=1)
                eventVT = st.plotly_chart(
                    figVT, width="stretch", on_select="rerun", key="s3twist_fig"
                )
                ptsVT = eventVT.selection.points if (eventVT and eventVT.selection) else []
                if ptsVT:
                    p0t = ptsVT[0]
                    ra_sel_t, rb_sel_t = p0t.get("x"), p0t.get("y")
                    if ra_sel_t is not None and rb_sel_t is not None:
                        ident_t = (round(float(ra_sel_t), 1), round(float(rb_sel_t), 1))
                        if st.session_state.get("s3twist_fig_prev") != ident_t:
                            row_t = df_vdw_t[
                                np.isclose(df_vdw_t["Ra"], ident_t[0]) & np.isclose(df_vdw_t["Rb"], ident_t[1])
                            ]
                            if len(row_t):
                                rr = row_t.iloc[0]
                                st.session_state["s3twist_current"] = {
                                    "cx": float(rr["Ra"]), "cy": float(rr["Rb"]), "cz": float(rr["cz"]),
                                    "label": f"clicked: cx={rr['Ra']} cy={rr['Rb']} cz={rr['cz']:.2f} V={rr['V']:.1f}",
                                }
                                st.session_state["s3twist_fig_prev"] = ident_t

                st.caption(f"{len(min_ra_t)} local minima marked (gold squares) out of {len(df_vdw_t)} grid points.")

            with col_3dt:
                st.markdown("**Structure preview**")
                current3t = st.session_state.get("s3twist_current")
                if current3t is None and len(df_vdw_t):
                    best_t3 = df_vdw_t.loc[df_vdw_t["z"].idxmin()]
                    current3t = {
                        "cx": float(best_t3["Ra"]), "cy": float(best_t3["Rb"]), "cz": float(best_t3["cz"]),
                        "label": f"default (min z): cx={best_t3['Ra']} cy={best_t3['Rb']} cz={best_t3['cz']:.2f}",
                    }
                if mol3t is not None and current3t is not None:
                    st.caption(current3t["label"])
                    syms3t, coords3t = bilayer_preview(
                        mol3t, params3t["a"], params3t["b"], params3t["theta"],
                        params3t["Rt"], 0.0,
                        current3t["cx"], current3t["cy"], current3t["cz"],
                        A2=params3t["A2"],
                        radii_overrides=st.session_state.get("vdw_radii_overrides"),
                    )
                    render_molecule_3d(
                        syms3t, coords3t,
                        f"{mol3t.name} twisted bilayer cx={current3t['cx']} "
                        f"cy={current3t['cy']} cz={current3t['cz']:.2f}",
                        key_suffix="s3twist", style_key="s3twist_style",
                    )
                    st.download_button(
                        "Download step3_twist_init_params.csv (1 row for the DFT step below)",
                        data=pd.DataFrame([{
                            "a": params3t["a"], "b": params3t["b"], "theta": params3t["theta"],
                            "Rt": params3t["Rt"], "A2": params3t["A2"],
                            "cx": current3t["cx"], "cy": current3t["cy"], "cz": current3t["cz"],
                            "status": "NotYet",
                        }]).to_csv(index=False),
                        file_name="step3_twist_init_params.csv", mime="text/csv",
                        key="dl_s3twist_init",
                    )

        st.divider()
        st.subheader("How to run (CLI, i.e. commands typed in a terminal)")
        cli_howto(
            what=(
                "Step 3 (twist variant): interlayer c-vector optimization "
                f"for the twisted (Type III) packing (`{LEGACY_DIR}/step3_twist.py`)."
            ),
            prepare=(
                "- `step3_twist_init_params.csv`: rows of `a, b, alpha, Rt, "
                "A2, cx, cy, cz` (starting points). `a, b, alpha, Rt, A2` "
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
            scripts_files="`step3_twist.py`, `make_step3_twist.py`, and `utils.py`",
        )
        st.divider()
        st.subheader("Fig. 7(d)-style: Total E vs theta_twist (A2)")
        st.caption(
            "The paper's actual Fig. 7(d): for each theta_twist (A2) in the "
            "promising range narrowed down from Tab 3 twist's Fig. 7(c) map, "
            "(a, b, Rt) is re-optimized by hill-climbing (a separate, denser "
            "CSV from Tab 3 twist's one-point-per-grid-cell CSV), then the "
            "interlayer energy Eint(near) is computed on top of that "
            "structure. Total E = Eintra(6) + Eint(near), taking the minimum "
            "of each at every A2. **Needs two CSVs, one in each box below:** "
            "(1) the re-optimized intralayer grid -- columns `A2, Rt, a, b, "
            "E` (alpha optional); (2) the interlayer scan, i.e. "
            "step3_twist.csv -- columns `A2, cx, cy, cz, E`."
        )
        df_intra_d, src_intra_d = load_results_csv(
            "s4d_results", "step2_twist_refined.csv",
            label="1) Re-optimized intralayer grid (A2, Rt, a, b, E) -- "
                  "sample shown until you drop your own",
        )
        df, src = load_results_csv(
            "s4b_results", "step3_twist.csv",
            label="2) Interlayer scan, step3_twist.csv (A2, cx, cy, cz, E) -- "
                  "sample shown until you drop your own",
        )
        if df_intra_d is not None and df is not None:
            if src == "sample" or src_intra_d == "sample":
                st.info(
                    "This sample is Ono's real naphthalene data (both CSVs "
                    "share the same re-optimized (A2, Rt, a, b) per point, "
                    "confirmed by joining them here). Select naphthalene in "
                    "the picker above for the 3D preview to match the "
                    "plotted lattice constants."
                )
            cl_id = {c.lower(): c for c in df_intra_d.columns}
            cl_er = {c.lower(): c for c in df.columns}
            need_id = {"a2", "e", "a", "b", "rt"}
            need_er = {"a2", "e"}
            if not (need_id <= set(cl_id) and need_er <= set(cl_er)):
                st.warning(
                    f"Expected A2, Rt, a, b, E in the intralayer CSV (got "
                    f"{list(df_intra_d.columns)}) and A2, E in step3_twist.csv "
                    f"(got {list(df.columns)})."
                )
            else:
                a2_id, e_id, a_id, b_id, rt_id = (
                    cl_id["a2"], cl_id["e"], cl_id["a"], cl_id["b"], cl_id["rt"]
                )
                a2_er, e_er = cl_er["a2"], cl_er["e"]
                best_intra = df_intra_d.loc[df_intra_d.groupby(a2_id)[e_id].idxmin()].set_index(a2_id)
                best_inter = df.loc[df.groupby(a2_er)[e_er].idxmin()].set_index(a2_er)[e_er]
                common_a2 = sorted(set(best_intra.index) & set(best_inter.index))
                if not common_a2:
                    st.warning("No shared A2 values between the two CSVs -- can't join.")
                else:
                    df_total = pd.DataFrame({
                        "A2": common_a2,
                        "Eintra": [float(best_intra.loc[a2, e_id]) for a2 in common_a2],
                        "Einter": [float(best_inter.loc[a2]) for a2 in common_a2],
                        "a": [float(best_intra.loc[a2, a_id]) for a2 in common_a2],
                        "b": [float(best_intra.loc[a2, b_id]) for a2 in common_a2],
                        "Rt": [float(best_intra.loc[a2, rt_id]) for a2 in common_a2],
                    })
                    thc_d = cl_id.get("alpha", cl_id.get("theta"))
                    df_total["theta"] = (
                        [float(best_intra.loc[a2, thc_d]) for a2 in common_a2] if thc_d
                        else _twist_theta
                    )
                    df_total["Total"] = df_total["Eintra"] + df_total["Einter"]
                    # Plot every curve relative to its own A2=0 (untwisted)
                    # value -- the paper's Fig. 7(d) shows stabilization
                    # *from* the untwisted reference, not absolute E, and the
                    # three quantities' absolute offsets aren't otherwise
                    # comparable on one axis.
                    ref_row = df_total[np.isclose(df_total["A2"], 0.0)]
                    if len(ref_row):
                        ref = ref_row.iloc[0]
                    else:
                        ref = df_total.iloc[(df_total["A2"] - 0.0).abs().argmin()]
                    for col in ("Total", "Eintra", "Einter"):
                        df_total[f"{col}_rel"] = df_total[col] - ref[col]
                    col_d_plot, col_d_3d = st.columns([2, 1])
                    with col_d_plot:
                        figD = go.Figure()
                        figD.add_trace(go.Scatter(
                            x=df_total["A2"], y=df_total["Total_rel"], mode="lines+markers",
                            name="Total (Eintra+Einter)", line=dict(color="gold", width=3),
                        ))
                        figD.add_trace(go.Scatter(
                            x=df_total["A2"], y=df_total["Eintra_rel"], mode="lines+markers",
                            name="Eintra(6)", line=dict(dash="dot"),
                        ))
                        figD.add_trace(go.Scatter(
                            x=df_total["A2"], y=df_total["Einter_rel"], mode="lines+markers",
                            name="Eint(near)", line=dict(dash="dash"),
                        ))
                        figD.update_layout(
                            xaxis_title="theta_twist / A2 (deg)",
                            yaxis_title="E relative to A2=0 (kcal/mol)",
                            margin=dict(l=20, r=20, t=30, b=20),
                        )
                        event_d = st.plotly_chart(figD, width="stretch", on_select="rerun", key="s4d_chart")
                        pts_d = event_d.selection.points if (event_d and event_d.selection) else []
                        if pts_d:
                            x_sel_d = pts_d[0].get("x")
                            if x_sel_d is not None:
                                match_d = df_total[np.isclose(df_total["A2"], x_sel_d)]
                                if len(match_d):
                                    r = match_d.iloc[0]
                                    st.session_state["s4d_current"] = {
                                        "a2": float(r["A2"]), "rt": float(r["Rt"]),
                                        "a": float(r["a"]), "b": float(r["b"]),
                                        "theta": float(r["theta"]),
                                        "total": float(r["Total_rel"]), "eintra": float(r["Eintra_rel"]),
                                        "einter": float(r["Einter_rel"]),
                                        "label": f"clicked: A2={r['A2']} Rt={r['Rt']}",
                                    }
                        min_row_d = df_total.loc[df_total["Total_rel"].idxmin()]
                        st.caption(
                            f"Minimum Total E at A2={min_row_d['A2']:.0f} deg "
                            f"(ΔTotal={min_row_d['Total_rel']:.2f} kcal/mol "
                            "relative to A2=0)."
                        )
                        with st.expander("Data table"):
                            st.dataframe(df_total, width="stretch")
                    with col_d_3d:
                        st.markdown("**Structure preview**")
                        if mol2 is None:
                            st.caption("Select a molecule using the picker above.")
                        else:
                            current_d = st.session_state.get("s4d_current")
                            if current_d is None:
                                best_d = df_total.loc[df_total["Total_rel"].idxmin()]
                                current_d = {
                                    "a2": float(best_d["A2"]), "rt": float(best_d["Rt"]),
                                    "a": float(best_d["a"]), "b": float(best_d["b"]),
                                    "theta": float(best_d["theta"]),
                                    "total": float(best_d["Total_rel"]), "eintra": float(best_d["Eintra_rel"]),
                                    "einter": float(best_d["Einter_rel"]),
                                    "label": f"default (min Total E): A2={best_d['A2']} Rt={best_d['Rt']}",
                                }
                                st.session_state["s4d_current"] = current_d
                            st.caption(current_d["label"])
                            st.metric("ΔTotal E vs A2=0 (kcal/mol)", f"{current_d['total']:.2f}")
                            st.caption(
                                f"ΔEintra(6)={current_d['eintra']:.2f}, "
                                f"ΔEint(near)={current_d['einter']:.2f} "
                                "(all relative to A2=0)"
                            )
                            c_i_d, c_j_d = dimer(
                                mol2, "t", current_d["a"], current_d["b"], current_d["theta"],
                                A2=current_d["a2"], z=current_d["rt"],
                            )
                            render_molecule_3d(
                                list(mol2.symbols) * 2, np.vstack([c_i_d, c_j_d]),
                                f"{mol2.name} twist dimer A2={current_d['a2']} Rt={current_d['rt']}",
                                key_suffix="s4d", style_key="s4d_style",
                            )

        st.divider()
        st.subheader("Interlayer scan detail  (step3_twist.csv)")
        if df is not None:
            source_badge(src)
            cols_lower_detail = {c.lower(): c for c in df.columns}
            a2_detail = cols_lower_detail.get("a2")
            e_detail = cols_lower_detail.get("e")
            if a2_detail:
                a2_values_detail = sorted(df[a2_detail].unique())
                default_a2_detail = (
                    df.loc[df[e_detail].idxmin(), a2_detail] if e_detail else a2_values_detail[0]
                )
                sel_a2_detail = st.selectbox(
                    "theta_twist / A2 (deg) -- each A2 has its own hill-climb "
                    "search, so pick one to view its own (cy, cz) landscape",
                    a2_values_detail,
                    index=(
                        a2_values_detail.index(default_a2_detail)
                        if default_a2_detail in a2_values_detail else 0
                    ),
                    key="s4b_detail_a2",
                )
                df_detail = df[df[a2_detail] == sel_a2_detail]
                st.caption(
                    f"The (cy, cz) interlayer c-vector hill-climb search at "
                    f"this single A2={sel_a2_detail} (cx is fixed at 0 for "
                    "this twisted packing, so it isn't a useful map axis). "
                    "Each A2 explores a different (cy, cz) region since "
                    "(a, b, Rt) also change with A2 -- overlaying every A2 "
                    "on one map would mix unrelated hill-climb trajectories "
                    "together."
                )
            else:
                df_detail = df
                st.caption(
                    "The raw (cy, cz) interlayer c-vector scan (cx is fixed "
                    "at 0 for this twisted packing, so it isn't a useful map "
                    "axis)."
                )
            heatmap_section(
                df_detail, "s4b_map",
                x_candidates=["cy", "x"],
                y_candidates=["cz", "y"],
                val_candidates=["e"],
                val_label="E (kcal/mol)",
            )


# ══════════════════════════════════════════════════════════
#  Tab 5: Transfer Integrals  (CLI; results shown here)
# ══════════════════════════════════════════════════════════
with tab_transfer:
    st.subheader("How to run (CLI, i.e. commands typed in a terminal)")
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
            "- `<auto-dir>/init_params.csv`: rows of `a, b, alpha, A2, z` — "
            "the representative arrangements (e.g. the R-form optimum from "
            "**Tab 2**'s DFT results, and any Type II/III/IV structures "
            "from **Tab 4**'s para/twist sub-tabs) you want J for."
        ),
        setup=(
            _SETUP_MONOMER_ENV.replace("Tab 1", "Tab 1 — this feeds `get_monomer_xyzR` in tcal_csv.py too") + "\n"
            + "- `job.sh` (next to `tcal_csv.py`) is a Fujitsu/PJM batch script "
            "template as uploaded; edit it (or `qsub_process()`'s `pjsub "
            "job.sh` call) to match whatever scheduler *your* cluster uses."
        ),
        command=(
            "python tcal_csv.py --init   --auto-dir /path/to/workdir --monomer-name pentacene   # build inputs\n"
            "python tcal_csv.py --qsub   --auto-dir /path/to/workdir --monomer-name pentacene   # submit Gaussian\n"
            "python tcal_csv.py --tcal   --auto-dir /path/to/workdir --monomer-name pentacene   # run tcal_1.py\n"
            "python tcal_csv.py --result --auto-dir /path/to/workdir --monomer-name pentacene   # collect result.txt"
        ),
        output="`<auto-dir>/result.txt` — space-separated HOMO transfer integrals per row (T-shaped, then slipped-parallel).",
        show_route=False,
        tree=_TREE_TCAL,
        scripts_dir="legacy/ono_scripts/tcal_csv",
        scripts_files="`tcal_csv.py`, `tcal_1.py`, `utils.py`, and `job.sh`",
    )

    st.divider()
    st.subheader("Results — J per arrangement / contact type")
    st.caption(
        "Drop a CSV (e.g. result.txt converted to CSV with columns like "
        "a, b, alpha, J_t, J_p)."
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

    st.divider()
    st.subheader("Fig. 11(b)-style: J vs alpha (HB / PS branches)")
    st.caption(
        "Transfer integrals along the Fig. 2(b) branches (A2=0, z=0 at "
        "each alpha): run tcal_csv.py with one init_params.csv row per "
        "(a, b, alpha) point from Tab 2's HB/PS branches, "
        "then combine init_params.csv with result.txt's J columns (same "
        "row order) into one CSV. Expected columns: alpha (or theta), J_t, "
        "J_p, and optionally kind (e.g. HB/PS) to color "
        "separate branches, and E (E_intra(8), kcal/mol) for the upper "
        "total-energy panel -- plotted on its own axis since its scale "
        "(kcal/mol) isn't comparable to J (meV). Click a J_t/J_p point "
        "for its 2-molecule dimer, or a Total point for the full "
        "9-molecule cluster; the clicked point's own value is always "
        "shown as text next to the structure."
    )
    df_b, src_b = load_results_csv("s5b_results", "transfer_integrals_alpha.csv")
    if df_b is not None:
        source_badge(src_b)
        cols_lower_b = {c.lower(): c for c in df_b.columns}
        xcand_b = [c for c in ("alpha", "theta") if c in cols_lower_b]
        if not xcand_b or not {"j_t", "j_p"} <= set(cols_lower_b):
            st.warning(f"Expected columns alpha/theta, J_t, J_p (got {list(df_b.columns)}).")
        else:
            xb = cols_lower_b[xcand_b[0]]
            jtb, jpb = cols_lower_b["j_t"], cols_lower_b["j_p"]
            kindb = cols_lower_b.get("kind")
            ab_ = cols_lower_b.get("a")
            bb_ = cols_lower_b.get("b")
            eb_ = cols_lower_b.get("e")
            col_b_plot, col_b_3d = st.columns([2, 1])
            with col_b_plot:
                has_total = eb_ is not None
                if has_total:
                    figb = make_subplots(
                        rows=2, cols=1, shared_xaxes=True, row_heights=[0.3, 0.7],
                        vertical_spacing=0.06,
                    )
                else:
                    figb = go.Figure()
                trace_info_b = []  # ("total"/"t"/"p", kind), aligned with trace order
                groups_b = list(df_b.groupby(kindb)) if kindb else [(None, df_b)]
                for kind, grp in groups_b:
                    grp = grp.sort_values(xb)
                    suffix = f" ({kind})" if kind is not None else ""
                    if has_total:
                        figb.add_trace(go.Scatter(
                            x=grp[xb], y=grp[eb_], mode="lines+markers", name=f"Total{suffix}",
                            line=dict(dash="dot"),
                        ), row=1, col=1)
                        trace_info_b.append(("total", kind))
                        figb.add_trace(go.Scatter(
                            x=grp[xb], y=grp[jtb], mode="lines+markers", name=f"J_t{suffix}",
                        ), row=2, col=1)
                        trace_info_b.append(("t", kind))
                        figb.add_trace(go.Scatter(
                            x=grp[xb], y=grp[jpb], mode="lines+markers", name=f"J_p{suffix}",
                            line=dict(dash="dash"),
                        ), row=2, col=1)
                        trace_info_b.append(("p", kind))
                    else:
                        figb.add_trace(go.Scatter(
                            x=grp[xb], y=grp[jtb], mode="lines+markers", name=f"J_t{suffix}",
                        ))
                        trace_info_b.append(("t", kind))
                        figb.add_trace(go.Scatter(
                            x=grp[xb], y=grp[jpb], mode="lines+markers", name=f"J_p{suffix}",
                            line=dict(dash="dash"),
                        ))
                        trace_info_b.append(("p", kind))
                if has_total:
                    figb.update_yaxes(title_text="E_intra(8) (kcal/mol)", row=1, col=1)
                    figb.update_yaxes(title_text="J (meV)", row=2, col=1)
                    figb.update_xaxes(title_text=xb, row=2, col=1)
                    figb.update_layout(margin=dict(l=20, r=20, t=30, b=20))
                else:
                    figb.update_layout(
                        xaxis_title=xb, yaxis_title="J (meV)",
                        margin=dict(l=20, r=20, t=30, b=20),
                    )
                event_b = st.plotly_chart(figb, width="stretch", on_select="rerun", key="s5b_chart")
                pts_b = event_b.selection.points if (event_b and event_b.selection) else []
                if pts_b and ab_ and bb_:
                    p0b = pts_b[0]
                    x_sel_b, curve_num_b = p0b.get("x"), p0b.get("curve_number", 0)
                    if x_sel_b is not None and 0 <= curve_num_b < len(trace_info_b):
                        contact_sel, kind_sel = trace_info_b[curve_num_b]
                        mask_b = np.isclose(df_b[xb], x_sel_b)
                        if kindb:
                            mask_b &= (df_b[kindb] == kind_sel)
                        match_b = df_b[mask_b]
                        if len(match_b):
                            r = match_b.iloc[0]
                            a_r, b_r = float(r[ab_]), float(r[bb_])
                            kind_label = f" ({kind_sel})" if kind_sel is not None else ""
                            if contact_sel == "total":
                                st.session_state["s5b_current"] = {
                                    "kind": "total", "a": a_r, "b": b_r,
                                    "theta": float(r[xb]),
                                    "value": float(r[eb_]) if eb_ else None,
                                    "value_label": "E_intra(8) (kcal/mol)",
                                    "label": f"clicked: Total{kind_label} at alpha={r[xb]}",
                                }
                            else:
                                # tcal_csv.py's own rule: the slipped-parallel
                                # contact shifts along whichever of a/b is
                                # shorter (p2 along a when b>a, else p1 along b).
                                dimer_kind = "t" if contact_sel == "t" else ("p2" if b_r > a_r else "p1")
                                jval = float(r[jtb]) if contact_sel == "t" else float(r[jpb])
                                st.session_state["s5b_current"] = {
                                    "kind": dimer_kind, "a": a_r, "b": b_r,
                                    "theta": float(r[xb]),
                                    "value": jval, "value_label": "J (meV)",
                                    "label": (
                                        f"clicked: {'J_t' if contact_sel == 't' else 'J_p'}"
                                        f"{kind_label} at alpha={r[xb]}"
                                    ),
                                }
                with st.expander("Data table"):
                    st.dataframe(df_b, width="stretch")
            with col_b_3d:
                st.markdown("**Structure preview**")
                if mol2 is None:
                    st.caption("Select a molecule using the picker above.")
                elif not (ab_ and bb_):
                    st.caption("Need `a`, `b` columns in the CSV for a 3D preview.")
                else:
                    current_b = st.session_state.get("s5b_current")
                    if current_b is None and len(df_b):
                        first_b = df_b.iloc[0]
                        kind_label0 = f" ({first_b[kindb]})" if kindb else ""
                        current_b = {
                            "kind": "t", "a": float(first_b[ab_]), "b": float(first_b[bb_]),
                            "theta": float(first_b[xb]),
                            "value": float(first_b[jtb]), "value_label": "J (meV)",
                            "label": f"default: J_t{kind_label0} at alpha={first_b[xb]}",
                        }
                        st.session_state["s5b_current"] = current_b
                    st.caption(current_b["label"])
                    if current_b.get("value") is not None:
                        st.metric(current_b["value_label"], f"{current_b['value']:.2f}")
                    if current_b["kind"] == "total":
                        syms_b, coords_b = cluster9(mol2, current_b["a"], current_b["b"], current_b["theta"])
                        render_molecule_3d(
                            syms_b, coords_b,
                            f"{mol2.name} cluster9 alpha={current_b['theta']}",
                            key_suffix="s5b", style_key="s5b_style",
                        )
                    else:
                        c_i_b, c_j_b = dimer(
                            mol2, current_b["kind"], current_b["a"], current_b["b"], current_b["theta"],
                        )
                        render_molecule_3d(
                            list(mol2.symbols) * 2, np.vstack([c_i_b, c_j_b]),
                            f"{mol2.name} {current_b['kind']} dimer alpha={current_b['theta']}",
                            key_suffix="s5b", style_key="s5b_style",
                        )

    st.divider()
    st.subheader("Fig. 11(c)-style: J vs theta_incl")
    st.caption(
        "Transfer integrals along the uniform long-axis inclination "
        "direction (Step 2 para's z / theta_incl, glide-symmetric so a "
        "single z per row is enough — the same z used in Tab 3's "
        "Et(z)/Ep(z) scan feeds tcal_csv.py's init_params.csv `z` "
        "column). Expected columns: theta_incl (or z), J_t, J_p, and "
        "optionally a, b, alpha (fixed at the Step 1 optimum) for the "
        "3D structure preview on click."
    )
    df_c, src_c = load_results_csv("s5c_results", "transfer_integrals_thetaincl.csv")
    if df_c is not None:
        source_badge(src_c)
        cols_lower_c = {c.lower(): c for c in df_c.columns}
        xcand_c = [c for c in ("theta_incl", "z") if c in cols_lower_c]
        if not xcand_c or not {"j_t", "j_p"} <= set(cols_lower_c):
            st.warning(f"Expected columns theta_incl/z, J_t, J_p (got {list(df_c.columns)}).")
        else:
            xc_ = cols_lower_c[xcand_c[0]]
            jtc_, jpc_ = cols_lower_c["j_t"], cols_lower_c["j_p"]
            ac_ = cols_lower_c.get("a")
            bc_ = cols_lower_c.get("b")
            thc_ = cols_lower_c.get("alpha", cols_lower_c.get("theta"))
            grp_c = df_c.sort_values(xc_)
            col_c_plot, col_c_3d = st.columns([2, 1])
            with col_c_plot:
                figc = go.Figure()
                figc.add_trace(go.Scatter(
                    x=grp_c[xc_], y=grp_c[jtc_], mode="lines+markers", name="J_t (T-shaped)",
                ))
                figc.add_trace(go.Scatter(
                    x=grp_c[xc_], y=grp_c[jpc_], mode="lines+markers", name="J_p (slipped-parallel)",
                ))
                figc.update_layout(
                    xaxis_title=xc_, yaxis_title="J (meV)",
                    margin=dict(l=20, r=20, t=30, b=20),
                )
                event_c = st.plotly_chart(figc, width="stretch", on_select="rerun", key="s5c_chart")
                pts_c = event_c.selection.points if (event_c and event_c.selection) else []
                if pts_c and ac_ and bc_ and thc_:
                    p0c = pts_c[0]
                    x_sel_c, curve_num_c = p0c.get("x"), p0c.get("curve_number", 0)
                    if x_sel_c is not None:
                        match_c = df_c[np.isclose(df_c[xc_], x_sel_c)]
                        if len(match_c):
                            r = match_c.iloc[0]
                            a_r, b_r = float(r[ac_]), float(r[bc_])
                            contact_sel = "t" if curve_num_c == 0 else "p"
                            dimer_kind = "t" if contact_sel == "t" else ("p2" if b_r > a_r else "p1")
                            st.session_state["s5c_current"] = {
                                "kind": dimer_kind, "a": a_r, "b": b_r,
                                "theta": float(r[thc_]), "z": float(x_sel_c),
                                "label": (
                                    f"clicked: {'J_t' if contact_sel == 't' else 'J_p'} "
                                    f"at {xc_}={x_sel_c}"
                                ),
                            }
                with st.expander("Data table"):
                    st.dataframe(df_c, width="stretch")
            with col_c_3d:
                st.markdown("**Structure preview**")
                if mol2 is None:
                    st.caption("Select a molecule using the picker above.")
                elif not (ac_ and bc_ and thc_):
                    st.caption("Need `a`, `b`, `alpha` columns in the CSV for a 3D preview.")
                else:
                    current_c = st.session_state.get("s5c_current")
                    if current_c is None and len(grp_c):
                        first_c = grp_c.iloc[0]
                        current_c = {
                            "kind": "t", "a": float(first_c[ac_]), "b": float(first_c[bc_]),
                            "theta": float(first_c[thc_]), "z": float(first_c[xc_]),
                            "label": f"default: J_t at {xc_}={first_c[xc_]}",
                        }
                        st.session_state["s5c_current"] = current_c
                    st.caption(current_c["label"])
                    c_i_c, c_j_c = dimer(
                        mol2, current_c["kind"], current_c["a"], current_c["b"], current_c["theta"],
                        z=current_c["z"],
                    )
                    render_molecule_3d(
                        list(mol2.symbols) * 2, np.vstack([c_i_c, c_j_c]),
                        f"{mol2.name} {current_c['kind']} dimer {xc_}={current_c['z']}",
                        key_suffix="s5c", style_key="s5c_style",
                    )
