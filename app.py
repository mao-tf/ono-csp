#!/usr/bin/env python3
"""csp — Crystal Structure Prediction UI (Streamlit)

Usage:
  streamlit run app.py

Reproduces the structure-search pipeline of Ono et al. (JACS, submitted).
See spec.md for the full method-to-tab mapping. Tab 1 (Molecule Setup) is
functional; Tabs 2-7 are UI skeletons for the calculation engines described
in spec.md, which are not yet implemented (see each tab's caption for what
is pending and why).
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from csp.structure.molecule import load_molecule, _VDW  # noqa: E402
from csp.plot.viewer3d import to_xyz_string, render_3d_html  # noqa: E402

MOLECULE_DIR = Path(__file__).resolve().parent / "data" / "molecules"
PRESET_MOLECULES = ["naphthalene", "anthracene", "tetracene", "pentacene", "hexacene"]

st.set_page_config(page_title="csp — Crystal Structure Prediction", layout="wide")
st.title("csp — Crystal Structure Prediction")
st.caption(
    "Reproduces Ono et al., \"Origin of Layered Herringbone Packing and "
    "Polymorphism in Polyacenes\" (JACS, submitted). See spec.md for the "
    "full method-to-tab mapping."
)

tab_setup, tab_step1_vdw, tab_step1_dft, tab_step2, tab_step3, tab_step4, tab_step5 = st.tabs([
    "1. Molecule Setup",
    "2. Step 1 – Intralayer vdW Scan",
    "3. Step 1 – DFT-D Optimization",
    "4. Step 2 – Inclination Map",
    "5. Step 3 – Interlayer Stacking",
    "6. Step 4 – Refinement",
    "7. Transfer Integrals",
])


def _not_implemented(step: str, spec_section: str) -> None:
    st.info(
        f"{step} calculation engine is not yet implemented. "
        f"See spec.md {spec_section} for the algorithm, or "
        f"legacy/ono_scripts/ once Ono's existing code is integrated."
    )


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
#  Tab 2: Step 1 – Intralayer vdW Scan
# ══════════════════════════════════════════════════════════
with tab_step1_vdw:
    st.subheader("Settings")
    c1, c2, c3 = st.columns(3)
    alpha_min = c1.number_input("alpha min (deg)", value=5.0, key="s1vdw_amin")
    alpha_max = c2.number_input("alpha max (deg)", value=85.0, key="s1vdw_amax")
    alpha_step = c3.number_input("alpha step (deg)", value=5.0, min_value=0.1, key="s1vdw_astep")
    vdw_tol = st.number_input("VdW tolerance (Å, overlap threshold)", value=0.0, key="s1vdw_tol")

    if st.button("Run vdW Scan", key="s1vdw_run"):
        _not_implemented("Step 1a (vdW scan)", "§Step 1a")

    st.divider()
    col_l, col_r = st.columns(2)
    with col_l:
        st.caption("alpha vs S=a×b plot (Plotly) — click a point to preview its structure")
    with col_r:
        st.caption("Layer structure preview (2D) — T-shaped / SP contacts")


# ══════════════════════════════════════════════════════════
#  Tab 3: Step 1 – DFT-D Optimization
# ══════════════════════════════════════════════════════════
with tab_step1_dft:
    st.subheader("HPC settings")
    c1, c2 = st.columns(2)
    hpc_host = c1.text_input("Host", value="miyoshi@133.11.68.31", key="s1dft_host")
    hpc_workdir = c2.text_input(
        "Work dir", value="/home/miyoshi/Working/ono_paper_dir", key="s1dft_workdir"
    )

    st.subheader("Gaussian keywords")
    st.code("# B3LYP empiricaldispersion=gd3bj 6-311g** counterpoise=2 nosymm", language="text")

    c3, c4 = st.columns(2)
    if c3.button("Generate & Submit Jobs", key="s1dft_submit"):
        _not_implemented("Step 1b (DFT-D refinement)", "§Step 1b")
    if c4.button("Check Status", key="s1dft_status"):
        _not_implemented("Step 1b job status check", "§Step 1b")

    st.divider()
    col_l, col_r = st.columns(2)
    with col_l:
        st.caption("alpha vs E_intra(8) plot (kcal/mol) — click to preview structure")
    with col_r:
        st.caption("3D structure viewer — R-form highlighted (alpha ≈ 25°)")


# ══════════════════════════════════════════════════════════
#  Tab 4: Step 2 – Inclination Map
# ══════════════════════════════════════════════════════════
with tab_step2:
    st.subheader("Settings")
    c1, c2 = st.columns(2)
    with c1:
        theta_min = st.number_input("theta_incl min (deg)", value=0.0, key="s2_tmin")
        theta_max = st.number_input("theta_incl max (deg)", value=40.0, key="s2_tmax")
        theta_step = st.number_input("theta_incl step (deg)", value=1.0, key="s2_tstep")
    with c2:
        phi_min = st.number_input("phi_incl min (deg)", value=0.0, key="s2_pmin")
        phi_max = st.number_input("phi_incl max (deg)", value=360.0, key="s2_pmax")
        phi_step = st.number_input("phi_incl step (deg)", value=5.0, key="s2_pstep")
    st.caption("Starting structure: R-form from Step 1 (alpha = 25°)")

    if st.button("Run DFT Jobs", key="s2_run"):
        _not_implemented("Step 2 (inclination map)", "§Step 2")

    st.divider()
    st.caption(
        "2D heatmap of E_intra(6): x = theta_incl·cos(phi_incl), "
        "y = theta_incl·sin(phi_incl); markers: R-form / G-form / N-form"
    )
    c5, c6, c7 = st.columns(3)
    c5.text_input("theta_incl (N-form)", key="s2_theta_n", disabled=True)
    c6.text_input("phi_incl (N-form)", key="s2_phi_n", disabled=True)
    c7.text_input("theta_incl (G-form)", key="s2_theta_g", disabled=True)


# ══════════════════════════════════════════════════════════
#  Tab 5: Step 3 – Interlayer Stacking
# ══════════════════════════════════════════════════════════
with tab_step3:
    st.subheader("Settings")
    form = st.radio("Form", ["R-form", "G-form", "N-form"], horizontal=True, key="s3_form")
    c1, c2 = st.columns(2)
    x_step = c1.number_input("x step (Å, range -a/2..a/2)", value=0.1, key="s3_xstep")
    y_step = c2.number_input("y step (Å, range -b/2..b/2)", value=0.1, key="s3_ystep")
    st.caption("z: auto-computed from the VdW contact distance")

    c3, c4 = st.columns(2)
    if c3.button("Run vdW Volume Scan", key="s3_vdw_run"):
        _not_implemented("Step 3 (vdW volume map)", "§Step 3")
    if c4.button("Run DFT E_inter(7)", key="s3_dft_run"):
        _not_implemented("Step 3 (DFT-D interlayer energy)", "§Step 3")

    st.divider()
    col_l, col_r = st.columns(2)
    with col_l:
        st.caption("V(x,y) map (vdW volume) — markers: N1 (Type II) / N2 (Type IV)")
    with col_r:
        st.caption("E_inter(7) map (DFT-D)")
    st.selectbox("Select stacking for Step 4", ["N1", "N2"], key="s3_stacking_choice")


# ══════════════════════════════════════════════════════════
#  Tab 6: Step 4 – Refinement
# ══════════════════════════════════════════════════════════
with tab_step4:
    sub_twist, sub_incl = st.tabs([
        "6a. Twist (G-form → Type III)",
        "6b. Non-uniform Inclination (N-form → Type IV)",
    ])

    with sub_twist:
        c1, c2, c3 = st.columns(3)
        c1.number_input("theta_twist min (deg)", value=0.0, key="s4a_min")
        c2.number_input("theta_twist max (deg)", value=20.0, key="s4a_max")
        c3.number_input("theta_twist step (deg)", value=1.0, key="s4a_step")
        st.caption("Optimize: a, b, ΔzT, x, y, z")
        if st.button("Run DFT Jobs", key="s4a_run"):
            _not_implemented("Step 4a (twist refinement)", "§Step 4a")
        st.caption("Plot: theta_twist vs E_int(near); 3D structure with optimal theta_twist marked")

    with sub_incl:
        c1, c2, c3 = st.columns(3)
        c1.number_input("theta'_incl min (deg)", value=0.0, key="s4b_min")
        c2.number_input("theta'_incl max (deg)", value=5.0, key="s4b_max")
        c3.number_input("theta'_incl step (deg)", value=0.5, key="s4b_step")
        st.caption("phi'_incl: auto from Step 4a's E_intra(6) minimum direction")
        if st.button("Run DFT Jobs", key="s4b_run"):
            _not_implemented("Step 4b (non-uniform inclination)", "§Step 4b")
        st.caption(
            "Plot: theta'_incl vs E_int(near); table: Type II vs Type IV "
            "comparison for tetracene/pentacene/hexacene"
        )


# ══════════════════════════════════════════════════════════
#  Tab 7: Transfer Integrals
# ══════════════════════════════════════════════════════════
with tab_step5:
    st.subheader("Select structures")
    c1, c2, c3, c4 = st.columns(4)
    c1.checkbox("Type I (R-form)", value=True, key="s5_type1")
    c2.checkbox("Type II (N1)", value=True, key="s5_type2")
    c3.checkbox("Type III (G+twist)", value=True, key="s5_type3")
    c4.checkbox("Type IV (N2)", value=True, key="s5_type4")

    st.code("B3LYP/6-31g*  (HOMO -> transfer integral J)", language="text")
    c5, c6 = st.columns(2)
    if c5.button("Run Gaussian (B3LYP/6-31g*)", key="s5_run"):
        _not_implemented("Step 5 (transfer integrals)", "§Step 5")
    if c6.button("Check Status", key="s5_status"):
        _not_implemented("Step 5 job status check", "§Step 5")

    st.divider()
    st.caption("Bar chart: J (meV) for T-type and SP-type contacts per polymorph")
    st.caption("Plot: J vs alpha (R-form, varying alpha)")
    st.caption("Plot: J vs theta_incl (N-form, varying theta_incl)")
    st.caption(
        "Uses Prof. Matsui's HOMO-HOMO overlap integral code "
        "(not yet integrated — pending legacy/ono_scripts/)."
    )
