# csp — Crystal Structure Prediction

Reproduces the crystal structure search pipeline of Ono et al., *"Origin of
Layered Herringbone Packing and Polymorphism in Polyacenes: A Quantum
Chemical Optimization Approach"* (JACS, submitted).

See `spec.md` for the full method-to-code mapping (paper section ↔ package
step ↔ UI tab), and `CLAUDE.md` for background on what is reused from the
companion `auto_opt` package.

## Status

Working prototype:

- **Tab 1 (Molecule Setup)** — preset/custom molecule loading, editable vdW
  radii, 3D viewer, and export of the monomer CSV consumed by the CLI scripts.
- **Tab 2 (Intralayer vdW Scan)** — runs in the GUI: rigid-sphere contact
  scan over the herringbone half-angle, alpha-vs-S plot, 9-molecule cluster
  preview, and download of the initial candidates for the DFT step.
- **Tabs 3–7 (DFT steps / transfer integrals)** — these run from the command
  line with the scripts in `legacy/ono_scripts/` (Gaussian16 on an HPC); each
  tab shows the commands and visualizes the result CSVs (drag & drop, with
  bundled samples as fallback).

See `spec.md` (section "大野コード対応表") for the mapping between the legacy
scripts, the `src/csp` modules, and the GUI tabs.

## Install

```bash
conda create -n csp-py310 python=3.10
conda activate csp-py310
pip install -r requirements.txt
pip install -e .
```

## Run

```bash
streamlit run app.py
```

## Layout

```
src/csp/
├── structure/   molecule loading, glide-symmetric placement
├── vdw/         rigid-sphere VdW contact distance + scans
├── dft/         Gaussian16 input generation, log parsing, SGE job submission
├── transfer/    transfer-integral (HOMO-HOMO overlap) calculation
└── plot/        py3Dmol / Plotly viewers for the Streamlit UI
data/molecules/  B3LYP-D3/6-311G** optimized polyacene monomer XYZ files
legacy/          staging area for Ono's existing CUI scripts
```

## License

MIT (see `LICENSE`). The bundled `legacy/ono_scripts/tcal_csv/tcal_1.py` is a
modified copy of [tcal](https://github.com/matsui-lab-yamagata/tcal) by
Matsui Lab. at Yamagata University, also MIT-licensed (notice included in
`LICENSE`).
