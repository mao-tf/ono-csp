# csp — Crystal Structure Prediction

Reproduces the crystal structure search pipeline of Ono et al., *"Origin of
Layered Herringbone Packing and Polymorphism in Polyacenes: A Quantum
Chemical Optimization Approach"* (JACS, submitted).

See `spec.md` for the full method-to-code mapping (paper section ↔ package
step ↔ UI tab), and `CLAUDE.md` for background on what is reused from the
companion `auto_opt` package.

## Status

Early scaffolding. Tab 1 (Molecule Setup) of the UI is functional; the
Step 1–5 calculation engines (Tabs 2–7) are UI skeletons pending
implementation — see each tab's caption, and `legacy/` for Ono's existing
scripts once integrated.

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

MIT
