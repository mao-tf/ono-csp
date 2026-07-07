# csp — Crystal Structure Prediction for Polyacenes

A Python package + Streamlit GUI that reproduces the stepwise crystal
structure prediction pipeline of Ono et al., *"Origin of Layered
Herringbone Packing and Polymorphism in Polyacenes: A Quantum Chemical
Optimization Approach"* (JACS, submitted).

Instead of a brute-force search over every structural parameter at once,
the method separates variables by which physical interaction dominates
them (intralayer packing first, then interlayer stacking) and optimizes
them stepwise, starting from a cheap rigid-sphere van der Waals model and
refining with dispersion-corrected DFT (B3LYP-D3/6-311G\*\*). This
reproduces the observed crystal structures of naphthalene, anthracene,
tetracene, pentacene, and hexacene from first principles, without
fitting to experimental unit cells.

## What the pipeline does

| Step | What it finds | Method |
|------|----------------|--------|
| 1 | Intralayer packing (herringbone half-angle α, lattice constants a, b) | rigid-sphere vdW rough scan → DFT-D refinement |
| 2 | Uniform long-axis inclination (θ_incl, φ_incl) — the R/G/N-form branch point | DFT-D energy map, reconstructed from 1D scans |
| 3 | Interlayer stacking offset (x, y, z between layers) | rigid-sphere vdW map → DFT-D refinement |
| 4 | Long-axis twist / non-uniform inclination refinements (Type III/IV polymorphs) | DFT-D |
| 5 | Intermolecular transfer integrals (charge-transport properties) | Gaussian16 MOs, via [tcal](https://github.com/matsui-lab-yamagata/tcal) (Matsui Lab.) |

The rigid-sphere vdW scans (Steps 1 and 3) run directly in the GUI in
under a second. The DFT-D refinement steps need Gaussian16 and are meant
to run on an HPC cluster; the GUI shows the exact CLI commands for each
step and displays whatever result CSV you drag & drop back in.

## Install

```bash
conda create -n csp python=3.10
conda activate csp
pip install -r requirements.txt
pip install -e .
```

Requires Python ≥ 3.10. For just the vdW-scan/structure-building parts of
the package (no GUI), the core dependencies (`numpy`, `pandas`, `scipy`,
`pyyaml`) are enough — see `pyproject.toml`'s `[project.optional-dependencies]`
for the `viz` extras (`streamlit`, `plotly`, `py3Dmol`, `matplotlib`)
needed for the GUI itself.

## Quickstart

```bash
streamlit run app.py
```

This opens a 5-tab GUI:

1. **Molecule Setup** — pick a preset polyacene (naphthalene, anthracene,
   tetracene, pentacene, hexacene) or upload your own monomer XYZ, adjust
   van der Waals radii, and preview it in 3D. Download the monomer CSV
   needed by the CLI scripts below.
2. **Step 1 – Intralayer** — runs the vdW rough scan live in the GUI
   (a few seconds), showing candidate (a, b, α) combinations and a
   9-molecule cluster preview. For the DFT-D refinement, drop the
   resulting `step1.csv` back in to see the optimized α (a
   branch-classified energy-vs-α plot, reconstructed to cover the full
   0–90° range even from a partial scan).
3. **Step 2 – Long-Axis Shift** — `para` and `twist` sub-tabs. `para`
   shows the T-shaped/slipped-parallel dimer energy scans and
   reconstructs the 2D (θ_incl, φ_incl) landscape (R/G/N-form map) from
   them, with a clickable 3D preview of any point. `twist` reconstructs
   the analogous 2D map for the long-axis twist refinement.
4. **Step 3 – Interlayer Stacking** — `para` and `twist` sub-tabs; the
   vdW interlayer scan runs live in the GUI, showing the unit-cell-volume
   map and a two-layer 3D preview.
5. **Transfer Integrals** — shows the CLI command for computing
   HOMO-HOMO transfer integrals at the T-shaped and slipped-parallel
   contacts, and plots the resulting values.

Every result-CSV section works the same way: drop your own CSV to plot
it, or leave it empty to see a bundled pentacene sample
(`example/pentacene/`) so you can see what a tab looks like with real
data before running anything yourself.

## Running the DFT steps

The DFT-D refinement steps (everything past the in-GUI vdW scans) run
from the command line using the scripts in `legacy/ono_scripts/`, on
whatever HPC cluster you have Gaussian16 access to. Each GUI tab shows
the exact command plus what input file to prepare and what the output
looks like. Two things to set up once:

- `CSP_MONOMER_DIR`: point this at the folder holding your monomer CSV
  (`export CSP_MONOMER_DIR=/path/to/data/monomer`) — download it from
  Tab 1.
- Job scheduler: the scripts are uploaded as-is, targeting the Fujitsu
  PJM scheduler (`pjsub`) they were originally written for. Adapt the
  batch-script generation and submission call to whatever your own
  cluster uses (`qsub` for SGE/PBS, `sbatch` for Slurm, ...);
  `src/csp/dft/job_cluster.py` has a worked SGE example if useful as a
  reference.

### Getting a single-point energy instead of a hill-climb (Step 2 twist)

`step2_twist.py`'s `get_opt_params_dict` normally hill-climbs (a, b) in
±0.1 Å steps at each fixed (θ, Rt, A2). To instead get just the energy
*at* the (a, b) you supply — e.g. to build the paper's Fig. 7(c)-style
landscape from vdW-derived starting structures without also
re-optimizing (a, b) at every grid point — change the `±0.1` step sizes
in that function to `0`. With a zero step, the "neighbor grid" collapses
to the single input point, so the hill-climb reports convergence
immediately and returns the input (a, b) unchanged; only its energy gets
computed. This isn't documented anywhere else, so it's easy to miss.

## Layout

```
src/csp/
├── structure/   molecule loading, glide-symmetric layer/cluster construction
├── vdw/         rigid-sphere vdW contact scans (Step 1/3 pre-scans)
├── dft/         Gaussian16 input generation, log parsing, HPC job submission
└── plot/        py3Dmol/Plotly viewers and figure-reconstruction logic

data/molecules/  B3LYP-D3/6-311G** optimized polyacene monomer XYZ files
example/pentacene/  bundled sample result CSVs (shown until you drop your own)
legacy/          Ono's original CUI scripts, kept as uploaded
```

## License

MIT (see `LICENSE`). The bundled `legacy/ono_scripts/tcal_csv/tcal_1.py` is a
modified copy of [tcal](https://github.com/matsui-lab-yamagata/tcal) by
Matsui Lab. at Yamagata University, also MIT-licensed (notice included in
`LICENSE`).
