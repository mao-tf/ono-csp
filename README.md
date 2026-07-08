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

## Getting the code (no git experience needed)

You don't need to know git to try this out — just download a snapshot of
the code as a zip file:

1. Go to <https://github.com/mao-tf/ono-csp> in your browser.
2. Click the green **Code** button, then **Download ZIP**.
3. Unzip it wherever you like (double-click it on Mac/Windows). You'll
   get a folder named `ono-csp-main`.
4. Open a terminal (**Terminal** app on Mac, **PowerShell** or
   **Command Prompt** on Windows) and move into that folder, e.g.:
   ```bash
   cd ~/Downloads/ono-csp-main
   ```
5. Continue with **Install** below.

If you later want to fetch updates without re-downloading the whole zip
each time, it's worth installing git instead (one-time setup):

- Mac: open Terminal and run `git --version` — if it's not installed,
  macOS will prompt you to install the Xcode Command Line Tools
  automatically.
- Windows: install [Git for Windows](https://git-scm.com/download/win),
  then use **Git Bash** as your terminal.
- Then, once, download the code with:
  ```bash
  git clone https://github.com/mao-tf/ono-csp.git
  cd ono-csp
  ```
- Any time later, get the newest updates with:
  ```bash
  git pull
  ```
  (run this from inside the `ono-csp` folder)

## Install

**Run every command below from inside the project folder** (the one
from **Getting the code** above — `ono-csp-main` or `ono-csp`, whichever
you ended up with). If you're not sure you're still there, check with:

```bash
ls pyproject.toml    # should print "pyproject.toml", not an error
```

If that errors, `cd` into the folder first (e.g. `cd ~/Downloads/ono-csp-main`).

Requires [Python](https://www.python.org/) 3.10+. Either way below gives
you an isolated environment so this project's dependencies don't clash
with anything else on your machine.

**Have neither Python nor conda?** Just install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) — its
installer sets up a matching Python *and* conda together, so it's the one
thing you need. Download the installer for your OS, run it, then open a
**new** terminal window (so the install takes effect) before continuing
with **With conda** below.

**With conda** (install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) first if you
don't have it). **On Windows**, use the **Anaconda Prompt** app (search
for it in the Start menu) instead of the regular Command Prompt or
PowerShell — those don't recognize the `conda` command:

```bash
conda create -n csp python=3.10
conda activate csp
pip install -e ".[viz]"
```

**Without conda**, using Python's built-in `venv` instead (install
[Python](https://www.python.org/downloads/) 3.10+ first if you don't have
it). First check what `python3 --version` gives you — many systems have
an older Python 3 (e.g. 3.8/3.9) as the default `python3`, in which case
use whatever command your 3.10+ install actually goes by instead
(`python3.10`, `python3.11`, ...):

```bash
python3 --version                # confirm it says 3.10 or higher
python3 -m venv csp-env
source csp-env/bin/activate      # Windows: csp-env\Scripts\activate
pip install -e ".[viz]"
```

If that fails with `requires a different Python: ... not in '>=3.10'`,
that's this exact mismatch — the `venv` was built from a too-old
`python3`. Delete the `csp-env` folder, rerun the first command with the
right interpreter (e.g. `python3.11 -m venv csp-env`), and try again.

`pip install -e ".[viz]"` pulls in everything the GUI needs (including
Streamlit itself) in one step. (`requirements.txt` also exists and
installs the same set of packages — either is fine, but don't mix the
two up: running only `pip install -e .` *without* `[viz]` installs just
the core numpy/pandas/scipy library, not Streamlit/Plotly/py3Dmol, and
the GUI will fail to start.)

For just the vdW-scan/structure-building parts of the package (no GUI),
the core dependencies (`numpy`, `pandas`, `scipy`, `pyyaml`) are enough —
see `pyproject.toml`'s `[project.optional-dependencies]` for the `viz`
extras (`streamlit`, `plotly`, `py3Dmol`, `matplotlib`) needed for the GUI
itself.

## Quickstart

Still in the same terminal as the Install step above (same folder, same
activated environment)? Just run:

```bash
streamlit run app.py
```

**Next time** (after closing the terminal and coming back later): you
don't need to reinstall anything, but you do need to re-activate the
environment first, since that doesn't carry over between terminal
sessions —

```bash
conda activate csp               # or: source csp-env/bin/activate
cd /path/to/ono-csp              # wherever you put the code
streamlit run app.py
```

Forgetting the activate step is the usual reason `streamlit run app.py`
suddenly says `streamlit: command not found`.

This opens a 5-tab GUI. A molecule picker (preset polyacene or your own
uploaded XYZ) sits above the tabs and applies everywhere — no need to
switch tabs just to change molecules.

1. **Molecule Setup** — vdW radius table, a 3D sanity-check view of the
   currently picked molecule, and a monomer CSV download needed by the
   CLI scripts below.
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

### What a working directory (`--auto-dir`) looks like

Each command creates/uses a working directory (the path you pass to
`--auto-dir`) with the same shape, adapted from Ono's own
`legacy/ono_scripts/stepwise_optimization/readme.txt`:

```
your-working-dir/
├── stepX_init_params.csv   # starting points -- from a GUI download, or hand-written
├── stepX.csv               # results, created/updated as jobs finish
├── gaussian/               # created automatically
│   ├── *.inp               # Gaussian16 input files
│   ├── *.log               # Gaussian16 output files
│   └── *.r1                # batch job scripts
└── gaussview/              # created automatically
    └── *.xyz               # structures, for visualization
```

For the transfer-integral workflow (`legacy/ono_scripts/tcal_csv/`), the
same idea but one subdirectory per arrangement (from
`legacy/ono_scripts/tcal_csv/readme.txt`):

```
your-working-dir/
├── init_params.csv                # arrangements to compute (a, b, theta, A2, z)
├── result.txt                     # summarized transfer integrals, after --result
└── a=.._b=.._theta=.._A2=.._z=../ # one folder per arrangement
    ├── job.sh, tcal_1.py          # copied in automatically
    ├── test_t/test_p*.gjf         # Gaussian16 inputs (T-shaped / slipped-parallel)
    ├── test_t/test_p*.log         # Gaussian16 outputs
    └── test_t/test_p.txt          # per-arrangement transfer integral
```

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
