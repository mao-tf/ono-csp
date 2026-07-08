Overview
The Python scripts in this workflow generate molecular arrangements (.xyz files), Gaussian16 input files, and batch job scripts. These files are created by the make_*.py scripts.
The molecular optimization based on a hill-climbing algorithm is performed by the stepX_*.py scripts. These scripts determine the optimization parameters, execute Gaussian16 calculations, and summarize the calculated DFT energies.

Optimization Workflow
Step 1: Optimize the arrangement of eight adjacent molecules within a single layer without molecular displacement along the molecular long axis.
Step 2: Optimize the molecular displacement and torsion.
    step2_twist.py is used for Type-3 packing with glide symmetry.
    step2_para.py is used for Type-2 and Type-4 packing.
Step 3: Optimize the interlayer molecular arrangement (the c-axis vector).
Separate optimization workflows are provided for both the twist and parallel (para) packing structures.

Utility Functions
utils.py provides several utility functions.
For example:
    get_E() extracts intermolecular interaction energies from Gaussian output (.log) files.
    Rod rotates molecules for generating molecular arrangements.

Notes on step2_para.py
The contents of the make_*.py and stepX_*.py scripts are largely identical, except for step2_para.py.
In step2_para.py, the intermolecular potential-energy map is calculated using a simplified procedure in which the lattice parameters are fixed to the optimized structure obtained in Step 1.

Directory Setup
Prepare the following directories before running the workflow:
    source_directory: contains the Python scripts.
    Working_directory: stores optimization results.
During execution, the calculated energies are saved in stepX.csv in the working directory.
The program also creates the following subdirectories automatically:
    gaussian/ (Gaussian16 input, output, and batch files)
    gaussview/ (.xyz files for visualization)

Command-line Options
--auto-dir Specifies the working directory where the optimization is performed. The CSV files are written to this directory, and the gaussian/ and gaussview/ subdirectories are created automatically.
--monomer-name Specifies the name of the monomer to be optimized.
--num-nodes Specifies the maximum number of Gaussian jobs that can be executed simultaneously.
--num-init Specifies the number of initial structures (parameter sets) optimized simultaneously.

Example of Command-line
nohup python /path/to/stepX/source_directory/step_X.py --monomer-name pentacene --auto-dir /path/to/step_X/Working_directory --num-init N1 --num-nodes N2 &
Using nohup together with & allows the optimization to continue running in the background after logging out.

Directory Structures
stepX/
├── source_directory/
│   ├── make_step_X.py
│   ├── step_X.py
│   └── utils.py
└── Working_directory/
    ├── stepX_init_params.csv # Can be created from scratch in Step 1
    ├── stepX.csv # Generated automatically (an existing file may also be reused)
    ├── gaussian/ # Generated during execution
    │   ├── input files (.inp)
    │   ├── output files (.log)
    │   └── batch job scripts (.r1)
    └── gaussview/ # Generated during execution
        └── visuallization files (.xyz)