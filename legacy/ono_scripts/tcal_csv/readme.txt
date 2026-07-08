This workflow automates the calculation of transfer integrals for arbitrary molecular arrangements using Gaussian 16 and the existing tcal program.
First, molecular arrangements are extracted from the output files of the stepwise optimization. Next, molecular orbital information is obtained from Gaussian 16 calculations, and transfer integrals are calculated from the Gaussian 16 output files.
This workflow utilizes the tcal program, which enables analysis of atomic-pair contributions to the transfer integrals. The primary purpose of this workflow is to automate transfer integral calculations over a wide range of molecular arrangements and calculation parameters. However, it can also be used as a general workflow for transfer integral analysis.

More information about the tcal program is available at:
https://github.com/matsui-lab-yamagata/tcal

Please note that tcal_1.py included in this workflow is a slightly modified version of the original tcal.py.

Example of command-line
nohup python /path/to/tcal_csv/tcal_csv.py --monomer-name pentacene --auto-dir pentacene_HB --init --qsub --tcal --reult & 

Directory structures
tcal_csv/
├── tcal_1.py
├── tcal_csv.py
├── utils.py
├── job.sh ##you should use your batch job script and Gaussian16 settings.
└── Working directory/
    ├── result_params.csv ##parameter sets of molecular arrangements
    ├── result.txt (summarized calculation results)
    (generated after execution)
    └── dir_with_each_parameters /
        ├── job.sh (copied)
        ├── tcal_1.py (copied)
        ├── test_t or _p (_m1 or _m2).gjf (Gaussian16 input files for dimer and monomer of T-shaped and slipped parallel contact)
        ├── test_t or _p (_m1 or _m2).out (Gaussian16 output files)
        └── test_t or _p.txt (result file of transfer integral calculations)