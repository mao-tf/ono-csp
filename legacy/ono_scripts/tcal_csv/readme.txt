With this workflow, we can calculate transfer integrals of arbitrary molecular arrangements with Gaussian16 calculations and existing programm.
We first identify the molecular arrangements from the result files of stepwise optimization.
Then we obtain molecular orbital data from Gaussian16 calculations and calculate transfer integrals utilizing Gaussian16 output files.    

We utilize tcal.py program in this workflow and this program enables us to analyze the atomic pair contiribution.
Our main focus is the automation of the transfer integral calculation with a broad parameter space, but you can use this workflow as transfer integral analysis.

The detail information about tcal program is in https://github.com/matsui-lab-yamagata/tcal.
Please note that tcal_1.py has small change from older version of tcal.py.


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