import subprocess
import numpy as np
import shutil
from utils import Rod, R2atom
import time
import os
import pandas as pd
import argparse
import subprocess

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

MONOMER_LIST = ["naphthalene","anthracene","tetracene","pentacene","hexacene"]
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(os.path.join(os.path.expanduser(os.environ.get('CSP_MONOMER_DIR', '~/path/to/monomer/')), f'{monomer_name}.csv'))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])
    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    return np.concatenate([xyz_array,R_array],axis=1)
    
def pre_process(args):
    work_dir = args.auto_dir
    monomer_name = args.monomer_name

    df=pd.read_csv(f'{work_dir}/init_params.csv')##parameters of molecular arrangements to be calculated
    a_list=df['a'].values.tolist();b_list=df['b'].values.tolist();alpha_list=df['alpha'].values.tolist()
    A2_list=df['A2'].values.tolist();z_list=df['z'].values.tolist()
    for a,b,alpha,A2,z in zip(a_list,b_list,alpha_list,A2_list,z_list):
        
        os.makedirs(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}",exist_ok=True)## calculation directory
        shutil.copy(os.path.join(_SCRIPT_DIR,"job.sh"),rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}") ## batch job script
        shutil.copy(os.path.join(_SCRIPT_DIR,"tcal_1.py"),rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}") ## tcal program for calculating the transfer integral  
        
    for a,b,alpha,A2,z in zip(a_list,b_list,alpha_list,A2_list,z_list):  ## Gaussian16 input files without molecuar coordinates  
        lines1 = [     
            '%mem=30GB\n',
            '%nproc=40\n',
            '%chk=test.chk\n',
            '#b3lyp/6-31g*\n',
            '# symmetry=none\n',
            '\n',
            'alpha={}_t\n'.format(alpha),
            '\n',
            '0 1\n'
        ]
        lines2 = [
            '\n',
            '--Link1--\n',
            '%mem=30GB\n',
            '%nproc=40\n',
            '%chk=test.chk\n',
            '#b3lyp/6-31g*\n',
            '# symmetry=none\n',
            '# geom=allcheck\n',
            '# guess=read\n',
            '# pop=full\n',
            '# iop(3/33=4,5/33=3)\n\n\n'
        ]

        monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,alpha)## central molecule
        monomer_array_t = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-alpha)
        if b>a:
            monomer_array_p = get_monomer_xyzR(monomer_name,a,0,0,A2,alpha)
        else:
            monomer_array_p = get_monomer_xyzR(monomer_name,0,b,2*z,A2,alpha)###molecular shift under glide symmetry
        
        line_i=[];line_t=[];line_p=[]
        for X,Y,Z,R in monomer_array_i:
            atom = R2atom(R)
            line = '{} {} {} {}\n'.format(atom,X,Y,Z)     
            line_i.append(line)
        for X,Y,Z,R in monomer_array_t:
            atom = R2atom(R)
            line = '{} {} {} {}\n'.format(atom,X,Y,Z)     
            line_t.append(line)
        for X,Y,Z,R in monomer_array_p:
            atom = R2atom(R)
            line = '{} {} {} {}\n'.format(atom,X,Y,Z)     
            line_p.append(line)

        ### 6 files consisting of dimer monomer1 monomer2 of t-shaped and slipped parallel
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t.gjf","w") as f:
            f.writelines(lines1)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t.gjf","a") as f:
            f.writelines(line_i)
            f.writelines(line_t)
            f.writelines(lines2)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t_m1.gjf","w") as f:
            f.writelines(lines1)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t_m1.gjf","a") as f:
            f.writelines(line_i)
            f.writelines(lines2)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t_m2.gjf","w") as f:
            f.writelines(lines1)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t_m2.gjf","a") as f:
            f.writelines(line_t)
            f.writelines(lines2)

        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p.gjf","w") as f:
            f.writelines(lines1)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p.gjf","a") as f:
            f.writelines(line_i)
            f.writelines(line_p)
            f.writelines(lines2)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p_m1.gjf","w") as f:
            f.writelines(lines1)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p_m1.gjf","a") as f:
            f.writelines(line_i)
            f.writelines(lines2)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p_m2.gjf","w") as f:
            f.writelines(lines1)
        with open(rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p_m2.gjf","a") as f:
            f.writelines(line_p)
            f.writelines(lines2)
        
def qsub_process(args):## submission of Gaussian16 calculations
    work_dir = args.auto_dir
    
    df=pd.read_csv(f'{work_dir}/init_params.csv')
    a_list=df['a'].values.tolist();b_list=df['b'].values.tolist();alpha_list=df['alpha'].values.tolist()
    A2_list=df['A2'].values.tolist();z_list=df['z'].values.tolist()
    for a,b,alpha,A2,z in zip(a_list,b_list,alpha_list,A2_list,z_list):
        path=rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}"
        os.chdir(path)
        subprocess.run(['pjsub job.sh'],shell=True)
        time.sleep(10)##we do not prepare controlling system so we submit the job every 10 seconds

def tcal_process(args): ## calculate transfer integral using tcal_1.py published by Prof. Matsui
    work_dir = args.auto_dir
    monomer_name = args.monomer_name
    df=pd.read_csv(f'{work_dir}/init_params.csv')
    a_list=df['a'].values.tolist();b_list=df['b'].values.tolist();alpha_list=df['alpha'].values.tolist()
    A2_list=df['A2'].values.tolist();z_list=df['z'].values.tolist()

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,0,0)
    n1=len(monomer_array_i)

    for a,b,alpha,A2,z in zip(a_list,b_list,alpha_list,A2_list,z_list):
        path=rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}"
        check_path=rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test.fch"
        while not os.path.exists(check_path):
            time.sleep(1)
        os.chdir(path)
        subprocess.run([f'python tcal_1.py -r --n_mono1 {n1} test_t.gjf > test_t.txt &'],shell=True)
        subprocess.run([f'python tcal_1.py -r --n_mono1 {n1} test_p.gjf > test_p.txt &'],shell=True)

def result_process(args):
    work_dir = args.auto_dir
    
    df=pd.read_csv(f'{work_dir}/init_params.csv')##pathは適宜決定
    a_list=df['a'].values.tolist();b_list=df['b'].values.tolist();alpha_list=df['alpha'].values.tolist()
    A2_list=df['A2'].values.tolist();z_list=df['z'].values.tolist()
    with open(f"{work_dir}/result.txt",'w')as f1:
        for a,b,alpha,A2,z in zip(a_list,b_list,alpha_list,A2_list,z_list):
            path_t=rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_t.txt"
            path_p=rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}/test_p.txt"
            with open(path_t) as f:
                for line in f:
                    s=line.split()
                    if (len(s) == 3) and (s[0] =='HOMO'):
                        f1.write('{} '.format(s[1]))
            with open(path_p) as f:
                for line in f:
                    s=line.split()
                    if (len(s) == 3) and (s[0] =='HOMO'):
                        f1.write('{} '.format(s[1]))
            f1.write('\n')

def rm_process(args):## log files of this workflow are large, so you can remove them if you want.   
    work_dir = args.auto_dir
    df=pd.read_csv(f'{work_dir}/init_params.csv')##
    a_list=df['a'].values.tolist();b_list=df['b'].values.tolist();alpha_list=df['alpha'].values.tolist()
    A2_list=df['A2'].values.tolist();z_list=df['z'].values.tolist()
    for a,b,alpha,A2,z in zip(a_list,b_list,alpha_list,A2_list,z_list):
        path=rf"{work_dir}/a={a}_b={b}_alpha={alpha}_A2={A2}_z={z}"
        os.chdir(path)
        subprocess.run(['rm *log'],shell=True)## you can add 'rm *chk' and 'rm *fch'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--init',action='store_true')
    parser.add_argument('--test',action='store_true')
    parser.add_argument('--qsub',action='store_true')
    parser.add_argument('--tcal',action='store_true')
    parser.add_argument('--result',action='store_true')
    parser.add_argument('--rm',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    args = parser.parse_args()

    ## each steps can be separately exected with proper usage of options.
    if args.init:
        print("----initial process----")
        pre_process(args)

    if args.test:
        print("----test has finished----")
    else:
        if args.qsub:
            print("----qsub process----")
            qsub_process(args)
        if args.tcal:
            print("----tcal process----")
            tcal_process(args)
        if args.result:
            print("----result process----")
            result_process(args)
        if args.rm:
            print("----remove process----")
            rm_process(args)
        print("----all processes have finished----")
    