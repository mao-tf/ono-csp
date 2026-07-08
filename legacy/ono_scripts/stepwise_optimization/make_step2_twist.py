import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ["BTBT","naphthalene","anthracene","tetracene","pentacene","hexacene"]
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):## rotation and translation operations
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(os.path.join(os.path.expanduser(os.environ.get('CSP_MONOMER_DIR', '~/path/to/monomer/')), '{}.csv'.format(monomer_name)))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)## A2 is a torsion angle in step2
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)## A3 is a half of dihedral angle
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    if monomer_name in MONOMER_LIST:
        return np.concatenate([xyz_array,R_array],axis=1)
    
    else:
        raise RuntimeError('invalid monomer_name={}'.format(monomer_name))
        
def get_xyzR_lines(xyzR_array,file_description):## Gaussian input files
    lines = [     
        '%mem=24GB\n',
        '%nproc=48\n',
        '#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2\n',
        '\n',
        file_description+'\n',
        '\n',
        '0 1 0 1 0 1\n'
    ]
    mol_len = len(xyzR_array)//2
    atom_index = 0
    mol_index = 0
    for x,y,z,R in xyzR_array:
        atom = R2atom(R)
        mol_index = atom_index//mol_len + 1
        line = '{}(Fragment={}) {} {} {}\n'.format(atom,mol_index,x,y,z)     
        lines.append(line)
        atom_index += 1
    return lines

def get_one_exe(file_name):## batch job scripts
    file_basename = os.path.splitext(file_name)[0]
    cc_list=[
        ## you can use your environment for Gaussian16 calculations
        '\n',
        'g16 < {}.inp > {}.log \n'.format(file_basename,file_basename),
        '\n'
        ## you can use your environment for Gaussian16 calculations
            ]

    return cc_list

def make_xyzfile(monomer_name,params_dict):#.xyz files for visualization
    a_ = params_dict['a']; b_ = params_dict['b']
    A2 = params_dict.get('A2',0.0); A3 = params_dict['alpha']
    
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3) # central molecule
    monomer_array_p1 = get_monomer_xyzR(monomer_name,0,b_,0,A2,A3) #p1 - p4 are slipped parallel contact molecules
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,-b_,0,A2,A3)
    monomer_array_p3 = get_monomer_xyzR(monomer_name,a_,0,0,A2,A3)
    monomer_array_p4 = get_monomer_xyzR(monomer_name,-a_,0,0,A2,A3)
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2,b_/2,0,A2,-A3) # t1 - t4 are t-shaped contact molecules
    monomer_array_t2 = get_monomer_xyzR(monomer_name,a_/2,-b_/2,0,A2,-A3)
    monomer_array_t3 = get_monomer_xyzR(monomer_name,-a_/2,-b_/2,0,A2,-A3)
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a_/2,b_/2,0,A2,-A3)
    
    monomers_array = np.concatenate([monomer_array_i,monomer_array_p1,monomer_array_p3,monomer_array_p2,monomer_array_p4,monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4],axis=0)
    xyz_list=[f'{len(monomers_array)} \n','polyacene9 \n']##xyz_file with 9 molecules
    
    for x,y,z,R in monomers_array:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a','b','cx','cy','cz','alpha']:
            val = np.round(val,2)
        elif key in ['A2']:#,'theta']:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):#.xyz files for visualization
    a_ = params_dict['a']; b_ = params_dict['b']
    Rt = params_dict['Rt']; A2 = params_dict['A2']; A3 = params_dict['alpha']
    
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b_,2*Rt,A2,A3)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a_,0,0,A2,A3)
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2,b_/2,Rt,A2,-A3)
    
    dimer_array_t1 = np.concatenate([monomer_array_i,monomer_array_t1])
    dimer_array_p1 = np.concatenate([monomer_array_i,monomer_array_p1])
    dimer_array_p2 = np.concatenate([monomer_array_i,monomer_array_p2])
    
    file_description = '{}_Rt={}_A2={}_A3={}'.format(monomer_name,round(Rt,1),int(A2),round(A3,2))
    line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1')
    line_list_dimer_p2 = get_xyzR_lines(dimer_array_p2,file_description+'_p2')
    line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1')
    
    gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer_t1 + ['\n\n--Link1--\n'] + line_list_dimer_p1 + ['\n\n\n']#+ ['\n\n--Link1--\n'] + line_list_dimer_p2 + ['\n\n\n']
    
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    inp_dir = os.path.join(auto_dir,'gaussian')
    gij_xyz_path = os.path.join(inp_dir,file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,paras_dict):##basic file name for input, output and job script files
    file_name = ''
    file_name += monomer_name
    for key,val in paras_dict.items():
        if key in ['a','b','Rt','alpha']:
            val = np.round(val,2)
        elif key in ['A2']:#,'theta']:
            val = np.round(val,2)
        file_name += '_{}={}'.format(key,val)
    return file_name + '.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict,isTest=True):##submission of batch job 
    inp_dir = os.path.join(auto_dir,'gaussian')
    xyz_dir = os.path.join(auto_dir,'gaussview')
    
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    cc_list = get_one_exe(file_name)
    sh_filename = os.path.splitext(file_name)[0]+'.r1'
    sh_path = os.path.join(inp_dir,sh_filename)
    with open(sh_path,'w') as f:
        f.writelines(cc_list)
    if not(isTest):
        subprocess.run(['pjsub',sh_path])
    log_file_name = os.path.splitext(file_name)[0]+'.log'
    return log_file_name
    
############################################################################################