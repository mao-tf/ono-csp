import numpy as np
import os
import pandas as pd
from utils import Rod, R2atom
import subprocess

MONOMER_LIST = ["BTBT","naphthalene","anthracene","tetracene","pentacene","hexacene","demo"]

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

def make_gjf_xyz(auto_dir,monomer_name,params_dict):## atomic coordinates for 41 * 2 = 82 pairs to be calculate
    a_ = params_dict['a']; b_ = params_dict['b']
    A2 = 0; A3 = params_dict['theta']
    gij_xyz_lines1 = ['$ RunGauss\n']
    gij_xyz_lines2 = ['$ RunGauss\n']
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    z_list=[np.round(z,1) for z in np.linspace(np.round(0,1),np.round(4,1),int(np.round(np.round(4,1)/0.1))+1)]
    ##molecular shift along molecular long axis
    for z in z_list:
        if a_>b_:## slipped parallel contact
            monomer_array_p1 = get_monomer_xyzR(monomer_name,0,b_,z,A2,A3)
        else:
            monomer_array_p1 = get_monomer_xyzR(monomer_name,a_,0,z,A2,A3)
    
        monomer_array_t1 = get_monomer_xyzR(monomer_name,a_/2,b_/2,z,A2,-A3) ##T-shape contact

        dimer_array_t1 = np.concatenate([monomer_array_i,monomer_array_t1])
        dimer_array_p1 = np.concatenate([monomer_array_i,monomer_array_p1])
    
        file_description = '{}_A2={}_A3={}'.format(monomer_name,int(A2),round(A3,2))##
        line_list_dimer_p1 = get_xyzR_lines(dimer_array_p1,file_description+'_p1')
        line_list_dimer_t1 = get_xyzR_lines(dimer_array_t1,file_description+'_t1')
        
        gij_xyz_lines1 = gij_xyz_lines1 + line_list_dimer_t1 + ['\n\n--Link1--\n']##T-shpe
        gij_xyz_lines2 = gij_xyz_lines2 + line_list_dimer_p1 + ['\n\n--Link1--\n']##slipped parallel
    
    gij_xyz_lines1 = gij_xyz_lines1 + ['\n\n\n']
    gij_xyz_lines2 = gij_xyz_lines2 + ['\n\n\n']
    file_base_name = ''
    file_base_name += monomer_name
    file_base_name += '_step2_'
    file_base_name += 'a={}_b={}_theta={}_'.format(a_,b_,A3)
    file_name1 = file_base_name
    file_name2 = file_base_name
    file_name1 +='1.inp'##T-shape
    file_name2 +='2.inp'##slipped parallel
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)
    gij_xyz_path1 = os.path.join(auto_dir,'gaussian',file_name1)
    gij_xyz_path2 = os.path.join(auto_dir,'gaussian',file_name2)
    with open(gij_xyz_path1,'w') as f1:
        f1.writelines(gij_xyz_lines1)
    with open(gij_xyz_path2,'w') as f2: 
        f2.writelines(gij_xyz_lines2)
    return file_base_name

def get_one_exe(file_basename):## batch job scripts
    
    cc_list=[
        ## you can use your environment for Gaussian16 calculations
        '\n',
        'g16 < {}.inp > {}.log \n'.format(file_basename,file_basename),
        '\n'
        ## you can use your environment for Gaussian16 calculations
            ]

    return cc_list

def exec_gjf(auto_dir, monomer_name, params_dict,isTest):
    inp_dir = os.path.join(auto_dir,'gaussian')
    print(params_dict)

    file_base_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)

    file_basename1 = file_base_name
    file_basename2 = file_base_name
    file_basename1 +='1'##T-shaped
    file_basename2 +='2'##slipped parallel

    cc_list1 = get_one_exe(file_basename1)
    sh_filename1 = file_basename1 + '.r1'
    sh_path1 = os.path.join(inp_dir,sh_filename1)
    with open(sh_path1,'w') as f:
        f.writelines(cc_list1)
    if not(isTest):## subission of the batch job fot T-shape calculations
        subprocess.run(['pjsub',sh_path1])
    log_file_name1 = file_basename1 + '.log'

    cc_list2 = get_one_exe(file_basename2)
    sh_filename2 = file_basename2 + '.r1'
    sh_path2 = os.path.join(inp_dir,sh_filename2)
    with open(sh_path2,'w') as f:
        f.writelines(cc_list2)
    if not(isTest):## subission of the batch job fot slipped parallel calculations
        subprocess.run(['pjsub',sh_path2])
    log_file_name2 = file_basename2 + '.log'
    return log_file_name1,log_file_name2
    
############################################################################################