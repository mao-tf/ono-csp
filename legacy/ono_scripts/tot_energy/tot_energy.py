import numpy as np
import os
import pandas as pd
from utils import Rod, R2atom
from L1_norm import L1norm2,L1norm
import argparse
import subprocess

MONOMER_LIST = ["naphthalene","anthracene","tetracene","pentacene","hexacene"]

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):## rotation and translation operations
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('~/path/to/monomer/{}.csv'.format(monomer_name))
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
        
def translation(xyzR_array,vx,vy,vz):
    vec = np.array([vx,vy,vz])
    xyz_array = xyzR_array[:,:3]
    xyz_array = xyz_array + vec
    R_array = xyzR_array[:,3].reshape((-1,1))
    return np.concatenate([xyz_array,R_array],axis=1)
        
def mk_xyz(auto_dir,monomer_name,params_dict,I):## I = 1,2,3
    a_ = params_dict['a']; b_ = params_dict['b']; cx=params_dict['cx'];cy=params_dict['cy'];cz=params_dict['cz']
    A2 = params_dict['A2']; A3 = params_dict['theta']
    Rt = params_dict['Rt']; Rp = params_dict['Rp']
    
    monomer_array_c1 = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_c2 = get_monomer_xyzR(monomer_name,0,0,0,A2,-A3)
    a =np.array([a_,0,2*Rt-Rp]);b =np.array([0,b_,Rp])
    layer_list_list=[]
    for i in range(I+1):
        layer_list_list.append(L1norm2(i))
    layer_list_list_i=[]
    for i in range(I):
        if i == 0:
            layer_list=np.concatenate([L1norm(0),L1norm(1)])
            layer_list_list_i.append(layer_list)
        else:
            layer_list_list_i.append(L1norm(i+1))
    #print(len(layer_list))
    intra_array=monomer_array_c1    
    layer_list = layer_list_list[I]
    for x,y,arr in layer_list:
        x=float(x);y=float(y)
        v=x*a+y*b
        vx=v[0];vy=v[1];vz=v[2]
        if arr =='p':
            monomer_array=translation(monomer_array_c1,vx,vy,vz)
        if arr =='t':
            monomer_array=translation(monomer_array_c2,vx,vy,vz)
        intra_array=np.concatenate([intra_array,monomer_array])
        
    inter_array_1=monomer_array_c1
    for i in range(1,I+1):
        layer_list=layer_list_list_i[I-i]
        for x,y,arr in layer_list:
            x=float(x);y=float(y)
            v=x*a+y*b
            vx=v[0]+i*cx;vy=v[1]+i*cy;vz=v[2]+i*cz
            if arr =='p':
                monomer_array=translation(monomer_array_c1,vx,vy,vz)
            if arr =='t':
                monomer_array=translation(monomer_array_c2,vx,vy,vz)
            inter_array_1=np.concatenate([inter_array_1,monomer_array])
                
    inter_array_2=monomer_array_c2
    for i in range(1,I+1):
        layer_list=layer_list_list_i[I-i]
        for x,y,arr in layer_list:
            x=float(x);y=float(y)
            v=x*a+y*b
            vx=v[0]+i*cx;vy=v[1]+i*cy;vz=v[2]+i*cz
            if arr =='p':
                monomer_array=translation(monomer_array_c2,vx,vy,vz)
            if arr =='t':
                monomer_array=translation(monomer_array_c1,vx,vy,vz)
            inter_array_2=np.concatenate([inter_array_2,monomer_array])        
    
    os.makedirs(auto_dir,exist_ok=True)
    
    lines=['10000 \n',f'{monomer_name} \n']##4分子のxyzファイルを作成

    intra_path = os.path.join(auto_dir,f'test{I}_intra.xyz')
    with open(intra_path,'w') as f:
        f.writelines(lines)
        for x,y,z,R in intra_array:
            atom = R2atom(R)
            line = f'{atom} {x} {y} {z}\n'
            f.write(line)
    
    inter_path_1 = os.path.join(auto_dir,f'test{I}_inter_1.xyz')
    with open(inter_path_1,'w') as f:
        f.writelines(lines)
        for x,y,z,R in inter_array_1:
            atom = R2atom(R)
            line = f'{atom} {x} {y} {z}\n'
            f.write(line)
    
    inter_path_2 = os.path.join(auto_dir,f'test{I}_inter_2.xyz')
    with open(inter_path_2,'w') as f:
        f.writelines(lines)
        for x,y,z,R in inter_array_2:
            atom = R2atom(R)
            line = f'{atom} {x} {y} {z}\n'
            f.write(line)
    
def init_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    os.makedirs(auto_dir, exist_ok=True)
    params_path=os.path.join(auto_dir,'params.csv')
    df_init=pd.read_csv(params_path)
    params_dict = df_init.loc[0,['cx','cy','cz','a','b','theta','Rt','Rp','A2']].to_dict()
    for I in range(1,4):
        mk_xyz(auto_dir,monomer_name,params_dict,I)

def make_gjf_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    os.makedirs(auto_dir, exist_ok=True)
    df_mono=pd.read_csv('~/path/to/monomer/{}.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    mol=len(atoms_array_xyzR)

    lines1 = [
        '%mem=24GB\n',
        '%nproc=48\n',
        '#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2\n',
        '\n',
        'interaction calculation\n',
        '\n',
        '0 1 0 1 0 1\n'
    ]

    ## intralayer##
    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I)}_intra.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        N = int(len(o_1) / mol)
        for i in range(1,N):
            with open(os.path.join(auto_dir,f'{int(I)}_intra_{i}.inp'), 'w') as f1:
                f1.writelines(lines1)
                for k in range(mol):
                    f1.write('{} {} {} {} 1\n'.format(o_1[k][0], o_1[k][1], o_1[k][2], o_1[k][3]))
                for k in range(mol):
                    f1.write('{} {} {} {} 2\n'.format(o_1[mol * i + k][0], o_1[mol * i + k][1], o_1[mol * i + k][2],o_1[mol * i + k][3]))
                f1.write('\n\n\n')
    
    ## interlayer##
    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I)}_inter_1.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        N = int(len(o_1) / mol)
        for i in range(1,N):
            with open(os.path.join(auto_dir,f'{int(I)}_inter_1_{i}.inp'), 'w') as f1:
                f1.writelines(lines1)
                for k in range(mol):
                    f1.write('{} {} {} {} 1\n'.format(o_1[k][0], o_1[k][1], o_1[k][2], o_1[k][3]))
                for k in range(mol):
                    f1.write('{} {} {} {} 2\n'.format(o_1[mol * i + k][0], o_1[mol * i + k][1], o_1[mol * i + k][2],o_1[mol * i + k][3]))
                f1.write('\n\n\n')

    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I)}_inter_2.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        N = int(len(o_1) / mol)
        for i in range(1,N):
            with open(os.path.join(auto_dir,f'{int(I)}_inter_2_{i}.inp'), 'w') as f1:
                f1.writelines(lines1)
                for k in range(mol):
                    f1.write('{} {} {} {} 1\n'.format(o_1[k][0], o_1[k][1], o_1[k][2], o_1[k][3]))
                for k in range(mol):
                    f1.write('{} {} {} {} 2\n'.format(o_1[mol * i + k][0], o_1[mol * i + k][1], o_1[mol * i + k][2],o_1[mol * i + k][3]))
                f1.write('\n\n\n')

def job_sub_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    df_mono=pd.read_csv('~/path/to/monomer/{}.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    mol=len(atoms_array_xyzR)
    cc_list1=['batch job scripts and Gaussian source in your environment']
    cc_list3=['the end of batch job scripts']

    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I)}_inter_1.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        N = int(len(o_1) / mol)
        for i in range(1,N):
            with open(os.path.join(auto_dir,f'job_{int(I)}_inter_1_{i}.sh'),'w')as f1:
                f1.writelines(cc_list1)
                f1.write(f'g16 < {int(I)}_inter_1_{i}.inp > {int(I)}_inter_1_{i}.log \n\n')
                f1.writelines(cc_list3)
        for i in range(N-1):
            sh_path=os.path.join(auto_dir,f'job_{int(I)}_inter_1_{i}.sh')
            subprocess.run(['pjsub',sh_path])  

    for I in range(1,4):
        o_1 = []
        with open(f'test{int(I)}_inter_2.xyz') as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        N = int(len(o_1) / mol)
        for i in range(N-1):
            with open(os.path.join(auto_dir,f'job_{int(I)}_inter_2_{i}.sh'),'w')as f1:
                f1.writelines(cc_list1)
                f1.write(f'g16 < {int(I)}_inter_2_{i}.inp > {int(I)}_inter_2_{i}.log \n\n')
                f1.writelines(cc_list3)
        for i in range(1,N):
            sh_path=os.path.join(auto_dir,f'job_{int(I)}_inter_2_{i}.sh')
            subprocess.run(['pjsub',sh_path])  

def get_E(path_file):
    with open(path_file,'r') as f:
        lines=f.readlines()
    lines_E=[]
    for line in lines:
        if line.find('E(R')>-1 and len(line.split())>5:
            lines_E.append(float(line.split()[4])*627.510)
    E_list=[lines_E[5*i]-lines_E[5*i+1]-lines_E[5*i+2] for i in range(int(len(lines_E)/5))]
    return E_list[0]

def result_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    df_mono=pd.read_csv('~/path/to/monomer/{}.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    mol=len(atoms_array_xyzR)
    
    E_list=[]
    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I+1)}_intra.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        E_list_=[]
        N = int(len(o_1) / mol)
        for i in range(N - 1):
            j = int(i + 1)
            file=os.path.join(auto_dir,f'{int(I+1)}_intra_{j}.log')
            E=get_E(file)
            E_list_.append(E)
        E_list.append(E_list_)

    with open(os.path.join(auto_dir,f'result_intra_sub.txt','w'))as f:
        for e_list in E_list:
            for E in e_list:
                f.write(f'{E} ')
            f.write('\n')
    with open(os.path.join(auto_dir,f'result_intra.txt'),'w')as f:
        for e_list in E_list:
            E_tot=2*sum(e_list)
            f.write(f'{E_tot}')
            f.write('\n')

    E_list=[]
    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I+1)}_inter_1.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        E_list_=[]
        N = int(len(o_1) / mol)
        for i in range(N - 1):
            j = int(i + 1)
            file=os.path.join(auto_dir,f'{int(I+1)}_inter_1_{j}.log')
            E=get_E(file)
            E_list_.append(E)
        E_list.append(E_list_)

    with open(os.path.join(auto_dir,f'result_inter_sub_1.txt','w'))as f:
        for e_list in E_list:
            for E in e_list:
                f.write(f'{E} ')
            f.write('\n')
    with open(os.path.join(auto_dir,f'result_inter_1.txt'),'w')as f:
        for e_list in E_list:
            E_tot=2*sum(e_list)
            f.write(f'{E_tot}')
            f.write('\n')

    E_list=[]
    for I in range(1,4):
        o_1 = []
        with open(os.path.join(auto_dir,f'test{int(I+1)}_inter_2.xyz')) as f1:
            for line in f1:
                s = line.split()
                if len(s) > 3:
                    o_1.append([s[0], s[1], s[2], s[3]])
        E_list_=[]
        N = int(len(o_1) / mol)
        for i in range(N - 1):
            j = int(i + 1)
            file=os.path.join(auto_dir,f'{int(I+1)}_inter_2_{j}.log')
            E=get_E(file)
            E_list_.append(E)
        E_list.append(E_list_)

    with open(os.path.join(auto_dir,f'result_inter_sub_2.txt','w'))as f:
        for e_list in E_list:
            for E in e_list:
                f.write(f'{E} ')
            f.write('\n')
    with open(os.path.join(auto_dir,f'result_inter_2.txt'),'w')as f:
        for e_list in E_list:
            E_tot=2*sum(e_list)
            f.write(f'{E_tot}')
            f.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--init',action='store_true')
    parser.add_argument('--job-sub',action='store_true')
    parser.add_argument('--result',action='store_true')
    args = parser.parse_args()

    print("----main process----")
    if args.init:
        init_process(args)
        make_gjf_process(args)
    if args.job_sub:
        job_sub_process(args)
    if args.result_process(args):
        result_process(args)
    print("----finish process----")