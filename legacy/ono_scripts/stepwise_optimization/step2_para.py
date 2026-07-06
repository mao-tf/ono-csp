import os
import pandas as pd
import time
from utils import get_E
import argparse
import numpy as np
from make_step2_para import exec_gjf
import math
import cmath
import matplotlib.pyplot as plt

def submit_process(args):## submission of the energy calculations jobs
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    isTest= args.isTest
    isMain= args.isMain
    if isMain:
        return
    isEnd= args.isEnd
    if isEnd:
        return
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)
    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    log_file= exec_gjf(auto_dir, monomer_name, params_dict,isTest)
    time.sleep(2)
    print(log_file)


def main_process(args):
    isEnd= args.isEnd
    if isEnd:
        return
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)

    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    file_base_name = get_file_base_name(monomer_name,params_dict)
    file_name_1 = file_base_name
    file_name_2 = file_base_name
    file_name_1 += '1.log'
    file_name_2 += '2.log'##log file name
    log_filepath_1 = os.path.join(*[auto_dir,'gaussian',file_name_1])
    log_filepath_2 = os.path.join(*[auto_dir,'gaussian',file_name_2])

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(log_filepath_1,log_filepath_2)
        time.sleep(1)

def listen(log_filepath_1,log_filepath_2):## checking whether the energy calculatios are finished or not
    E_list1=get_E(log_filepath_1)
    E_list2=get_E(log_filepath_2)
    if len(E_list1)!=41 or len(E_list2)!=41:## 41 molecular pairs are calculated in each file 
        isOver =False
    else:
        isOver=True
    return isOver

def get_file_base_name(monomer_name,params_dict):
    a_ = params_dict['a']; b_ = params_dict['b']; theta = params_dict['theta']
    file_base_name = ''
    file_base_name += monomer_name
    file_base_name += '_step2_'
    file_base_name += 'a={}_b={}_theta={}_'.format(a_,b_,theta)
    return file_base_name

def end_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    file_base_name = get_file_base_name(monomer_name,params_dict)
    file_name_1 = file_base_name
    file_name_2 = file_base_name
    file_name_1 += '1.log'
    file_name_2 += '2.log'
    while True:
        log_filepath_1 = os.path.join(*[auto_dir,'gaussian',file_name_1])
        log_filepath_2 = os.path.join(*[auto_dir,'gaussian',file_name_2])
        E_list1=get_E(log_filepath_1)##t
        E_list2=get_E(log_filepath_2)##p
        if len(E_list1)==41 & len(E_list2)==41:
            break
        time.sleep(60)
    E_list_1=[]
    E_list_2=[]
    for i in range(len(E_list1)):
        E_list_1.append(E_list1[len(E_list1)-i-1])
        E_list_2.append(E_list2[len(E_list2)-i-1])
    for i in range(len(E_list1)-1):
        E_list_1.append(E_list1[i+1])
        E_list_2.append(E_list2[i+1])
    z_list=[np.round(z,1) for z in np.linspace(-np.round(4,1),np.round(4,1),int(np.round(np.round(8,1)/0.1))+1)]
    para_list=[]
    for E1,E2,z in zip (E_list_1,E_list_2,z_list):
        para_list.append(z,E1,E2)
    df=pd.DataFrame(para_list,columns=['z','Et','Ep'])
    result_params_csv=os.path.join(auto_dir, 'step2_para.csv')
    df.to_csv(result_params_csv,index=False)
    
def plot2d(args):
    auto_dir = args.auto_dir
    result_params_csv=os.path.join(auto_dir, 'step2_para.csv')
    df=pd.read_csv(result_params_csv)
    E_1d_list=[];z_list_t=[];z_list_p=[]
    z_list=[np.round(z,1) for z in np.linspace(-np.round(4,1),np.round(4,1),int(np.round(np.round(8,1)/0.1))+1)]
    for z1 in z_list:##t
        for z2 in z_list:##p            
            if ((z1-z2)>4.01) or (-4.01>(z1-z2)):
                continue
            else:
                z_list_t.append(z1);z_list_p.append(z2)
                Et1 = df[df['z']==z1]['Et'].values.tolist()[0]
                Et2 = df[df['z']==round(z1-z2,1)]['Et'].values.tolist()[0]
                Ep = df[df['z']==z2]['Ep'].values.tolist()[0]
                E=2*(Et1+Et2+Ep)
            E_1d_list.append(E)

    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    a = params_dict['a']; b = params_dict['b']

    x_list=[];y_list=[];para_list=[]
    for zt,zp,E in zip(z_list_t,z_list_p,E_1d_list):
        za=2*zt-zp;zb=zp
        Z=1/math.sqrt(1+(za/a)**2+(zb/b)**2)
        theta=np.degrees(math.acos(Z))
        phi=cmath.phase(za/a+zb/b*1j)
        x_list.append(theta*math.cos(phi))
        y_list.append(theta*math.sin(phi))
        para_list.append([])

    plt.tick_params(labelsize=16)
    plt.gca().set_aspect('equal')
    plt.scatter(x_list,y_list,s=11,c=E_1d_list,cmap='coolwarm')
    plt.xticks([-45,-30,-15,0,15,30,45])
    plt.yticks([-30,-20,-10,0,10,20,30])
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    plt.savefig(os.path.join(auto_dir,'result.png'),dpi=500)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    #parser.add_argument('--init',action='store_true')
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--isEnd',action='store_true')
    parser.add_argument('--isMain',action='store_true')
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    args = parser.parse_args()

    
    print("----main process----")
    submit_process(args)
    main_process(args)
    end_process(args)
    if args.plot:
        plot2d(args)
    print("----finish process----")
    