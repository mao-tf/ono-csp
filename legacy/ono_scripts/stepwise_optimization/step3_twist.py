import os
import pandas as pd
import time
from make_step3_twist import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse
import numpy as np

def main_process(args):
    auto_dir = args.auto_dir
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    auto_csv_path = os.path.join(auto_dir,'step3_twist.csv')
    if not os.path.exists(auto_csv_path):        
        df_E = pd.DataFrame(columns = ['cx','cy','cz','a','b','alpha','Rt','A2','E','E_i01','E_ip1','E_ip2','E_it1','E_it2','E_it3','E_it4','status','file_name'])##いじる
        df_E.to_csv(auto_csv_path,index=False)

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.isTest)
        time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,isTest):
    num_init = args.num_init
    fixed_param_keys = ['a','b','alpha','Rt','A2']
    opt_param_keys = ['cx','cy','cz']
    
    auto_csv = os.path.join(auto_dir,'step3_twist.csv')
    df_E = pd.read_csv(auto_csv)
    df_queue = df_E.loc[df_E['status']=='InProgress',['file_name']]
    len_queue = len(df_queue)
    
    for idx,row in zip(df_queue.index,df_queue.values):
        file_name = row[0]
        log_filepath = os.path.join(*[auto_dir,'gaussian',file_name])
        if not(os.path.exists(log_filepath)):
            continue
        E_list=get_E(log_filepath)
        if len(E_list)!=5:
            continue
        else:
            len_queue-=1
            Ei01=float(E_list[0]);Eip1=float(E_list[1]);Eip2=float(E_list[2]);Eit1=float(E_list[3]);Eit2=float(E_list[4]);Eit3=float(E_list[4]);Eit4=float(E_list[3])##ここも計算する分子数に合わせて調整
            E = Ei01+Eip1+Eip2+Eit1+Eit2+Eit3+Eit4 ## inter action energy of adjacent 7 molecular pair 
            #### TODO
            df_E.loc[idx, ['E_i01','E_ip1','E_ip2','E_it1','E_it2','E_it3','E_it4','E','status']] = [Ei01,Eip1,Eip2,Eit1,Eit2,Eit3,Eit4,E,'Done']
            df_E.to_csv(auto_csv,index=False)
            break
    isAvailable = len_queue < num_nodes 
    if isAvailable:
        dict_matrix = get_params_dict(auto_dir,num_init, fixed_param_keys, opt_param_keys)
        if len(dict_matrix)!=0:
            for i in range(len(dict_matrix)):
                params_dict=dict_matrix[i]
                alreadyCalculated = check_calc_status(auto_dir,params_dict)
                if not(alreadyCalculated):
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict},isTest=isTest)
                    df_newline = pd.Series({**params_dict,'E':0.,'E_i01':0.,'E_ip1':0.,'E_ip2':0.,'E_it1':0.,'E_it2':0.,'E_it3':0.,'E_it4':0.,'status':'InProgress','file_name':file_name})
                    df_E=pd.concat([df_E,pd.DataFrame([df_newline])],ignore_index=True)
                    df_E.to_csv(auto_csv,index=False)
    
    init_params_csv=os.path.join(auto_dir, 'step3_twist_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params_done = filter_df(df_init_params,{'status':'Done'})
    isOver = True if len(df_init_params_done)==len(df_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    df_E= pd.read_csv(os.path.join(auto_dir,'step3_twist.csv'))
    if len(df_E)==0:
        return False
    df_E_filtered = filter_df(df_E, params_dict)
    df_E_filtered = df_E_filtered.reset_index(drop=True)
    try:
        status = get_values_from_df(df_E_filtered,0,'status')
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_init,fixed_param_keys,opt_param_keys):
    init_params_csv=os.path.join(auto_dir, 'step3_twist_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step3_twist.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='InProgress']

    if len(df_init_params_inprogress) < num_init:
        df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
        for index in df_init_params_notyet.index:
            df_init_params = update_value_in_df(df_init_params,index,'status','InProgress')
            df_init_params.to_csv(init_params_csv,index=False)
            params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
            return [params_dict]
    dict_matrix=[]
    for index in df_init_params_inprogress.index:
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'cx':opt_params_matrix[0][0],'cy':opt_params_matrix[0][1],'cz':opt_params_matrix[0][2]}
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            if np.max(df_init_params.index) < index+1:
                status = 'Done'
            else:
                status = get_values_from_df(df_init_params,index+1,'status')
            df_init_params.to_csv(init_params_csv,index=False)
            
            if status=='NotYet':           
                opt_params_dict = get_values_from_df(df_init_params,index+1,opt_param_keys)
                df_init_params = update_value_in_df(df_init_params,index+1,'status','InProgress')
                df_init_params.to_csv(init_params_csv,index=False)
                dict_matrix.append({**fixed_params_dict,**opt_params_dict})
            else:
                continue

        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict={'cx':opt_params_matrix[i][0],'cy':opt_params_matrix[i][1],'cz':opt_params_matrix[i][2]}
                df_inprogress = filter_df(df_cur, {**fixed_params_dict,**opt_params_dict,'status':'InProgress'})
                if len(df_inprogress)>=1:
                    continue
                else:
                    d={**fixed_params_dict,**opt_params_dict}
                    dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    cx_init_prev = init_params_dict['cx']; cy_init_prev = init_params_dict['cy']; cz_init_prev = init_params_dict['cz']
    a = init_params_dict['a'];b = init_params_dict['b'];Rt = init_params_dict['Rt'];A2 = init_params_dict['A2']
    
    while True:
        E_list=[];heri_list=[]
        para_list=[]
        for cx in [cx_init_prev]:## cx is kept to 0 because of glide symmetry 
            for cy in [cy_init_prev-0.1,cy_init_prev,cy_init_prev+0.1]:
                for cz in [cz_init_prev-0.1,cz_init_prev,cz_init_prev+0.1]:
                    cx = np.round(cx,1);cy = np.round(cy,1);cz = np.round(cz,1)
                    df_val_ab = df_val[
                        (df_val['cx']==cx)&(df_val['cy']==cy)&(df_val['cz']==cz)
                        &(df_val['Rt']==Rt)&(df_val['A2']==A2)&(df_val['a']==a)&(df_val['b']==b)&
                        (df_val['status']=='Done')]
                    if len(df_val_ab)==0:
                        para_list.append([cx,cy,cz])
                        continue
                    heri_list.append([cx,cy,cz]);E_list.append(df_val_ab['E'].values[0])
        if len(para_list) != 0:
            return False,para_list
        cx_init,cy_init,cz_init = heri_list[np.argmin(np.array(E_list))]
        if cx_init==cx_init_prev and cy_init==cy_init_prev and cz_init==cz_init_prev:
            return True,[[cx_init,cy_init,cz_init]]
        else:
            cx_init_prev=cx_init;cy_init_prev=cy_init;cz_init_prev=cz_init

def get_values_from_df(df,index,key):
    return df.loc[index,key]

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

def filter_df(df, dict_filter):
    for k, v in dict_filter.items():
        if type(v)==str:
            df=df[df[k]==v]
        else:
            df=df[df[k]==v]
    df_filtered=df
    return df_filtered

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--init',action='store_true')
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    parser.add_argument('--num-init',type=int,help='number of parameters in progress at init_params.csv')
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    print("----finish process----")
    