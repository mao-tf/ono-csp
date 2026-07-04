##エネルギー等計算の対象を確認　R3:t-shaped R4:slipped-parallel
import numpy as np
from tqdm import tqdm
from make_step3_para import get_monomer_xyzR
from scipy.ndimage.filters import maximum_filter

def get_c_vec_vdw(monomer_name,R3,R4,a_,b_,theta):#,name_csv
    
    i=np.zeros(3); a=np.array([a_,0,2*R3-R4]); b=np.array([0,b_,R4]); t1=(a+b)/2;t2=(a-b)/2 ##ずらし方の定義
    monomer_array_i = get_monomer_xyzR(monomer_name,0.,0.,0.,0.,0.,theta)
    
    monomer_array_t = get_monomer_xyzR(monomer_name,0.,0.,0.,0.,0.,-theta)
    
    arr_list1=[[i,'p'],[b,'p'],[-b,'p'],[a,'p'],[-a,'p'],[t1,'t'],[-t1,'t'],[t2,'t'],[-t2,'t']]##層内の分子の座標とp or t
    arr_list2=[[i,'t'],[b,'t'],[-b,'t'],[a,'t'],[-a,'t'],[t1,'p'],[-t1,'p'],[t2,'p'],[-t2,'p']]##層内の分子の座標とp or t
    Rb_list=[np.round(Rb,1) for Rb in np.linspace(-np.round(b_/2,1),np.round(b_/2,1),int(np.round(2*np.round(b_/2,1)/0.1))+1)]
    Ra_list=[np.round(Ra,1) for Ra in np.linspace(-np.round(a_/2,1),np.round(a_/2,1),int(np.round(2*np.round(a_/2,1)/0.1))+1)]
    z_list=[];V_list=[];z_2dlist=[]
    for Ra in tqdm(Ra_list):
        z_maxlist=[]
        for Rb in Rb_list:
            z_max1=0
            for R,arr in arr_list1:
                if arr=='t':
                    monomer_array1=monomer_array_t
                elif arr=='p':
                    monomer_array1=monomer_array_i
                for x1,y1,z1,R1 in monomer_array1:#層内
                    x1,y1,z1=np.array([x1,y1,z1])+R
                    for x2,y2,z2,R2 in monomer_array_i:#i0
                        x2+=Ra
                        y2+=Rb
                        z2+=(2*R3-R4)*Ra/a_+R4*Rb/b_
                        z_sq=(R1+R2)**2-(x1-x2)**2-(y1-y2)**2
                        if z_sq<0:
                            z_clps=0.0
                        else:
                            z_clps=np.sqrt(z_sq)+z1-z2
                        z_max1=max(z_max1,z_clps)
            z_max2=0
            for R,arr in arr_list2:
                if arr=='t':
                    monomer_array1=monomer_array_t
                elif arr=='p':
                    monomer_array1=monomer_array_i
                for x1,y1,z1,R1 in monomer_array1:#層内
                    x1,y1,z1=np.array([x1,y1,z1])+R
                    for x2,y2,z2,R2 in monomer_array_t:#i0
                        x2+=Ra
                        y2+=Rb
                        z2+=(2*R3-R4)*Ra/a_+R4*Rb/b_
                        z_sq=(R1+R2)**2-(x1-x2)**2-(y1-y2)**2
                        if z_sq<0:
                            z_clps=0.0
                        else:
                            z_clps=np.sqrt(z_sq)+z1-z2
                        z_max2=max(z_max2,z_clps)
            z_max=max(z_max1,z_max2)
            z_list.append(z_max)
            V_list.append(a_*b_*z_max)
            z_maxlist.append(-z_max)
        z_2dlist.append(z_maxlist)
    return z_2dlist##ここでは層間距離を計算している

def detect_peaks(image, filter_size):
    local_max = maximum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant')
    detected_peaks = np.ma.array(image, mask=~(image == local_max))
    return detected_peaks

    

