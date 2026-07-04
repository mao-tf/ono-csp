import numpy as np
from make_step1 import get_monomer_xyzR

def vdw_R(A2,A3,theta,dimer_mode,monomer_name):##2分子が接しないようにするにはtheta方向にどれだけずらせばいいのか
    monomer_1=get_monomer_xyzR(monomer_name,Ta=0.,Tb=0.,Tc=0.,A2=A2,A3=A3)
    if dimer_mode=='t':
        monomer_2=get_monomer_xyzR(monomer_name,Ta=0.,Tb=0.,Tc=0.,A2=A2,A3=-A3)#convertor(monomer,-A2,-A3+glide)
    elif dimer_mode=='a' or dimer_mode=='b':
        monomer_2=get_monomer_xyzR(monomer_name,Ta=0.,Tb=0.,Tc=0.,A2=A2,A3=A3)
    R_clps=0
    for x1,y1,z1,rad1 in monomer_1:
        for x2,y2,z2,rad2 in monomer_2:
            eR=np.array([np.cos(np.radians(theta)),np.sin(np.radians(theta)),0.0])
            R_1b=np.dot(eR,np.array([x1,y1,z1]))
            R_2b=np.dot(eR,np.array([x2,y2,z2]))
            R_12=np.array([x2-x1,y2-y1,z2-z1])
            R_12b=np.dot(eR,R_12)
            R_12a=np.linalg.norm(R_12-R_12b*eR)
            if (rad1+rad2)**2-R_12a**2<0:
                continue
            else:
                R_clps=max(R_clps,R_1b-R_2b+np.sqrt((rad1+rad2)**2-R_12a**2))
    return R_clps
    