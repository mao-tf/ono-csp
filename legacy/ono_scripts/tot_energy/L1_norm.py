import numpy as np

def L1norm(N):
    para_list=[]
    for i in range(N+1):
        if i == 0:
            para_list.append([0,N,'p'])
            if N>0:
                para_list.append([0,-N,'p'])
        elif i == N:
            para_list.append([N,0,'p'])
            para_list.append([-N,0,'p'])
        else:
            M=abs(N-i)
            para_list.append([i,M,'p']);para_list.append([-i,M,'p'])
            para_list.append([i,-M,'p']);para_list.append([-i,-M,'p'])
    for i in range(N):
        X=round(i+0.5,1)
        Y=round(N-i-0.5,1)
        para_list.append([X,Y,'t']);para_list.append([-X,Y,'t'])
        para_list.append([X,-Y,'t']);para_list.append([-X,-Y,'t'])
    return(para_list)

def L1norm2(N):
    para_list=[]
    for i in range(N):
        para_list.append([i-N,i,'p'])
    for i in range(N,2*N):
        para_list.append([i-N,2*N-i,'p'])
    for i in range(N):
        X=round(i+0.5,1)
        Y=round(N-i-0.5,1)
        para_list.append([X,Y,'t'])
        para_list.append([-Y,X,'t'])
    return(para_list)