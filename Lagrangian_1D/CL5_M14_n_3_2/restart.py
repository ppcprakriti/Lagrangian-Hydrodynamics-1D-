# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:09:35 2017

@author: prakriti
"""

from __future__ import division
import numpy as np
import os
import parameters as para

i_dmp = para.i_t_dmp

def restart():
    
    flag = True
    t=0
    while (flag==True):
        if os.path.exists("./data/t_%d"%(i_dmp*t)+".dat"):
            t = t+1
        else:
            flag=False
            
            
    data = np.loadtxt("./data/t_%d"%(i_dmp*(t-1))+".dat")
    tme = data[0,0]
    L = data[1,0]
    dt = data[2,0]
    Ncool = data[3,0]
    r = data[:,1]
    v = data[:,2]
    rho = data[:,3]
    prs = data[:,4]
    u = data[:,5]
    K = data[:,6]
    T = data[:,7]
    dm = data[:,8]
    N = np.size(data[:,1])
    return tme,r,v,rho,prs,u,K,T,dm, N,(t-1),L,dt,Ncool

    
