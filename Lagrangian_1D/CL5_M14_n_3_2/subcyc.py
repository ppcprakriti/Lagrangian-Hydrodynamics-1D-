# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:31:58 2017

@author: prakriti
"""

import sys
import mypath as mypath
import numpy as np
from numba import jit
sys.path.insert(0,mypath.path1)
import constants as cons
import WRT_t as wrtt
import cooling as cl
sys.path.insert(0,mypath.path2)
import compton_ht as ch
from scipy.integrate import quad
import L_BH_comp as Lbh
import find_mindt as fmd
import heating as ht
import cl_time as ct


@jit
def subcyc(r,v,u,rho,T,t_i,t_f,dt,L,M_solm0,Ncl):
    dtc = fmd.find_min_dtc(r,v,rho,T,u,t_i,dt,L,M_solm0,Ncl)
    dtc = min(dtc,dt)
#    Nt = (t_f-t_i)/dtc
    tt=t_i 
#    for tt in np.arange(1,Nt):
#        t = t_i + tt*dtc
    n = np.zeros(len(r))
    ne = np.zeros(len(r))
    ni = np.zeros(len(r))
    p = np.zeros(len(r))
    n[:]=rho[:]/(cons.mp*cons.mu)
    ne[:] = rho[:]/(cons.mp*cons.mue)
    ni[:] = n[:]-ne[:]
    m = 0
    while (tt<t_f):
        m = m+1
        tt= tt+dtc
        rvir = wrtt.r200(tt,M_solm0)
        u[np.where(r<=rvir)] = u[np.where(r<=rvir)]/(1.+((ne[np.where(r<=rvir)]*ni[np.where(r<=rvir)]*cl.lam_2(T[np.where(r<=rvir)]) - (ht.H(r,v,rho,T,tt,dtc,L)[0][np.where(r<=rvir)]))*dtc/(rho[np.where(r<=rvir)]*u[np.where(r<=rvir)]) ))
        p[:] = (cons.gamma - 1)*rho[:]*u[:]
        T[:] = p[:]*cons.mu*cons.mp/(cons.kB*rho[:])
        #dtc = fmd.find_min_dt(r,v,rho,T,u,tt,dtc,L,M_solm0,Ncl)[5]
 
    return u[:],p[:],T[:],dtc, m
    
