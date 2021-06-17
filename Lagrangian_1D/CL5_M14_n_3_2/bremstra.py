# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 16:29:34 2017

@author: prakriti
"""

from __future__ import division
from numba import jit
import numpy as np
import mypath as mypath
import sys
from scipy.integrate import quad
sys.path.insert(0,mypath.path1)
import constants as cons
sys.path.insert(0,mypath.path2)
import f_x as fx
import compton_ht as ch

@jit
def F_p(y,t_p,t_n):
    return np.exp(-y)*fx.KN_cs(y*t_p)*fx.E_param(y*t_p,t_n)
#F_pn = np.vectorize(F_p)

@jit    
def F(t_p,t_n):
    I = quad(F_p,0.0,np.inf,(t_p,t_n), limit=1000)[0]
    return I 
#F_n = np.vectorize(F)        
@jit
def E_dot_br(r,rho,T):
    ne = rho/(cons.mue*cons.mp)
    n = rho/(cons.mu*cons.mp)
    constant = (2**(11/2))*((np.pi/3)**(3./2))*(cons.Z**2)*(cons.qe**6)/(cons.h_pl*cons.me*cons.c*cons.c)
    ni=n-ne
    return constant*ne*ni*(cons.kB*T/(cons.me*cons.c*cons.c))**0.5 
    
@jit
def Intr(r1,rho1,T1):
    N = np.size(r1)
    M = np.zeros(np.size(r1))
    dr = np.zeros(np.size(r1))
    t1 = np.zeros(np.size(r1))
    t1[:]=cons.kB*T1[:]/(cons.me*cons.c*cons.c)
    t2 = cons.kB*T1[N-1]/(cons.me*cons.c*cons.c)
    dr[0:-1] = r1[1:] - r1[0:-1]
    for i in np.arange(np.size(r1)-1):
        M[i] = r1[i]*r1[i]*E_dot_br(r1[i],rho1[i],T1[i])*F(t1[i],t2)*dr[i]
    #M[N-1] = r1[N-1]*r1[N-1]*E_dot_br(r1[N-1],rho1[N-1],T1[N-1])*F(t1[N-1],t2)*dr[N-2]
    return np.sum(M),np.size(dr)

@jit   
def dE_br(r,rho,T):
    ne = np.zeros(np.size(r))
    Tot = np.zeros(np.size(r))
    ne[:] = rho[:]/(cons.mp*cons.mue)
    for i in np.arange(1,len(r)+1):
        r1=r[0:i]
        rho1=rho[0:i]
        T1=T[0:i]
        #Tot[i-1] = Intr(r1,rho1,T1)[0]
        Tot[i-1] = (8./3.)*(1.5*cons.kB*T[i-1])*ne[i-1]*ch.Sigma_T*Intr(r1,rho1,T1)[0]/(r[i-1]*r[i-1]*cons.me*cons.c*cons.c) 
    return Tot
        
 #-(8./3.)(1.5*cons.kB*T[i-1])*ne[i-1]*ch.Sigma_T*Intr(r1,rho1,T1)/(r[i-1]*r[i-1]*cons.me*cons.c*cons.c)   
    
#r = np.zeros(5)
#r[0] = 3.086e21
#r[1] = 6.1e21
#r[2] = 1.e22
#r[3] = 2.e22
#r[4] = 3.e22
#rho = np.zeros(5)
#rho[0] = 5.e-22
#rho[1] = 2.e-22
#rho[2] = 1.e-25
#rho[3] = 5.e-26
#rho[4] = 2.e-26
#T = np.zeros(5)
#T[0] = 1.e4
#T[1] = 8.e3
#T[2] = 9.e3
#T[3] = 9.1e6
#T[4] = 2.e7
#print dE_br(r,rho,T)
