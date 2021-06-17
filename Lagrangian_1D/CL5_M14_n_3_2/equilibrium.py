# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:44:42 2017

@author: prakriti
"""

from __future__ import division
import numpy as np
import mypath as mypath
import sys
from scipy import integrate
import matplotlib.pyplot as plt
sys.path.insert(0,mypath.path1)
import constants as cons
import WRT_t as wrtt
import initial as ini
import make_initial as mi


C = cons.kB*cons.fT/(cons.mu*cons.mp*(cons.mue*cons.mp)**(cons.gamma-1))

#def K(r):
#    N = np.size(r)
#    K = np.zeros(N)
#    K0=10.0
#    K100=0.0
#    alpha=1.1
#    for i in np.arange(N):
#        K[i] = K0 + K100*(r[i]/(100*cons.Kpc))**alpha
#    return K

def rho_hse(r,K,t,M_solm0,rho_out):
    N = np.size(r)
    r_rev = r[::-1]
    K_rev = K[::-1]
    integr = np.zeros(N)
    rho = np.zeros(N)
    for i in np.arange(N):
        integr[i] = wrtt.g(r_rev[i],t,M_solm0)*(K_rev[i]**(-1./cons.gamma))
        
    y = -integrate.cumtrapz(integr,r_rev,initial=0)
    for i in np.arange(N):
        rho[i] = (((rho_out**(cons.gamma-1))*((K_rev[0]/K_rev[i])**((cons.gamma-1)/cons.gamma))) + (((cons.gamma-1)/(C*cons.gamma))*y[i]/(K_rev[i]**((cons.gamma-1)/cons.gamma))) )**(1./(cons.gamma-1))
    return rho[::-1]
    
#def T_hse(r,K,t,M_solm0,rho_out):
#    N = np.size(r)
#    T = np.zeros(N)
#    for i in np.arange(N):
#        T[i] = cons.fT*K[i]*(rho_hse(r,K,t,M_solm0,rho_out)[i]/(cons.mue*cons.mp))**(cons.gamma-1)
#    return T
#
#def prs_hse(r,K,t,M_solm0,rho_out):
#    N = np.size(r)
#    prs = np.zeros(N)
#    for i in np.arange(N):
#        prs[i] = (rho_hse(r,K,t,M_solm0,rho_out)[i]/(cons.mu*cons.mp))*cons.kB*T_hse(r,K,t,M_solm0,rho_out)[i]
#    return prs
#    
#def u_hse(r,K,t,M_solm0,rho_out):
#    N = np.size(r)
#    u = np.zeros(N)
#    for i in np.arange(N):
#        u[i] = prs_hse(r,K,t,M_solm0,rho_out)[i]/((cons.gamma-1)*rho_hse(r,K,t,M_solm0,rho_out)[i])
#    return u
    
def dm_hse(r,rho):
    N = np.size(r)
    dm = np.zeros(N)
    for i in np.arange(N-1):
        dm[i] = rho[i]*((4./3.)*np.pi)*(r[i+1]**3 - r[i]**3)
    dm[N-1] = dm[N-2]
    return dm
  
        
if __name__ == '__main__':   
    t=ini.t0
    M_solm0 = 5.e13
    rho_out = 32.5*wrtt.crit_rho(t)
    #Nri = 50
    #r = mi.ri_half1(1.0*cons.Kpc,1.1,Nri)[1]
    r = np.linspace(1*cons.Kpc,wrtt.r200(t,M_solm0),10)
    #rho = np.zeros(20)
    #for i in np.arange(20):
    #    rho[i] = rho_hse(r,K(r[:]),t,M_solm0,rho_out)[i]
    
