# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 16:34:00 2017

@author: prakriti
"""

from __future__ import division
import numpy as np
import mypath as mypath
import sys
sys.path.insert(0,mypath.path2)
import heating as ht
sys.path.insert(0,mypath.path1)
import constants as cons
import WRT_t as wrtt
import cooling as cl

def find_min_dt(r,v,rho,T,u,t,dt,L,M_solm0,Ncl):
    cd = .1
    c_C = .3
    cv = 0.1
    c_c = 0.3
    n = rho/(cons.mp*cons.mu)
    ne = rho/(cons.mp*cons.mue) 
    ni = n - ne
    N = np.size(r)
    dt_dyn = np.zeros(N-1)
    dt_Cou = np.zeros(N-1)
    dt_v = np.zeros(N-1)
    dt_c = np.zeros(N)
#    for i in np.arange(N-1):
#        dt_dyn[i] = cd*(2.*(r[i+1] - (0.5*(r[i+1]-r[i])))/wrtt.g((r[i+1]- (0.5*(r[i+1]-r[i]))),t,M_solm0))**0.5
#        dt_Cou[i] = c_C*((r[i+1]-r[i])/(cons.gamma*(cons.gamma - 1.)*u[i])**0.5)
#        dt_v[i] = cv*np.abs((r[i+1]-r[i])/(v[i+1]-v[i]))
    dt_dyn[:] = cd*(2.*(r[1:] - (0.5*(r[1:] - r[0:-1])))/wrtt.g((r[1:]),t,M_solm0))**0.5
    dt_Cou[:] = c_C*((r[1:]-r[0:-1])/(cons.gamma*(cons.gamma - 1.)*u[0:-1])**0.5)
    dt_v[:] = cv*np.abs((r[1:]-r[0:-1])/(v[1:]-v[0:-1]))
    dt_c[:] = c_c*np.abs(u[:]*rho[:]/((ne[:]*ni[:]*cl.lam_2(T[:])))) 
    #ht.H(r,v,rho,T,t,dt,L)[0][:]
    #dt_dyn[N-2]=cd*(2.*(r[N-1]-r[N-2])/wrtt.g(r[N-1],t,M_solm0))**0.5 
    #dt_Cou[N-2]=c_C*((r[N-1]-r[N-2])/(cons.gamma*(cons.gamma - 1.)*u[N-1])**0.5)
    dt_d = np.min(dt_dyn[Ncl:])
    dt_C = np.min(dt_Cou[Ncl:])
    dtv = np.min(dt_v[Ncl:])
    dtc = np.min(dt_c[Ncl:])
    return dt_dyn,dt_Cou,dt_v,dt_c,min(dt_d,dt_C,dtv),min(dt_d,dt_C,dtv,dtc)



def find_min_dtc(r,v,rho,T,u,t,dt,L,M_solm0,Ncl):
    cd = .1
    c_C = .3
    cv = 0.1
    c_c = 0.3
    n = rho/(cons.mp*cons.mu)
    ne = rho/(cons.mp*cons.mue)
    ni = n - ne
    N = np.size(r)
    dt_dyn = np.zeros(N-1)
    dt_Cou = np.zeros(N-1)
    dt_v = np.zeros(N-1)
    dt_c = np.zeros(N)
    #dt_dyn[:] = cd*(2.*(r[1:] - (0.5*(r[1:] - r[0:-1])))/wrtt.g((r[1:]- (0.5*(r[1:]-r[0:-1]))),t,M_solm0))**0.5
    #dt_Cou[:] = c_C*((r[1:]-r[0:-1])/(cons.gamma*(cons.gamma - 1.)*u[0:-1])**0.5)
    #dt_v[:] = cv*np.abs((r[1:]-r[0:-1])/(v[1:]-v[0:-1]))
    dt_c[:] = c_c*np.abs(u[:]*rho[:]/((ne[:]*ni[:]*cl.lam_2(T[:]))))
    #dt_d = np.min(dt_dyn)
    #dt_C = np.min(dt_Cou)
    #dtv = np.min(dt_v)
    dtc = np.min(dt_c[Ncl:])
    return dtc
