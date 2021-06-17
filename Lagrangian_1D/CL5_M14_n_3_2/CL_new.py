# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 17:51:10 2017

@author: prakriti
"""
from __future__ import division
import numpy as np
import mypath as mypath
import sys
import os
import math
from scipy import optimize
import matplotlib.pyplot as plt

#PATH TO DIFERENT MODULES. NEEDS EDITING DEPENDING ON HOW IT IS RUN. THIS VERSION RUNS FROM TERMINAL.
sys.path.insert(0,mypath.path1)
import constants as cons
import cooling as cl
import WRT_t as wrtt
import Find_tz as ftz



sys.path.insert(0,mypath.path2)
import initial as ini
sys.path.insert(0,mypath.path2)
import make_initial as mi
sys.path.insert(0,mypath.path2)
import find_min_tme as fmd
sys.path.insert(0,mypath.path2)
import restart as rs
sys.path.insert(0,mypath.path2)
import parameters as para
sys.path.insert(0,mypath.path2)
import heating as ht
sys.path.insert(0,mypath.path2)
import subcycle as sc
sys.path.insert(0,mypath.path2)
import main_step as ms
#import compton_ht as ch
#import bremstra as brems
sys.path.insert(0,mypath.path2)
import interpol_be_se as ibs

M_solm0 = para.M_halo
t_f = optimize.root(ftz.func, 2.e17, 0.0).x
C = cons.kB*cons.fT/(cons.mu*cons.mp*(cons.mue*cons.mp)**(cons.gamma-1))
restart = False
cooling = True
I = 0.192
cf=.01
rd = 0.5*cons.Kpc
i_dmp = para.i_t_dmp


#Diemer-Kravtsov outer dnsity parameters
se = ibs.calc_interp()[0]
be = ibs.calc_interp()[1]
frac_rm = para.frac_rm


if (restart==True):
    L = 0.0
    t_i = rs.restart()[0]
    ri = rs.restart()[1]
    v_i = rs.restart()[2]
    rho_i = rs.restart()[3]
    prs_i = rs.restart()[4]
    u_i = rs.restart()[5]
    K_i = rs.restart()[6]
    T_i = rs.restart()[7]
    dmi = rs.restart()[8]
    Nri = rs.restart()[9]
    M_BH = rs.restart()[11]
    dts = rs.restart()[12]
    Ncool = int(rs.restart()[13])
    g1 = np.zeros(Nri)
    g1[:] = wrtt.g(ri[:],t_i,M_solm0)
    dt = fmd.find_min_dt(g1,ri,v_i,rho_i,T_i,u_i,t_i,dts,L,M_solm0,Ncool)   #time step, min dt
    dt_step = dt   #to create the effect dt = 0.5(dt_n-1/2   +   dt_n+1/2) and here we don't have dt_-1/2
    dt_nxt = dt_step





else:
    Ncool = 0
    t_i = ini.t0
    r200 = para.f_rad*wrtt.r200(t_i,M_solm0)
    rho_frac = para.f_mass
    Nri = para.N_r  #number of initial radial points
    r_in = para.r_in
    r_out = para.r_out
    ri = mi.ri_half(r_in,r_out,Nri)[0]   #array
    ri1by2 = mi.ri_half(r_in,r_out,Nri)[1]  #array

    dt_trial = 1.e13
    dmt = cons.fb*(M_solm0*cons.solmass - wrtt.M200(t_i,M_solm0))

    rho_i = mi.di_half(ri1by2,r200,rho_frac*wrtt.crit_rho(t_i),M_solm0)   #initial density array
    dmi = mi.dmi_half(rho_i, ri)   #dm array; this will remain a constant with time evolution

    prs_i = ini.K*rho_i**cons.gamma
    u_i = prs_i/((cons.gamma - 1.)*rho_i)
    K_i = prs_i/(C*rho_i**cons.gamma)
    T_i = prs_i*cons.mu*cons.mp/(cons.kB*rho_i)
    v_i = wrtt.H(t_i)*ri
    #v_i = np.zeros(Nri)
    #v_i[(ri>wrtt.r200(t_i,M_solm0))] = wrtt.H(t_i)*ri[(ri>wrtt.r200(t_i,M_solm0))]

    if (v_i[1]<0):
        L = -4.*np.pi*ri[0]**2*rho_i[0]*v_i[0]*para.eps*cons.c*cons.c
    else:
        L=0.0

    M_BH = para.M_BH
    g1 = np.zeros(Nri)
    g1[:] = wrtt.g(ri[:],t_i,M_solm0)
    dt = fmd.find_min_dt(g1,ri,v_i,rho_i,T_i,u_i,t_i,dt_trial,L,M_solm0,Ncool)   #time step, min dt
    dt_step = dt   #to create the effect dt = 0.5(dt_n-1/2   +   dt_n+1/2) and here we don't have dt_-1/2
    dt_nxt = dt_step



if (restart==True):
    rho_nxt_i = rho_i.copy()
    prs_nxt_i = prs_i.copy()
    u_nxt_i = u_i.copy()
    v_nxt_i = v_i.copy()
    r_nxt_i = ri.copy()
    T_nxt_i = T_i.copy()
    K_nxt_i = K_i.copy()
else:
    rho_nxt_i = np.zeros(Nri)
    prs_nxt_i = np.zeros(Nri)
    u_nxt_i = np.zeros(Nri)
    v_nxt_i = np.zeros(Nri)
    r_nxt_i = np.zeros(Nri)
    T_nxt_i = np.zeros(Nri)
    K_nxt_i = np.zeros(Nri)


q = np.zeros(Nri)
g_minus = np.zeros(Nri)
cq=4
m=0
r_c = ri[0]
t = t_i+dt_nxt
dmt = cons.fb*(wrtt.M200(t,M_solm0) - wrtt.M200(t_i,M_solm0))
dr200 = 2.*wrtt.r200(t, M_solm0) - ri[Nri-1]
if (restart==True):
    i_t = i_dmp*rs.restart()[10]
else:
    i_t=0
dt_init = dt_nxt


while (t<=t_f):
    if math.fmod(i_t,i_dmp)==0:
        data = np.zeros([Nri, 9])
        data[0,0] = t-dt_nxt
        data[1,0] = M_BH
        data[2,0] = dt_nxt
        data[3,0] = int(Ncool)
        data[:,1] = ri[:]
        data[:,2] = v_i[:]
        data[:,3] = rho_i[:]
        data[:,4] = prs_i[:]
        data[:,5] = u_i[:]
        data[:,6] = K_i[:]
        data[:,7] = T_i[:]
        data[:,8] = dmi[:]
        if not os.path.exists("./data"):
            os.makedirs("./data")
        np.savetxt("./data/t_%d"%(i_t)+".dat", data)
    i_t = i_t + 1
    print ("Timestep=%d"%(i_t)+", "+"Time (in sec) =%e"%(t))

    g_t = np.zeros(Nri)
    g_t[:] = wrtt.g(ri[:],t,M_solm0)

    rvdpu = ms.mainstp(g_t,dmi,ri,v_i,rho_i,prs_i,u_i,T_i,r_nxt_i,v_nxt_i,rho_nxt_i,prs_nxt_i,u_nxt_i,dt_nxt,dt_step,Ncool,cooling,r_c,1.e10,2.e4,1.05)
    r_nxt_i =  np.array(rvdpu[0])
    v_nxt_i = np.array(rvdpu[1])
    rho_nxt_i = np.array(rvdpu[2])
    prs_nxt_i = np.array(rvdpu[3])
    u_nxt_i = np.array(rvdpu[4])
    Ncool = rvdpu[5]

    if (cooling==True):
        u = u_nxt_i.copy()
        T = T_i.copy()
        uTp = sc.subcyc(r_nxt_i,v_nxt_i,u,rho_nxt_i,T,t-dt_nxt,t,dt_nxt,L,M_solm0,wrtt.r200(t,M_solm0), Ncool, M_BH, wrtt.M200(t,M_solm0))
        u_nxt_i[Ncool:] = np.array(uTp[0][Ncool:]).copy()
        #u_add = map(lambda u: u_flr if u<=0.0 else 0.0, u_nxt_i)
        #u_nxt_i[0:-1] = np.maximum(u_nxt_i[0:-1],u_add[0:-1])
        prs_nxt_i[Ncool:] = np.array(uTp[1][Ncool:]).copy()
        T_nxt_i[Ncool:] = np.array(uTp[2][Ncool:]).copy()

        #u_add = map(lambda u: u_flr if u<=0.0 else 0.0, u_nxt_i)
        #u_add = np.array(u_add)
        #u_nxt_i = np.maximum(u_nxt_i,u_add)
        #prs_nxt_i = np.maximum(prs_nxt_i, (cons.gamma - 1)*rho_nxt_i*u_add)
        #T_nxt_i = np.maximum(T_nxt_i, 2.5e4)
        K_nxt_i = prs_nxt_i/(C*rho_nxt_i**cons.gamma)

        #print uTp[4]
    else:
        K_nxt_i = prs_nxt_i/(C*rho_nxt_i**cons.gamma)
        T_nxt_i = prs_nxt_i*cons.mu*cons.mp/(cons.kB*rho_nxt_i)


    vdiff = v_nxt_i[1:]-v_nxt_i[0:-1]
    vdiff = np.append(vdiff,0.0)
    q = map(lambda vd,rhon,rhoi: -cq*(2./((1./rhon)+(1./rhoi)))*np.abs(vd)*(vd) if vd<0.0 else 0.0, vdiff,rho_nxt_i,rho_i)
    prs_nxt_i = prs_nxt_i + list(q)



    ri = r_nxt_i.copy()
    v_i = v_nxt_i.copy()
    rho_i = rho_nxt_i.copy()
    prs_i = prs_nxt_i.copy()
    u_i = u_nxt_i.copy()
    T_i = T_nxt_i.copy()
    K_i = K_nxt_i.copy()

    #L = ht.H(r_nxt_i,v_nxt_i,rho_nxt_i,T_i,t,dt_nxt,L)[1]
    dt_step = dt_nxt
    dt_full = fmd.find_min_dt(g_t,r_nxt_i,v_nxt_i,rho_nxt_i,T_nxt_i,u_nxt_i,t,dt_step,L,M_solm0,Ncool)
    dt_nxt = dt_full
    M_BH = sc.calc_BHgrth(t,dt_nxt,ri,rho_i,v_i,T_i,Ncool, M_BH)
    t = t+dt_nxt
