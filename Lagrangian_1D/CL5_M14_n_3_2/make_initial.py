# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 18:21:27 2017

@author: prakriti
"""

from __future__ import division
#from numba import jit
import sys
import mypath as mypath
import numpy as np
import matplotlib.pyplot as plt
import constants as cons
import initial as ini
import interpol_be_se as ibs
import parameters as para
sys.path.insert(0,mypath.path1)
import WRT_t as wrtt

se = ibs.calc_interp()[0]
be = ibs.calc_interp()[1]


def ri_half1(r_frac,r_in,N):
    r = np.zeros(N)
    dr = np.zeros(N-1)
    r_half = np.zeros(N)
    r = r_frac*(r_in)**(np.arange(N))
    dr[:] = r[1:]-r[0:-1]
    r_half[0:-1] = r[0:-1] + 0.5*dr[0:]
    r_half[N-1] = r[N-1] + 0.5*(r_frac*(r_in)**(N) - r[N-1] )
    return r,r_half
    
def ri_half(r_in,r_out,N):
    rgi = np.linspace(r_in,r_out,N)
    dr = (r_out - r_in)/N
    rgi_h = np.linspace(r_in+(0.5*dr),r_out+(0.5*dr),N)
    return rgi,rgi_h,dr

#@jit
def di_half(r,r_out,rho_out,M_solm0):
    N = np.size(r)
    d=np.zeros(N)
    d_outer = np.zeros(N)
    ftrans = np.zeros(N)
    #d_outer[:] = (cons.fb*wrtt.crit_rho(ini.t0)*cons.omeg_m0)*((be(wrtt.z(ini.t0))*(r[:]/(para.frac_rm*wrtt.r200m(ini.t0,M_solm0)))**(-se(wrtt.z(ini.t0))))  + 1.)
    #ftrans[:] = (1. + 0.2*(r[:]/wrtt.r200m(ini.t0,M_solm0))**4)**(-1.5)
    #d[:] = ini.rho(r[:],r_out,rho_out,M_solm0)
    #d_outer[r > wrtt.r200(ini.t0,M_solm0)] = (cons.fb*wrtt.crit_rho(ini.t0)*cons.omeg_m0)*((be(wrtt.z(ini.t0))*(r[r > wrtt.r200(ini.t0,M_solm0)]/(para.frac_rm*wrtt.r200m(ini.t0,M_solm0)))**(-se(wrtt.z(ini.t0))))  + 1.)
    #d[r <= wrtt.r200(ini.t0,M_solm0)] = ini.rho(r[r <= wrtt.r200(ini.t0,M_solm0)],r_out,rho_out,M_solm0)
    #d[(r > wrtt.r200(ini.t0,M_solm0))] = d_outer[(r > wrtt.r200(ini.t0,M_solm0))]
    #d[:] = ini.rho(r[:],r_out,rho_out,M_solm0)*ftrans[:] + d_outer[:]
#    d[r<=9.*wrtt.r200(ini.t0,M_solm0)] = ini.rho(r[r<=9.*wrtt.r200(ini.t0,M_solm0)],r_out,rho_out,M_solm0)*ftrans[r<=9.*wrtt.r200(ini.t0,M_solm0)] + d_outer[r<=9.*wrtt.r200(ini.t0,M_solm0)]
#    d[r>9.*wrtt.r200(ini.t0,M_solm0)]=cons.fb*wrtt.crit_rho(ini.t0)
    a1 = para.a_mid
    indx = para.indx
    indx1 = para.indx1
    
    d[r<=wrtt.r200(ini.t0,M_solm0)]= ini.rho(r[r<=wrtt.r200(ini.t0,M_solm0)],r_out,rho_out,M_solm0)
    dr200 = d[r<=wrtt.r200(ini.t0,M_solm0)]
    l1 = len(dr200)
    d[(r>wrtt.r200(ini.t0,M_solm0))&(r<=a1*wrtt.r200(ini.t0,M_solm0))] =((dr200[l1-1])*(r[(r>wrtt.r200(ini.t0,M_solm0))&(r<=a1*wrtt.r200(ini.t0,M_solm0))]/wrtt.r200(ini.t0,M_solm0))**(-indx))
  #  d1_5r200 = d[(r>wrtt.r200(ini.t0,M_solm0))&(r<=a1*wrtt.r200(ini.t0,M_solm0))]
  #  l2 = len(d1_5r200)
  #  d[(r>a1*wrtt.r200(ini.t0,M_solm0))&(r<=a2*wrtt.r200(ini.t0,M_solm0))] = ((d1_5r200[l2-1])*(r[(r>a1*wrtt.r200(ini.t0,M_solm0))&(r<=a2*wrtt.r200(ini.t0,M_solm0))]/wrtt.r200(ini.t0,M_solm0))**(-0.5))
    d2r200 = d[(r>wrtt.r200(ini.t0,M_solm0))&(r<=a1*wrtt.r200(ini.t0,M_solm0))]
    l = len(d2r200)
    d[(r>a1*wrtt.r200(ini.t0,M_solm0))]= d2r200[l-1]*(r[(r>a1*wrtt.r200(ini.t0,M_solm0))]/(2.*wrtt.r200(ini.t0,M_solm0)))**(-indx1)
    
    return d

#@jit
def dmi_half(di_half,r):
    N = len(r)
    dmi_h = np.zeros(N)
#    for i in np.arange(N-1):
    dmi_h[:-1] = di_half[0:-1]*((4.0/3.0)*np.pi*(r[1:]**3 - r[0:-1]**3))
    dr = r[N-1] - r[N-2]
    dmi_h[N-1] = di_half[N-1]*((4.0/3.0)*np.pi*((r[N-1]+dr)**3 - r[N-1]**3))
    
    return dmi_h

#@jit        
def pi_half(r,r_out,rho_out,M_solm0):
    N = len(r)
    p=np.zeros(N)
#    for i in np.arange(N):
    p[:] = ini.prs(r[:],r_out,rho_out,M_solm0)
    return p

#@jit    
def ui_half(r,r_out,rho_out,M_solm0):
    N = len(r)
    u=np.zeros(N)
#    for i in np.arange(N):
    u[:] = ini.u(r[:],r_out,rho_out,M_solm0)
    return u
        
def rho_dk(r,t,M_solm0):
    return (cons.fb*wrtt.crit_rho(t)*cons.omeg_m0)*((be(wrtt.z(t))*(r/(para.frac_rm*wrtt.r200m(t,M_solm0)))**(-se(wrtt.z(t))))  + 1.  )
