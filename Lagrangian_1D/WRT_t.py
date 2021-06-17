# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:06:20 2016

@author: prakriti

"""


from __future__ import division
import numpy as np
#import matplotlib.pyplot as plt
#import Find_tz as ftz
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import constants as cons
import solve_zf as s_zf
import Find_tz as ftz
from scipy import optimize

def a(t):
    c1 = (cons.omeg_m0/cons.omeg_lamb0)**(1./3.)
    c2 = 1.5*(cons.omeg_lamb0**0.5)*cons.H0
    return c1*(np.sinh(c2*t))**(2./3.)
    
def a_p(t):
    return a(t) - 1.0
        
def H(t):
    return cons.H0*(cons.omeg_lamb0 + cons.omeg_m0/(a(t)**3))**0.5

def z(t):
    return (1.-a(t))/a(t)
        
def M200(t, M_solm0):
    zfM0 = optimize.root(s_zf.func, .5, (M_solm0, cons.fr)).x
    nuM0 = s_zf.nu(zfM0,M_solm0)
    M200_0 = M_solm0*cons.solmass
    ztemp = z(t)
    if ztemp<=0.0:
        if np.abs(ztemp)<1.e-8:
            ztemp = 0.0
    return (M200_0)*(10.0**(-0.301*(np.log10(1+ztemp)/np.log10(1+zfM0))**nuM0))

def Mf(t, M_solm0):
    zfM0 = optimize.root(s_zf.func, .5, (M_solm0, cons.fr)).x
    nuM0 = s_zf.nu(zfM0,M_solm0)
    ztemp = z(t)
    if ztemp<=0.0:
        if np.abs(ztemp)<1.e-8:
            ztemp = 0.0
    return (10.0**(-0.301*(np.log10(1.+ ztemp)/np.log10(1.+zfM0))**nuM0)) - 0.04

def tf(M_solm0):
    return optimize.root(Mf, 1.e17 , M_solm0).x
    
def c(t, M_solm0):
    concent = 4.0*(1.0 + (t/(3.75*tf(M_solm0)))**8.4)**(1./8.)
    return concent
       
def F(c):
    return (np.log(1+c)-(c/(1+c)))
    
def crit_rho(t):
    return 3.0*H(t)**2/(8.0*np.pi*cons.G)
    
def r200(t, M_solm0):
    return ((3.0*M200(t,M_solm0)/(4.0*np.pi*200.0*crit_rho(t)))**(1./3.))

def r200m(t, M_solm0):
    return ((3.0*M200(t,M_solm0)/(4.0*np.pi*200.0*crit_rho(t)*cons.omeg_m0))**(1./3.))
    
def rs(t, M_solm0):
    return r200(t, M_solm0)/c(t, M_solm0)

def V200(t, M_solm0):
    return (cons.G*M200(t,M_solm0)/r200(t,M_solm0))**0.5

def g(r,t, M_solm0):
    b = M200(t, M_solm0)*cons.G/F(c(t,M_solm0))
    return b*(np.log(1 +r/rs(t,M_solm0))/(r*r) - (1/(r*(rs(t,M_solm0)+r)))) 

def phi(r,t,M_solm0):
    b = M200(t, M_solm0)*cons.G/F(c(t,M_solm0))
    return -b*np.log(1 +r/rs(t,M_solm0))/r

def rho_DM(r,t, M_solm0):
    del_c = (200./3.)*(c(t,M_solm0))**3/F(c(t,M_solm0))
    rsp = rs(t,M_solm0) 
    return del_c*crit_rho(t)/((r/rsp)*(1+(r/rsp))**2)   
def tff(r,t,M_solm0):
    return (2.0*r/g(r,t, M_solm0))**0.5

