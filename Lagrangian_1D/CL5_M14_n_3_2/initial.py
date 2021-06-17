# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 13:01:23 2016

@author: prakriti
"""

from __future__ import division
import sys
import mypath as mypath
import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import constants as cons
import Find_tz as ftz
from scipy import optimize
import WRT_t as wrtt
import cooling as cl
import matplotlib.pyplot as plt
sys.path.insert(0,mypath.path2)
import parameters as para
sys.path.insert(0,mypath.path1)
import WRT_t as wrtt

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

"""enter the initial redshift z0"""
z0 = 6.0
t0 = optimize.root(ftz.func, 2.e17, z0).x
"""enter initial pseudo-entropy"""
K_prime = para.KM13*(wrtt.M200(t0,para.M_halo)/wrtt.M200(t0,5.e13))**(2./3.)    #entropy in keVcm^2
#K_prime = para.KM13
K = K_prime*cons.kB*1.16e7/(cons.mu*cons.mp*(cons.mue*cons.mp)**(cons.gamma-1))

#"""enter current halo mass in units of solar mass"""
#M_solm0_1 = 5.0e13
#M_solm0_2 = 5.0e14


#t0 = 1.e17
def rho(r, r_out, rho_out,M_solm0):
    kons1 = (cons.gamma - 1.)/(cons.gamma*K)
    kons2 = (cons.gamma - 1.)
    return (rho_out**(kons2) + kons1*wrtt.phi(r_out,t0,M_solm0) - kons1*wrtt.phi(r,t0,M_solm0))**(1./kons2)
    
def prs(r, r_out, rho_out,M_solm0):
    return K*rho(r, r_out, rho_out,M_solm0)**cons.gamma
    
def Temp(r, r_out, rho_out,M_solm0):
    return prs(r, r_out, rho_out,M_solm0)*cons.mu*cons.mp/(cons.kB*rho(r, r_out, rho_out,M_solm0))

def u(r, r_out, rho_out,M_solm0):
    return prs(r, r_out, rho_out,M_solm0)/((cons.gamma - 1.)*rho(r, r_out, rho_out,M_solm0))

def dMdr(r,r_out,rho_out,M_solm0):
    return 4.0*np.pi*r*r*rho(r, r_out, rho_out,M_solm0)
    
def M(r_in,r_f, r_out,rho_out,M_solm0):
    return quad(dMdr,r_in,r_f,(r_out,rho_out,M_solm0))[0]
    
def tcool(r, r_out, rho_out,M_solm0):
    return 3.0*cons.mue*cons.mu*cons.mp*cons.kB*Temp(r, r_out, rho_out,M_solm0)/(2.0*rho(r, r_out, rho_out,M_solm0)*(cons.mue-cons.mu)*cl.lam(Temp(r, r_out, rho_out,M_solm0)))

  
#if __name__ == '__main__':    
#    r_g1 = np.linspace(0.1*cons.Kpc, wrtt.r200(t0,M_solm0_1), 200)
#    r_g2 = np.linspace(0.1*cons.Kpc, wrtt.r200(t0,M_solm0_2), 200)
#    tc1 = np.zeros(200)
#    tc2 = np.zeros(200)
#    for i in np.arange(200):
#        tc1[i] = tcool(r_g1[i], wrtt.r200(t0,M_solm0_1),28.*wrtt.crit_rho(t0),M_solm0_1)
#        tc2[i] = tcool(r_g2[i], wrtt.r200(t0,M_solm0_2),26.*wrtt.crit_rho(t0),M_solm0_2)
#    
#    plt.figure(figsize=[8,8])
#    plt.plot(r_g1[:]/wrtt.r200(t0,M_solm0_1), rho(r_g1[:], wrtt.r200(t0,M_solm0_1),28.*wrtt.crit_rho(t0),M_solm0_1)/(28.*wrtt.crit_rho(t0)),'b-',linewidth=3, label=r'$\mathbf{M_{halo0} = 5\times10^{13}M_{\odot}}$')
#    plt.plot(r_g1[:]/wrtt.r200(t0,M_solm0_1), Temp(r_g1[:], wrtt.r200(t0,M_solm0_1),28.*wrtt.crit_rho(t0),M_solm0_1)/1.16e7,'g-',linewidth=3) 
#    plt.plot(r_g1[:]/wrtt.r200(t0,M_solm0_1), tc1[:]/cons.Myr,'r-',linewidth=3) 
#    plt.plot(r_g1[:]/wrtt.r200(t0,M_solm0_1), wrtt.tff(r_g1[:],t0,M_solm0_1)/cons.Myr,'k-',linewidth=3)
#    plt.plot(r_g2[:]/wrtt.r200(t0,M_solm0_2), rho(r_g2[:], wrtt.r200(t0,M_solm0_2),26.*wrtt.crit_rho(t0),M_solm0_2)/(26.*wrtt.crit_rho(t0)),'b--',linewidth=3, label=r'$\mathbf{M_{halo0} = 5\times10^{14}M_{\odot}}$')
#    plt.plot(r_g2[:]/wrtt.r200(t0,M_solm0_2), Temp(r_g2[:], wrtt.r200(t0,M_solm0_2),26.*wrtt.crit_rho(t0),M_solm0_2)/1.16e7,'g--',linewidth=3) 
#    plt.plot(r_g2[:]/wrtt.r200(t0,M_solm0_2), tc2[:]/cons.Myr,'r--',linewidth=3)    
#    plt.plot(r_g2[:]/wrtt.r200(t0,M_solm0_2), wrtt.tff(r_g2[:],t0,M_solm0_2)/cons.Myr,'k--',linewidth=3)
#    plt.legend()
#    plt.xscale("log")
#    plt.yscale("log")
#    plt.xlim(0.1*cons.Kpc/wrtt.r200(t0,M_solm0_1), 1.0)
#    plt.xlabel(r'$\mathbf{r/r_{200}}$')
#    #plt.ylabel(r'$\boldsymbol{\rho/\rho_{200}}$')
#    plt.text(.06,2.1,r'$\boldsymbol{\rho/\rho_{200}}$')
#    plt.text(.054,0.25,r'$\mathbf{T(keV)}$')
#    plt.text(.6,6.6,r'$\mathbf{t_{cool}}$')
#    plt.text(.23,38.0,r'$\mathbf{t_{ff}}$')
#    #plt.title(r'$K_0 = 2 keVcm^2$')
##    
