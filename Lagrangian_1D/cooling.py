# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 15:42:42 2017

@author: prakriti
"""

from __future__ import division
import numpy as np
import constants as cons
# import matplotlib
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
# import constants as cons
# import Find_tz as ftz
# from scipy import optimize
# import WRT_t as wrtt
# import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
# import initial as ini


def lam0(T):
    
    if hasattr(T, "__len__"):
        lam00 = np.zeros(len(T))
        lam00[np.where((T<=1.3e5) & (T>=1.e4))] = 5.3547e-27*T[np.where((T<=1.3e5) & (T>=1.e4))]
        lam00[np.where((T<=1.e8) & (T>1.3e5))] = (2.1786e-18*T[np.where((T<=1.e8) & (T>1.3e5))]**(-0.6826)) + (2.7060e-47*T[np.where((T<=1.e8) & (T>1.3e5))]**2.976)
        return lam00 
    else:
       if (T<=1.3e5 and T>=1.e4):
           return 5.3547e-27*T
       if (T<=1.e8 and T>1.3e5):
           return (2.1786e-18*T**(-0.6826)) + (2.7060e-47*T**2.976)
       else:
           return 0.0
       
def B0(T):
    T8 = 1.e8
    return lam0(T8)*(T/T8)**0.5
    
#def io_para(r)
    
def lam(T):
    lamT = np.zeros(len(T)) 
    lamT[np.where(T>=1.e8)] = B0(T[np.where(T>1.e8)])
    lamT[np.where(T<1.e8)] = lam0(T[np.where(T<1.e8)])
    return lamT
#    if T>1.e8:
#        return B0(T)
#    if T<=1.e8:
#        return lam0(T)
        
def lam_noarr(T):
    if T>1.e8:
        return B0(T)
    if T<=1.e8:
        return lam0(T)
        
#T = np.array([2*1e06, 2*1e07,2*1e08])
#
#print lam0(T)
#print B0(T)
#print lam(T)
        
def lam_1(T):
    T1 = T/cons.keV
    lamT = np.zeros(len(T1))
    lamT[np.where(T1>.02)]=1.e-22*(8.6e-3*(T1[np.where(T1>.02)])**(-1.7) + 0.058*(T1[np.where(T1>.02)])**(0.5) + 0.063)
    lamT[np.where((T1<=.02) & (T1>=.0017235))] = 6.72e-22*(T1[np.where((T1<=.02) & (T1>=.0017235))]/0.02)**0.6
    lamT[np.where((T1<.0017235))] = 1.544e-22*(T1[np.where((T1<.0017235))]/.0017235)**6
    return lamT
    
    
def lam_2(T):#zero metallicity; Wang et al. 2014
    a = 4.86567e-13
    b = -2.21974
    c = 1.35332e-5
    d = 9.64775
    e = 1.11401e-9
    f = -2.66528
    g = 6.91908e-21
    h = -0.571255
    i = 2.45596e-27
    j = 0.49521
    lamT = np.zeros(len(T))
    lamT[np.where((T<1.e10) & (T>2.e4))] = (a*T[np.where((T<1.e10) & (T>2.e4))]**b + (c*T[np.where((T<1.e10) & (T>2.e4))])**d*(e*T[np.where((T<1.e10) & (T>2.e4))]**f + g*T[np.where((T<1.e10) & (T>2.e4))]**h))/(1 + (c*T[np.where((T<1.e10) & (T>2.e4))])**d) + i*T[np.where((T<1.e10) & (T>2.e4))]**j
    #lamT = lamT + 1.e-45
    return lamT

def lam_sam(T):
    a = 4.86567e-13
    b = -2.21974
    c = 1.35332e-5
    d = 9.64775
    e = 1.11401e-9
    f = -2.66528
    g = 6.91908e-21
    h = -0.571255
    i = 2.45596e-27
    j = 0.49521

    if((T<1.e10) & (T>2.e4)):
        lamT=(a*T**b + (c*T)**d*(e*T**f + g*T**h))/(1 + (c*T)**d) + i*T**j
    else:
        lamT=0.0
    return lamT 
