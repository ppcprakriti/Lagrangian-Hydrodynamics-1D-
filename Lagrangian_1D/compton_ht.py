# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 12:48:41 2016

@author: prakriti
"""
from __future__ import division
import numpy as np
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import constants as cons
import Find_tz as ftz
from scipy import optimize
import WRT_t as wrtt
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import initial as ini

CX = 0.137
TX = 1.04e9
Sigma_T = 6.65e-25 
nu_T = cons.me*cons.c*cons.c/cons.h

def L_BH(r, d, T, L_BH0):
    E = 1.5*cons.kB*T
    ne = d/(cons.mp*cons.mue)
    return L_BH0*np.exp(8.0*E*ne*Sigma_T*CX*(1. - (TX/T))/(3.0*cons.me*cons.c*cons.c*nu_T))
    
def dE_comp(r, d, T, L_BH0):
    E = 1.5*cons.kB*T
    ne = d/(cons.mp*cons.mue)
    return -8.0*E*ne*L_BH(r, d, T, L_BH0)*Sigma_T*CX*(1. - (TX/T))/(3.0*cons.me*cons.c*cons.c*nu_T*4.0*np.pi*r*r)
    
