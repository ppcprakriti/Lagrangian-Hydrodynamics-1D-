# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:06:20 2016

@author: prakriti

"""
from __future__ import division
import numpy as np

global h
global sigma8

sigma8 = 0.9
h = 0.73

def Tao(omeg_m0):
    omeg_b = 0.019
    return omeg_m0*h*np.exp(-omeg_b*(h**(-2))*(1. + (((2.*h)**0.5)/omeg_m0)))
    
def u(M_solm, omeg_m0):
    return 3.804e-4*Tao(omeg_m0)*((M_solm*h/omeg_m0)**(1./3.))
    
def u8(omeg_m0):
    return 32.*Tao(omeg_m0)
    
def f8(omeg_m0):
    return (64.087)*((1. + 1.074*(u8(omeg_m0))**0.3 - 1.581*(u8(omeg_m0))**0.4 + 0.954*(u8(omeg_m0))**0.5 - 0.185*(u8(omeg_m0))**0.6)**(-10.))
    
def f(M_solm, omeg_m0):
    return (64.087)*((1. + 1.074*(u(M_solm,omeg_m0))**0.3 - 1.581*(u(M_solm,omeg_m0))**0.4 + 0.954*(u(M_solm,omeg_m0))**0.5 - 0.185*(u(M_solm,omeg_m0))**0.6)**(-10.))
 
def sigma2(M_solm,omeg_m0):
    return (sigma8*f(M_solm,omeg_m0)/f8(omeg_m0))**2
    

   