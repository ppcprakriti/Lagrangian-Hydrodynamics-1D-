# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 18:01:53 2017

@author: prakriti
"""

from __future__ import division
import numpy as np
import mypath as mypath
import sys
import heating as ht
sys.path.insert(0,mypath.path1)
import constants as cons
import cooling as cl

def find_cltme(r,v,rho,T,u,t,dt,L):
    c_c = 0.8
    N = np.size(r)
    dt_c = np.zeros(N)
    dt_c[:] = c_c*np.abs(u[:]*rho[:]/(cl.lam(T[:]))) #ht.H(r,v,rho,T,t,dt,L)[0][:]-
    dtc = np.min(dt_c)
    return dtc