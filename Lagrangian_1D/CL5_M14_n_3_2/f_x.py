# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 15:29:45 2017

@author: prakriti
"""
from __future__ import division
import numpy as np


def KN_cs(x):
    if x<1:
        return 1.0
    else:
        return 0.75*(((1.+x)/x**2)*((2*(1.+x)/(1.+2*x)) - (np.log(1.+2*x)/x)) + (np.log(1.+2*x)/(2*x)) - ((1.+3*x)/(1.+2*x)**2))
KN_cs_n = np.vectorize(KN_cs)
    
def E_param(x, t_rat):
    return ((1.+ (3*x**2/8))/(1.+ x**3)) - (x*(1.+x**2)/(4.*t_rat*(1.+ x**3)))
E_param_n = np.vectorize(E_param)