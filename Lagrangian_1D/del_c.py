# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:06:20 2016

@author: prakriti

"""
from __future__ import division
import numpy as np

def del_c_woD(z,omeg_m0,omeg_lamb0):
    return 0.15*((12*np.pi)**(2./3.))*(((omeg_m0*(1+z)**3)/(omeg_lamb0 + omeg_m0*(1+z)**3))**0.0055)
    