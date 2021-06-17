# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 11:19:55 2016

@author: prakriti
Given a z, finds t for lambda cdm with lambda_m = 0.25, lambda_cos = 0.75
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import optimize



def func(t,z0):
    dens_pm = 0.25
    dens_plamd = 0.75
    Km = 1.e5
    Mpc = 1.e6*3.086e18
    H0 = 73*Km/Mpc
    c1 = (dens_pm/dens_plamd)**(1./3.)
    c2 = 1.5*(dens_plamd**0.5)*H0
    return np.sinh(c2*t) - ((1.0/(c1*(1.0+z0)))**(3./2.))


if __name__ == '__main__':
    sol = optimize.root(func, 2.e17, 6.0)
    print (sol.x, sol.x/(3600*24*30*12*10**9))  #in Gyr
