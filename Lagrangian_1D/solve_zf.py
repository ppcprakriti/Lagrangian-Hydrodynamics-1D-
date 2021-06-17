# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:06:20 2016

@author: prakriti

"""


from __future__ import division
import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import optimize
import del_c as dc
import sigma2 as s2
import constants as cons

def H_by_H0(z):
    return (cons.omeg_lamb0 + cons.omeg_m0*(1+z)**3)**(1./2.)

def H_by_H01(z):
    return (0.75 + 0.25*(1+z)**3)**(1./2.)

def integrand(z):
    return (1+z)/(cons.omeg_lamb0 + cons.omeg_m0*(1+z)**3)**(3./2.)

def integrand1(z):
    return (1+z)/(0.75+ 0.25*(1+z)**3)**(3./2.)

def integrate(z_l):
    return quad(integrand, z_l, np.inf)[0]

def integrate1(z_l):
    return quad(integrand1, z_l, np.inf)[0]

def func(z_l, M_solm, fr):
    return -(0.477*((2.*(s2.sigma2(fr*M_solm, cons.omeg_m0) - s2.sigma2(M_solm,cons.omeg_m0)))**0.5)) - (dc.del_c_woD(0.0, cons.omeg_m0,cons.omeg_lamb0)) + (dc.del_c_woD(z_l,cons.omeg_m0,cons.omeg_lamb0)/(H_by_H0(z_l)*integrate(z_l)/integrate(0)))

def func1(z_l):
    return -1.262 +  ((((0.25*(1+z_l)**3)/(0.75 + 0.25*(1+z_l)**3))**0.0055)/(H_by_H01(z_l)*integrate1(z_l)/integrate1(0.0)))

def nu(zf,M_solm):
    return 1.211 + 1.858*np.log10(1+zf) + 0.308*cons.omeg_lamb0**2 - 0.032*np.log10(M_solm*s2.h/1.e11)

if __name__ == '__main__':
    M_solm = 5.24e14
    omeg_m0 = 0.25
    omeg_lamb0 = 0.75
    f = 0.254
    sol = optimize.root(func, .5, (M_solm, omeg_m0, omeg_lamb0, f))
    print (sol.x, nu(sol.x,M_solm,omeg_lamb0))

#z_f = np.linspace(0,0.6, 400)
#f = np.zeros(400)
#for i in np.arange(400):
#    f[i] = func(z_f[i])
#
#
#plt.figure()
#plt.plot(z_f, f)
