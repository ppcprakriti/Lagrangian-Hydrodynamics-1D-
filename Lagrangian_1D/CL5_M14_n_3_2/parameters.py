# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 16:21:05 2017

@author: prakriti
"""
import sys
import mypath as mypath
sys.path.insert(0,mypath.path1)
import constants as cons
i_t_dmp = 100 

#heating
Heat = True
Compht = False
Thermht = True
eps=.01
Kappa = 5.e7

#Halo details
N_r = 180
M_halo = 5.e14   #in solar mass
f_mass = 20.0    #f_mass*crit_rho(t) is density at outer radius initially
f_rad = 2.0
KM13 = 2.0 
frac_rm = 5.0
r_in = 1.*cons.Kpc
r_out = 901.*cons.Kpc

#initial density profile in the outskirts

a_mid = 2.0 # times r200
indx = 4.0
indx1 = 0.0
#Blackhole

M_BH = 1.e8*cons.solmass
alpha = 1.0
rad_ef = 0.1
eff = .5
r_ht = 27.  #kpc

