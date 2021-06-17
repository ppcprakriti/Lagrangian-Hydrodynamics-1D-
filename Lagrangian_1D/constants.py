# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:05:39 2016

@author: prakriti
"""
from __future__ import division

global h
global c
global H0
global sigma8
global omeg_m0
global omeg_lamb0
global Kpc
global km
global Mpc
global G
global solmass
global fr    # used in Van den Bosch prescription fr*M0
global Gyr
global Myr
global mu
global mue
global mp
global me
global gamma
global epsilon #accretion efficiency
global fb #universal baryon fraction
global yr
global fT
global M_BH
global h_pl
global Z
global qe
global keV

h = 0.73
c = 3.e10
km = 1.e5
Mpc = 1.e6*3.086e18
H0 = 73*km/Mpc
G = 6.67e-8
omeg_m0 = 0.25
omeg_lamb0 = 0.75
Kpc = 3.086e21
sigma8 = 0.9
solmass = 2.e33
fr = 0.254
Gyr = 3600*24*30*12*10**9
Myr = 3600*24*30*12*10**6
yr = 3600*24*30*12
mu = 0.62
mue = 1.17
mp = 1.6726e-24
me = 9.1e-28
kB = 1.38e-16
gamma = 5./3.
epsilon = .001
fb = 0.17
fT = 1.16e7
M_BH = 1.e8*solmass
h_pl = 6.626e-27
Z = 1
qe = 4.803e-10
keV = 1.16e7