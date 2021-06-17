from __future__ import division
import numpy as np
import subcycle as sbc
import matplotlib.pyplot as plt



r = np.zeros(2)
u = np.zeros(2)
rho = np.zeros(2)
T = np.zeros(2)


r[0] = 3.086e21
r[1] = 5.e21
u[0] = 4.e14
u[1] = 5.e13
rho[0] = 4.e-24
rho[1] = 3.e-24
T[0] = 5.e8
T[1] = 3.e7

p = sbc.subcyc(r,u,rho,T,4.e16, 4.5e16, .5e16, 0.0, 5.e14, 1.5e22, 0)
print p
#dt = sbc.find_min_dtc(r, rho, T, u, 4.e16, 0.5e16, 0.0, 5.e14, 0)
#print dt 

