# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:35:41 2017

@author: prakriti
"""
import sys
import numpy as np
import compton_ht as ch
from scipy.integrate import quad

sys.path.insert(0,'/home/prakriti/Documents/project_at_MPA/New_plan')
import constants as cons


tX = cons.yr*9.7e-8*(cons.M_BH/cons.solmass)*(ch.TX/1.e9)**(-1.5)

def integrand(t,t1,ta):
    return np.exp(-(t1-t)/ta)
    
def L_BH_t(L,ti,tnxt,es,r1,rho1,v1,Kappa):
    if v1<0:
        F_BH = -4.*np.pi*r1**2*(rho1*v1)
    else:
        F_BH = 0.0
    t_accr = Kappa*tX
    I = quad(integrand,ti,tnxt,(tnxt,t_accr))[0]
    #return np.exp(-tnxt/t_accr),F_BH,t_accr/cons.Myr
    return L*np.exp(-(tnxt-ti)/t_accr) + es*(cons.c**2)*F_BH*np.exp(-tnxt/t_accr)*I/t_accr