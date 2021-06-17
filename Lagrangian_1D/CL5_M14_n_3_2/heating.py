# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 17:31:35 2017

@author: prakriti
"""

import sys
import mypath as mypath
import numpy as np
#from numba import jit
sys.path.insert(0,mypath.path2)
import compton_ht as ch
#import L_BH_comp as Lbh
#import bremstra as brems
import parameters as para
sys.path.insert(0,mypath.path1)
import constants as cons

#@jit
def H(r,v,d,T,t,dt,Lbhp):
    ti = t-dt
    tf=t
    #L_BH0 = Lbh.L_BH_t(Lbhp,ti,tf,para.eps,r[1],d[1],v[1],para.Kappa)
    L_BH0 = Lbhp
    dEdt = np.zeros(np.size(r))
    if (para.Heat==True):
        if (para.Thermht == True):
            p = (d/(cons.mu*cons.mp))*cons.kB*T
            rhoc = d[0]
            vc = v[1]
            pc = p[0]
            Mth_ht = (4.*np.pi*para.alpha*cons.G*cons.G*para.M_BH*para.M_BH*rhoc)/((cons.gamma*pc/rhoc) + vc*vc)**1.5
            Medd_lim = (4.*np.pi*cons.G*para.M_BH*cons.mp)/(para.rad_ef*ch.Sigma_T*cons.c)
            M_acc = min(Mth_ht,Medd_lim)
#            indx = np.where(r<=para.r_ht*cons.Kpc)
#            for i in np.arange(len(indx)):
#                if indx[i]!=0:
#                    dEdt[indx[i]-1] = para.eff*M_acc*cons.c*cons.c
            delr = r[1:]-r[0:-1]
            dEdt[np.where(r[1:]<=para.r_ht*cons.Kpc)] = 3.*para.eff*M_acc*cons.c*cons.c*r[1:][np.where(r[1:]<=para.r_ht*cons.Kpc)]*r[1:][np.where(r[1:]<=para.r_ht*cons.Kpc)]*(delr[np.where(r[1:]<=para.r_ht*cons.Kpc)])/((para.r_ht*cons.Kpc)**3)
       # if (para.Compht==True):
       #     dEdt[:] = ch.dE_comp(r[:],d[:],T[:],L_BH0) + brems.dE_br(r[:],d[:],T[:]) 
    else:
        dEdt[:] = 0.0
            
    return dEdt, L_BH0


    
    
    
    
    
    
