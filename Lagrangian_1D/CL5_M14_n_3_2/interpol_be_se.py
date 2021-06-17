import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def calc_interp():
    se = np.loadtxt("se.dat")
    be = np.loadtxt("be.dat")
    
    se_x = se[:,0]
    se_y = se[:,1]
    
    be_x = be[:,0]
    be_y = be[:,1]
    
    fse = interpolate.interp1d(se_x,se_y)
    fbe = interpolate.interp1d(be_x,be_y)
    
    return fse, fbe

#se = np.loadtxt("se.dat")
#be = np.loadtxt("be.dat")
#    
#se_x = se[:,0]
#se_y = se[:,1]
#    
#be_x = be[:,0]
#be_y = be[:,1]
#
#
#
#se0 = np.arange(0,7)
#se1 = calc_interp()[0](se0)
#
#be0 = np.arange(0,7)
#be1 = calc_interp()[1](be0)
#plt.plot(se_x, se_y, 'o', se0, se1, 'r-', label=r'$\mathbf{s_e}$')
#plt.plot(be_x,be_y,'o',be0,be1,'b-', label=r'$\mathbf{b_e}$')
#plt.legend()
#plt.xlabel(r'redshift', fontweight='bold')
#plt.show()
