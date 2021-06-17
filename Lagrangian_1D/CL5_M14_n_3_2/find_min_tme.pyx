import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport abs
from libc.stdlib cimport malloc, free
import sys
import mypath as mypath
import parameters as para
import compton_ht as ch
sys.path.insert(0,mypath.path1)
import constants as cons
#from cython.view cimport array as cvarray

#DTYPEi = np.int
#DTYPEf = np.float
#DTYPEd = np.double
#ctypedef np.int_t DTYPEi_t
#ctypedef np.float_t DTYPEf_t
#ctypedef np.double_t DTYPE_t
@cython.boundscheck(False)
@cython.cdivision(True)

def find_min_dt(double[:] g, double[:] r,double[:] v,double[:] rho,double[:] T,double[:] u,double t,double dt,double L,double M_solm0, int Ncl):
    cdef double mp,mu,mue,gamma
    cdef double cd, c_C, cv
    cdef int N = r.shape[0],i
    cdef double[:] n = np.zeros(N), ne = np.zeros(N), ni = np.zeros(N), dt_dyn = np.zeros(N-1), dt_Cou = np.zeros(N-1), dt_v = np.zeros(N-1)
    cdef double dt_d,dt_C,dtv
    mp = cons.mp
    mu = cons.mu
    mue = cons.mue
    gamma = cons.gamma
    cd = 0.1
    c_C = 0.3
    cv = 0.1
    for i in range(Ncl,N-1):
        dt_dyn[i] = cd*(2.*(r[i+1] - (0.5*(r[i+1] - r[i])))/g[i+1])**0.5
        dt_Cou[i] = c_C*((r[i+1]-r[i])/(gamma*(gamma - 1.)*u[i])**0.5)
        dt_v[i] = cv*np.abs((r[i+1]-r[i])/(v[i+1]-v[i]))

    dt_d = min(dt_dyn[Ncl:(N-1)])
    dt_C = min(dt_Cou[Ncl:(N-1)])
    dtv = min(dt_v[Ncl:(N-1)])
    
    return min(dt_d,dt_C,dtv)
