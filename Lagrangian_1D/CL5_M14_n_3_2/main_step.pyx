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


def mainstp(double[:] g, double[:] dm, double[:] r, double[:] v, double[:] d, double[:] p,
    double[:] u, double[:] T, double[:] r1, double[:] v1, double[:] d1, double[:] p1,
    double[:] u1, double dt_nxt, double dt_step, int Ncool, bint cooling,
    double r_c, double u_flr, double Tcool, double Rcool_by_r_c ):
    cdef int i
    cdef double dt_half, dm_h 
    cdef int N = r.shape[0]
    cdef double[:] u_add = np.zeros(N)
    cdef double gamma
    gamma = cons.gamma
    
    
    dt_half = 0.5*(dt_nxt + dt_step)
    
    for i in range(Ncool,N-1):
        dm_h = 0.5*(dm[1+i] + dm[i])
        v1[1+i] = v[1+i] - ((4.0*np.pi*(r[1+i]**2)*(p[1+i]-p[i])/dm_h) + g[1+i] )*dt_half
        r1[1+i] = r[1+i] +v1[1+i]*dt_nxt
    r1[Ncool] = r_c
    v1[Ncool] = 0.0

    if (cooling == True):
        if ((T[Ncool]<=Tcool)&(r1[1+Ncool]<= Rcool_by_r_c*r_c)):
            for i in range(1+Ncool):
                r1[i] = r_c
                v1[i] = 0.0 
            Ncool = Ncool+1
            r1[Ncool] = r_c
            v1[Ncool] = 0.0
      

    for i in range(Ncool,N-1):
        d1[i] = dm[i]/((4.*np.pi/3.0)*((r1[1+i]**3) - (r1[i]**3))) 
        u1[i] = u[i] - p[i]*((1./d1[i]) - (1./d[i]))
        if (u1[i] <= 0.0):
            u1[i] = u_flr
        p1[i] = (gamma - 1.0)*d1[i]*u1[i]
            
    delr = r1[N-1] - r1[N-2]
    d1[N-1] = dm[N-1]/((4.*np.pi/3.0)*(((r1[N-1]+delr)**3) - (r1[N-1]**3))) 
    if (u1[N-1]<=0.0):
        u1[N-1] = u_flr
    p1[N-1] = (gamma - 1.0)*d1[N-1]*u1[N-1]
    
    return [r1[i] for i in range(0,N) ], [v1[i] for i in range(0,N) ], [d1[i] for i in range(0,N) ],[p1[i] for i in range(0,N) ], [u1[i] for i in range(0,N) ], Ncool
