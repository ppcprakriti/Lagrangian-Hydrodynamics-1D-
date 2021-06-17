import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport abs
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
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

def lam(double T):
    cdef double a = 4.86567e-13
    cdef double b = -2.21974
    cdef double c = 1.35332e-5
    cdef double d = 9.64775
    cdef double e = 1.11401e-9
    cdef double f = -2.66528
    cdef double g = 6.91908e-21
    cdef double h = -0.571255
    cdef double i = 2.45596e-27
    cdef double j = 0.49521
    cdef double lamT    
  
    if ((T>5.e4)&(T<1.e10)):
        lamT = (a*T**b + (c*T)**d*(e*T**f + g*T**h))/(1 + (c*T)**d) + i*T**j
    else:
        lamT = 0.0 
    return lamT

def calc_BHgrth(double t_l, double dt,  double[:] r, double[:] rho, double[:] v, double[:] T, int Ncl, double M_BH):
    cdef int N = r.shape[0] 
    cdef double rhoc,vc,pc, Mth_ht, Medd_lim, M_acc
    cdef double pi,alpha,G,gamma,mp,rad_ef,Sigma_T,c,mu,kB, Kpc
    pi = np.pi
    alpha = 1.0
    G = cons.G
    gamma = cons.gamma
    mp = cons.mp
    rad_ef = 0.1
    Sigma_T = 6.65e-25
    c = cons.c
    mu = cons.mu
    kB = cons.kB
    Kpc = cons.Kpc
    rhoc = rho[Ncl]
    vc = v[Ncl+1]
    pc = (rhoc/(mu*mp))*kB*T[Ncl]
    Mth_ht = (4.*pi*alpha*G*G*M_BH*M_BH*rhoc)/((gamma*pc/rhoc) + vc*vc)**1.5
    Medd_lim = (4.*pi*G*M_BH*mp)/(rad_ef*Sigma_T*c)
    M_acc = min(Mth_ht,Medd_lim)
    M_BH = M_BH + M_acc*dt
    return M_BH

def Htng(double[:] r, double[:] v, double[:] rho, double[:] T, int Ncl, double M_BH, double M_cm):
    cdef int N = r.shape[0]
    cdef int i
    cdef double[:] dEdt = np.zeros(N)
    cdef double dr 
    cdef double rhoc,vc,pc, Mth_ht, Medd_lim, M_acc
    cdef double pi,alpha,G,gamma,mp,rad_ef,Sigma_T,c,mu,kB, r_ht, Kpc, eff
    pi = np.pi
    alpha = para.alpha
    G = cons.G
    gamma = cons.gamma
    mp = cons.mp
    rad_ef = para.rad_ef
    Sigma_T = ch.Sigma_T
    c = cons.c
    mu = cons.mu
    kB = cons.kB
    eff = para.eff
    r_ht = para.r_ht
    Kpc = cons.Kpc
    
    if (para.Heat == True):
     
        if (para.Thermht == True):
            rhoc = rho[Ncl]
            #if(v[1]<0):
            vc = v[Ncl+1]
            pc = (rhoc/(mu*mp))*kB*T[Ncl]
            Mth_ht = (4.*pi*alpha*G*G*M_BH*M_BH*rhoc)/((gamma*pc/rhoc) + vc*vc)**1.5
            Medd_lim = (4.*pi*G*M_BH*mp)/(rad_ef*Sigma_T*c)
            M_acc = min(Mth_ht,Medd_lim)
            for i in range(N-1):
                if (r[i+1]<=r_ht*Kpc):
                    dr = r[i+1] - r[i]
                    dEdt[i] = 3.*eff*M_acc*c*c/(4.*pi*(r_ht*Kpc)**3)
                   
               # else:
               #     dEdt[i] = 0.0
                
                
    #else:
    #    for i in range(N):
    #        dEdt[i] = 0.0
    
    return [dEdt[i] for i in range(0,N) ]
   



def find_min_dtc(double[:] r, 
    double[:] v, double[:] rho, 
    double[:] T, double[:] u, 
    double t, double dt, 
    double L, 
    double M_solm0, 
    int Ncl, double M_BH,double M_cm):
    cdef double c_c = 0.3
    cdef int N = r.shape[0]
    cdef double[:] n = np.zeros(N)
    cdef double[:] ne = np.zeros(N)
    cdef double[:] ni = np.zeros(N)
    cdef double[:] dt_c = np.zeros(N)
    cdef int il
    cdef double mp, mu, mue, temp
    mp = cons.mp
    mu = cons.mu
    mue = cons.mue

    for il in range(N):
        n[il] = rho[il]/(mp*mu)
        ne[il] = rho[il]/(mp*mue)
        ni[il] = n[il]-ne[il]
        temp = ((ne[il]*ni[il]*lam(T[il]))-(Htng(r,v,rho,T,Ncl,M_BH,M_cm)[il]) )
        if (temp== 0.0):
            dt_c[il] = 1.e40
        else:
            dt_c[il] = c_c*abs(rho[il]*u[il]/temp)
    dtc = min(dt_c[Ncl:N])
    #print dt/dtc
   
    return dtc
   
    

def subcyc(double[:] r, double[:] v, double[:] u, double[:] rho, double[:] T, double t_i, double t_f, double dt, double L, double M_solm0, double rvir, int Ncl, double M_BH, double M_cm):
    cdef double dtc, tt
    cdef int num, it
    cdef int N = r.shape[0]
    cdef double[:] n = np.zeros(N)
    cdef double[:] ne = np.zeros(N)
    cdef double[:] ni = np.zeros(N)
    cdef double[:] p = np.zeros(N)
    cdef int il
    cdef double mp, mu, mue, kB, gamma, temp
    mp = cons.mp
    mu = cons.mu
    mue = cons.mue
    kB = cons.kB
    gamma = cons.gamma

    for il in range(N):
        n[il] = rho[il]/(mp*mu)
        ne[il] = rho[il]/(mp*mue)
        ni[il] = n[il]-ne[il]
    dtc = find_min_dtc(r, v, rho, T, u, t_i, dt, L, M_solm0, Ncl, M_BH, M_cm)
    printf("minimum cooling time=%e\n",dtc)
    dtc = min(dtc, dt)
    printf("minimum of the hydro time and cooling time = %e\n", dtc)
    tt = t_i
    num = int(dt/dtc)
    printf("number of cool-heat cycles dt/min(dtc,dt) = %d\n", num)
    printf("current black hole mass = %e\n\n", M_BH)
    for it in range(num):
        tt = tt+dtc
        for il in range(N):
            #print Htng(r,v,rho,T,Ncl,M_BH)[il]
            temp = ((ne[il]*ni[il]*lam(T[il]))-(Htng(r,v,rho,T,Ncl,M_BH,M_cm)[il])) 
            if(r[il]<=rvir):
                #print Htng(r,v,rho,T,Ncl,M_BH)[il]
                if(temp>=0.0):
             
                    u[il] = u[il]/(1.+(((temp)*dtc)/(rho[il]*u[il])))
                else:
                    u[il] = u[il] - (temp*dtc/rho[il])     
                   
            
            p[il] = (gamma - 1.)*rho[il]*u[il]
            T[il] = p[il]*mu*mp/(kB*rho[il])
            
    
    return [u[il] for il in range(0,N) ], [p[il] for il in range(0,N) ], [T[il] for il in range(0,N) ], dtc
   
