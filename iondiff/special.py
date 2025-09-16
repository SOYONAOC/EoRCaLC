import numpy as np
import sys
from functools import lru_cache
import massfunc as mf 
from astropy.constants import m_p
from astropy import units as u
from scipy.optimize import brentq
from scipy.integrate import quad_vec,quad

cosmo = mf.SFRD()
m_H = (m_p.to(u.M_sun)).value #M_sun
omega_b = cosmo.omegab
omega_m = cosmo.omegam

def dST_dz(m,z):
    return (cosmo.dndmst(m,z+0.001*z) - cosmo.dndmst(m,z-0.001*z)) / (0.002*z)

def Xim(m,z):
    A,B,C,D,E,F,G = 4.4 , 0.334, 0.023, 0.199, -0.042, 0,1
    M7 = m/1e7
    F0 = 1.0
    return 1+A*M7**(B+C*np.log10(M7))*(F+G*(1+z)/10)* F0**(D+E*np.log10(F0))

def dNxi_dz_ST(m,z):
    return Xim(m,z)*m*dST_dz(m,z)/ m_H * omega_b / omega_m

def Nxi_dz_ST(z):
    Mj = cosmo.M_J(z)
    Mmax = cosmo.M_vir(0.61,1e4,z)
    if dST_dz(Mj,z)>0 and dST_dz(Mmax,z)>0:
        return 0
    if dST_dz(Mj,z)<0 and dST_dz(Mmax,z)<0:
        M = np.logspace(np.log10(Mj),np.log10(Mmax),12)
        ans = 0
        for i in range(len(M)-1):
            ans += quad(dNxi_dz_ST,M[i],M[i+1],args=(z))[0]
        return ans
    if dST_dz(Mj,z)>0 and dST_dz(Mmax,z)<0:
        Mmin = brentq(dST_dz, Mj, Mmax, args=(z))
        M = np.logspace(np.log10(Mmin),np.log10(Mmax),12)
        ans = 0
        for i in range(len(M)-1):
            ans += quad_vec(dNxi_dz_ST,M[i],M[i+1],args=(z))[0]
        return ans

def Nxi_ST_Interp():
    zstep = np.linspace(4,30,1000)
    return 