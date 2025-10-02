import numpy as np
import cupy as cp
import massfunc as mf 
from astropy.constants import m_p
from astropy import units as u
from scipy.optimize import brentq
from scipy.integrate import quad_vec,quad


cosmo = mf.SFRD()
m_H = (m_p.to(u.M_sun)).value #M_sun
omega_b = cosmo.omegab
omega_m = cosmo.omegam

def xim(m,z):
    A,B,C,D,E,F,G = 4.4 , 0.334, 0.023, 0.199, -0.042, 0,1
    M7 = m/1e7
    F0 = 1.0
    return 1+A*M7**(B+C*np.log10(M7))*(F+G*(1+z)/10)* F0**(D+E*np.log10(F0))

def fstar(Mh): # Donnan 25 
    eps0 = 0.16
    Mc = 10**11.7
    beta = 0.9
    gamma = 0.65
    return 2*eps0 * ((Mh / Mc)**-beta + (Mh / Mc)**gamma)**-1


def qion_sb99(Muv):
    if Muv > 20:
        return 1000
    elif Muv >= 16.5:
        return 4065000 * np.exp(-0.420 * Muv)
    elif Muv >= 13:
        return -32336 * np.log10(Muv) + 43720
    else:
        return 7700
    
#### GPU version of interp1d
def interp1d_gpu(x, y, x_new):
    idx = cp.searchsorted(x, x_new, side='right') - 1
    idx = cp.clip(idx, 0, x.size - 2)
    dx = x[idx + 1] - x[idx]
    dy = y[idx + 1] - y[idx]
    return y[idx] + dy * (x_new - x[idx]) / dx