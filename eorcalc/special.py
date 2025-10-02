import numpy as np
import cupy as cp
import sys
import massfunc as mf 
from astropy.constants import m_p
from astropy import units as u



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

def load_binary_data(filename,DIM, dtype=np.float32) -> np.ndarray:
    f = open(filename, "rb")
    data = f.read()
    f.close()
    _data = np.frombuffer(data, dtype)
    if sys.byteorder == 'big':
        _data = _data.byteswap()
    data_cubic = _data.reshape((DIM, DIM, DIM), order='F')
    data_cubic = np.clip(data_cubic, -0.99, 2)
    return data_cubic

def TopHat_filter(data_cubic_ffted,R,DIM,box_length):
    ksmooth_func=get_ksmooth_func(DIM, box_length)
    kR = ksmooth_func*R
    kR[0,0,0]=1.0
    ksf=3.0*(cp.sin(kR)-kR*cp.cos(kR) )/kR**3
    ksf[0,0,0]=1.0
    kspace_data_cubic_smoothed2=data_cubic_ffted*ksf
    data_smoothed = cp.fft.irfftn(kspace_data_cubic_smoothed2,norm="forward")
    return data_smoothed

def get_ksmooth_func(DIM, box_length):
    kxx = cp.fft.fftfreq(DIM) * (2*cp.pi*DIM/box_length)
    kyy = cp.fft.fftfreq(DIM) * (2*cp.pi*DIM/box_length)
    kzz = cp.fft.rfftfreq(DIM) * (2*cp.pi*DIM/box_length)
    ks = cp.array(cp.meshgrid(kxx, kyy, kzz, indexing='ij'))
    ksmooth_func = cp.sqrt(ks[0]**2 + ks[1]**2 + ks[2]**2)
    ksmooth_func[0,0,0] = 1.0
    return ksmooth_func

def xHII_field_update(xHII_field,partial_eff):
    xHII_field[xHII_field < 1.0] = partial_eff[xHII_field < 1.0]
    xHII_field[xHII_field < 0.0] = 0.0
    xHII_field[xHII_field > 1.0] = 1.0
    xHII_field = cp.nan_to_num(xHII_field, nan=1.0)
    del partial_eff
    return xHII_field