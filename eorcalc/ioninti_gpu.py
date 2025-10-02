from .powerspec import MassFunctions
import numpy as np
import cupy as cp
import astropy.units as u
from cosfunc import n_H,dtdz
from .special import qion_sb99,interp1d_gpu,xim,fstar

class Ion:

    def __init__(self,z,fesc=0.2,A2byA1=0.1,ktrans=200,alpha=2.0,beta=0.0):
        self.cosmo = MassFunctions(A2byA1=A2byA1,kMpc_trans=ktrans,alpha=alpha,beta=beta)
        self.z = z
        self.mH = (self.cosmo.mHu.to(u.M_sun)).value #M_sun
        self.ob0 = self.cosmo.omegab
        self.om0 = self.cosmo.omegam
        self.nH = self.cosmo.nH  #cm^-3
        self.fesc = fesc
        self.deltaV_interp = np.linspace(-0.951, 2, 100)  # delta_V
        self.deltav_gpu = cp.asarray(self.deltaV_interp)
        self.M_min = self.cosmo.M_vir(0.61, 1e4, z)  # M_sun
        self.M_max = self.cosmo.M_vir(0.61, 2e8, z)
        self.M_J = self.cosmo.M_Jeans(z)  # M_sun
        self.qion = qion_sb99(z)  # ionizing photon per baryon in stars

    ### Source term
    def nion_interp(self,Mv:float,deltav:cp.ndarray)->cp.ndarray:
        def diff(m,Mv,deltaR):
            return fstar(m)*m*self.cosmo.dndmeps(m,Mv,deltaR,self.z)/ self.mH * self.ob0 / self.om0
        return self._nion_trapez(diff,Mv,deltav)*self.qion*self.fesc
    
    def nion_st(self,z):
        def diff(m):
            return fstar(m)*m*self.cosmo.dndmst(m,z)/ self.mH * self.ob0 / self.om0
        x = np.linspace(np.log10(self.M_min),np.log10(self.M_max),1000)
        y = diff(10**x)*10**x*np.log(10)
        return np.trapezoid(y,x)*self.qion*self.fesc
    
    def nion_ps(self,z):
        def diff(m):
            return fstar(m)*m*self.cosmo.dndmps(m,z)/ self.mH * self.ob0 / self.om0
        x = np.linspace(np.log10(self.M_min),np.log10(self.M_max),1000)
        y = diff(10**x)*10**x*np.log(10)
        return np.trapezoid(y,x)*self.qion*self.fesc
    
    ### Minihalo term
    def nxi_interp(self, Mv, deltav):
        def diff(m, Mv, deltav):
            return xim(m,self.z) * m * self.cosmo.dndmeps(m, Mv, deltav, self.z) / self.mH * self.ob0 / self.om0
        return self._nxi_trapez(diff, Mv, deltav)

    def nxi_st(self, z):
        def diff(m):
            return xim(m,z) * m * self.cosmo.dndmst(m,z) / self.mH * self.ob0 / self.om0
        x = np.linspace(np.log10(self.M_J),np.log10(self.M_min),1000)
        y = diff(10**x)*10**x*np.log(10)
        return np.trapezoid(y,x)

    ### IGM term
    def n_HI(self,deltaV):
        return ((1/u.cm**3).to(u.Mpc**-3)).value*n_H(deltaV) ### coming number density of hydrogen atoms in cm^-3
    
    def CHII(self,z):
        return 2.9*((1+z)/6)**-1.1

    def dnrec_dz_path(self,deltaV,xHII_Field:cp.ndarray)->cp.ndarray:
        x_HE = 1.08
        CIGM = self.CHII(self.z)
        # CIGM = 3.0
        nh = self.n_HI(deltaV)*(u.Mpc**-3).to(u.cm**-3)
        Q_HII = xHII_Field
        alpha_A = 4.2e-13 #cm**3/s
        # alpha_B = 2.59e-13  # cm^3/s at 10^4 K
        differential_trans = dtdz(self.z).to(u.s).value
        return CIGM*x_HE*alpha_A*nh*Q_HII*(1+self.z)**3 * differential_trans
    #### medium function




    #### method functions
    def _nxi_trapez(self, diff_func, Mv, deltav):
        x = np.linspace(np.log10(self.M_J), np.log10(self.M_min), 1000)
        deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, x, indexing='ij')
        m_vals = 10**m_grid
        y = diff_func(m_vals, Mv, deltaV_grid) * m_vals * np.log(10)
        integrand = np.trapezoid(y, x, axis=1)
        integrand = cp.asarray(integrand)
        return interp1d_gpu(self.deltav_gpu, integrand, deltav)
    
    def _nion_trapez(self, diff_func, Mv,deltav):
        Mh_max = 0.9*min(self.M_max,Mv)
        x = np.linspace(np.log10(self.M_min), np.log10(Mh_max), 1000)
        deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, x, indexing='ij')
        m_vals = 10**m_grid
        y = diff_func(m_vals, Mv, deltaV_grid) * m_vals * np.log(10)
        integrand = np.trapezoid(y, x, axis=1)
        integrand = cp.asarray(integrand)
        return interp1d_gpu(self.deltav_gpu, integrand, deltav)
    
    