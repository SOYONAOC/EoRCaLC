from .from_power import MassFunctions
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import interp1d
import astropy.units as u
from cosfunc import n_H,dtdz,H
import os
from scipy.integrate import quad_vec,quad
from astropy.constants import m_p



class Ion:

    def __init__(self,z,fesc=0.2,qion=4000,A2byA1=0.1,ktrans=200,alpha=2.0,beta=0.0):
        self.cosmo = MassFunctions(A2byA1=A2byA1,kMpc_trans=ktrans,alpha=alpha,beta=beta)
        self.A2byA1 = A2byA1
        self.kMpc_trans = ktrans
        self.alpha = alpha
        self.beta = beta
        self.z = z
        self.mH = (self.cosmo.mHu.to(u.M_sun)).value #M_sun
        self.ob0 = self.cosmo.omegab
        self.om0 = self.cosmo.omegam
        self.nH = self.cosmo.nH  #cm^-3
        self.fesc = fesc
        self.deltaV_interp = np.linspace(-0.999, 2, 1000)  # delta_V
        self.M_min = self.cosmo.M_vir(0.61, 1e4, z)  # M_sun
        self.M_max = self.cosmo.M_vir(0.61, 2e8, z)
        self.M_J = self.cosmo.M_Jeans(z)  # M_sun
        self.qion = qion

    ### Source term
    def nion_interp(self,Mv:float,deltaV:np.ndarray)->np.ndarray:
        def diff(m,Mv,deltaR):
            return self.fstar(m)*m*self.cosmo.dndmeps(m,Mv,deltaR,self.z)/ self.mH * self.ob0 / self.om0
        return self._nion_trapez(diff,Mv,deltaV)*self.qion*self.fesc
    
    def nion_st(self,z):
        def diff(m):
            return self.fstar(m)*m*self.cosmo.dndmst(m,z)/ self.mH * self.ob0 / self.om0
        x = np.linspace(np.log10(self.M_min),np.log10(self.M_max),1000)
        y = diff(10**x)*10**x*np.log(10)
        return np.trapezoid(y,x)*self.qion*self.fesc
    
    ### Minihalo term
    def eps_dz(self,M,Mv,deltaV,z):
        return (self.cosmo.dndmeps(M,Mv,deltaV,z+0.001*z) - self.cosmo.dndmeps(M,Mv,deltaV,z-0.001*z)) / (0.002*z)
    
    def nxi_dz_interp(self, deltaV, Mv):
        def diff(m, Mv, deltaV):
            return self.xim(m) * m * self.eps_dz(m, Mv, deltaV, self.z) / self.mH * self.ob0 / self.om0
        return self._nxi_trapez(diff, Mv, deltaV)
        
    def nxi_interp(self, deltaV, Mv):
        def diff(m, Mv, deltaV):
            return self.xim(m, self.z) * m * self.cosmo.dndmeps(m, Mv, deltaV, self.z) / self.mH * self.ob0 / self.om0
        return self._nxi_trapez(diff, Mv, deltaV)
    
    ### Minihaolo ST_over_PS
    def Nxi_ST(self,z):
        m_H = (m_p.to(u.M_sun)).value #M_sun
        Min = self.cosmo.M_Jeans(z)
        Max = self.cosmo.M_vir(0.61,1e4,z)
        def Nxi_ST_diff(m):
            return (self.xim(m,z)*1 / m_H * self.cosmo.omegab / self.cosmo.omegam * m * self.cosmo.dndmst(m, z))
        mslice = np.logspace(np.log10(Min), np.log10(Max), 100)
        ans = 0
        for i in range(len(mslice)-1):
            ans += quad(Nxi_ST_diff, mslice[i], mslice[i+1], epsrel=1e-5)[0]
        return ans
    
    def dST_dz(self,m,z):
        return (self.cosmo.dndmst(m,z+0.001*z) - self.cosmo.dndmst(m,z-0.001*z)) / (0.002*z)

    def dNxi_dz_ST(self,m,z):
        m_H = (m_p.to(u.M_sun)).value #M_sun
        omega_b = self.cosmo.omegab
        omega_m = self.cosmo.omegam
        return self.Xim(m,z)*m*self.dST_dz(m,z)/ m_H * omega_b / omega_m

    def dNxi_ST_Simpson(self,z):
        Min = self.cosmo.M_Jeans(z)
        Max = self.cosmo.M_vir(0.61,1e4,z)
        ms = np.linspace(np.log10(Min), np.log10(Max), 1000)
        dNxist = self.dNxi_dz_ST(10**ms,z)
        dNxist[dNxist > 0] = 0
        integrand = dNxist * 10**ms * np.log(10)
        arr = simpson(integrand, ms, axis=0)
        return arr

    def Nxi_ST_AtZ(self,z):
        dz = -0.1
        zlin = np.arange(20, z + dz/2, dz)
        dNxi_ST = self.dNxi_ST_Simpson(zlin)
        Nxi_st = abs(np.trapezoid(dNxi_ST, zlin)) + self.Nxi_ST(20)
        return Nxi_st


    ### IGM term
    def n_HI(self,deltaV):
        return n_H(deltaV).to(u.Mpc**-3).value  ### coming number density of hydrogen atoms in cm^-3
    
    def CHII(self,z):
        return 2.9*((1+z)/6)**-1.1

    def dtdz(self,z):
        return ((-1/(H(z)*(1+z))).to(u.s)).value

    def dnrec_dz_path(self,deltaV,xHII_Field:np.ndarray)->np.ndarray:
        x_HE = 1.08
        CIGM = self.CHII(self.z)
        # CIGM = 3.0
        nh = self.n_HI(deltaV)*(u.Mpc**-3).to(u.cm**-3)
        Q_HII = xHII_Field
        alpha_A = 4.2e-13 #cm**3/s
        # alpha_B = 2.59e-13  # cm^3/s at 10^4 K
        differential_trans = self.dtdz(self.z)
        return -CIGM*x_HE*alpha_A*nh*Q_HII*(1+self.z)**3 * differential_trans
    

    #### medium function
    def xim(self,m):
        A,B,C,D,E,F,G = 4.4 , 0.334, 0.023, 0.199, -0.042, 0,1
        M7 = m/1e7
        F0 = 1.0
        return 1+A*M7**(B+C*np.log10(M7))*(F+G*(1+self.z)/10)* F0**(D+E*np.log10(F0))
    
    def fstar(self, Mh): # Donnan 25 
        eps0 = 0.16
        Mc = 10**11.7
        beta = 0.9
        gamma = 0.65
        return 2*eps0 * ((Mh / Mc)**-beta + (Mh / Mc)**gamma)**-1


    #### method functions
    def _nxi_trapez(self, diff_func, Mv, deltaV):
        x = np.linspace(np.log10(self.M_J), np.log10(self.M_min), 1000)
        deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, x, indexing='ij')
        m_vals = 10**m_grid
        y = diff_func(m_vals, Mv, deltaV_grid) * m_vals * np.log(10)
        integrand = np.trapezoid(y, x, axis=1)
        interp_func = interp1d(self.deltaV_interp, integrand, kind='cubic', bounds_error=False, fill_value=0)
        return interp_func(deltaV)
    
    def _nion_trapez(self, diff_func, Mv, deltaV):
        Mh_max = 0.9*min(self.M_max,Mv)
        x = np.linspace(np.log10(self.M_min), np.log10(self.M_max), 1000)
        deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, x, indexing='ij')
        m_vals = 10**m_grid
        y = diff_func(m_vals, Mv, deltaV_grid) * m_vals * np.log(10)
        integrand = np.trapezoid(y, x, axis=1)
        interp_func = interp1d(self.deltaV_interp, integrand, kind='cubic', bounds_error=False, fill_value=0)
        return interp_func(deltaV)