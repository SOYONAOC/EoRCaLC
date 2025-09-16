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

    def __init__(self,z,fesc=0.2,kakaka=1.4*1e-28,xi_ion=10**25.6,A2byA1=0.1,ktrans=200,alpha=2.0,beta=0.0):
        self.cosmo = MassFunctions(A2byA1=A2byA1,kMpc_trans=ktrans,alpha=alpha,beta=beta)
        self.A2byA1 = A2byA1
        self.kMpc_trans = ktrans
        self.alpha = alpha
        self.beta = beta
        self.z = z
        self.kakaka = kakaka
        self.xi_ion = xi_ion
        self.fesc = fesc
        self.deltaV_interp = np.linspace(-0.999, 2, 1000)  # delta_V
        self.M_min = self.cosmo.M_vir(0.61, 1e4, z)  # M_sun
        self.M_J = self.cosmo.M_Jeans(z)  # M_sun


    ### Source term
    def fstar(self, M):
        f0 = .14
        ylo = .46
        yhi = .82
        Mp = 10**12.3  # M_sun solmass
        fup = 2 * f0
        fdown = ((M / Mp)**-ylo + (M / Mp)**yhi)
        return fup / fdown

    def fduty(self, M):
        al = 1.5
        Mc = 6e7
        return (1 + (2.**(al / 3.) - 1) * (M / Mc)**-al)**(-3. / al)

    def dMdt(self, M, z):
        return 24.1 * (M / (1e12))**1.094 * (1 + 1.75 * z) * np.sqrt(self.cosmo.omegam * (1 + z)**3 + (1-self.cosmo.omegam))  # solmass/yr

    def rhosfrdiff(self, Mh, Mv,deltaV,z):
        diff = self.fstar(Mh) * self.cosmo.omegab / self.cosmo.omegam * self.dMdt(Mh, z) * self.cosmo.dndmeps(Mh,Mv,deltaV, z) * self.fduty(Mh)
        return diff
    
    def rhosfr(self,Mv,deltaV):
        try:
            Nxi_arr = np.load(f'.rhosfr_dz_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{Mv:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy')
        except FileNotFoundError:
            os.makedirs(f'.rhosfr_dz_Interp_init/NxiAtz{self.z:.2f}', exist_ok=True)
            M_min = self.cosmo.M_vir(0.61, 1e4, self.z)
            M_max = min(self.cosmo.M_vir(0.61, 1e8, self.z),0.9*Mv)
            Mh_interp = np.linspace(np.log10(M_min), np.log10(M_max), 1000)  # log10(M)
            deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, Mh_interp, indexing='ij')
            mh = self.rhosfrdiff(Mh=10**m_grid, Mv=Mv,deltaV=deltaV_grid, z=self.z)
            integrand = mh*10**m_grid * np.log(10)  # dN/dlogM -> dN/dM
            # Integrate using Simpson's rule
            Nxi_arr_save = simpson(integrand, x=Mh_interp, axis=1)
            np.save(f'.rhosfr_dz_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{Mv:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy', Nxi_arr_save)
            Nxi_arr = np.load(f'.rhosfr_dz_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{Mv:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy')
        Simpson_Interp = interp1d(self.deltaV_interp, Nxi_arr, kind='cubic')
        return Simpson_Interp(deltaV)

    def rho_sfr_st(self):
        M_min = self.cosmo.M_vir(0.61, 1e4, self.z)
        M_max = self.cosmo.M_vir(0.61, 1e8, self.z)
        Mh_interp = np.linspace(np.log10(M_min), np.log10(M_max), 1000) 
        def rhodiff(Mh):
            return self.fstar(Mh) * self.cosmo.omegab / self.cosmo.omegam * self.dMdt(Mh, self.z) * self.cosmo.dndmst(Mh,self.z) * self.fduty(Mh)
        mh = rhodiff(Mh=10**Mh_interp)
        integrand = mh*10**Mh_interp * np.log(10)  # dN/dlogM -> dN/dM
        return simpson(integrand, Mh_interp, axis=0)
    
    def rho_UV(self,Mv,deltaV):
        return self.rhosfr(Mv,deltaV=deltaV)/self.kakaka
    
    def nion_dot_st(self):
        return self.fesc * self.xi_ion * self.rho_sfr_st()/self.kakaka * u.s**-1 *u.Mpc**-3
    
    def nion_dot(self,Mv,deltaV):
        return self.fesc * self.xi_ion * self.rho_UV(Mv,deltaV) * u.s**-1 *u.Mpc**-3

    def nion_dz(self,Mv,deltaV):
        return (self.nion_dot(Mv,deltaV)/n_H(deltaV) * dtdz(self.z)).to(u.dimensionless_unscaled).value * (n_H(0).to(u.Mpc**-3)).value  # comoving

    def nion_dz_st(self):
        return (self.nion_dot_st()/n_H(0) * dtdz(self.z)).to(u.dimensionless_unscaled).value * (n_H(0).to(u.Mpc**-3)).value  # comoving
    
    ### IGM term
    def trec(self,deltaV):
        alpha_B = 2.59e-13  # cm^3/s at 10^4 K
        clumping_factor = 3.0
        Xe = 1.09
        return 1.0 / (clumping_factor * alpha_B * (1 + deltaV) * Xe*n_H(deltaV).to(1/u.cm**3).value * (1 + self.z)**3)  # s

    def n_HI(self,deltaV):
        return n_H(deltaV).to(u.Mpc**-3).value  ### coming number density of hydrogen atoms in cm^-3
    
    def CHII(self,z):
        return 2.9*((1+z)/6)**-1.1

    def dtdz(self,z):
        return ((-1/(H(z)*(1+z))).to(u.s)).value

    def dnrec_dz_path(self,deltaV,xHII_Field:np.ndarray)->np.ndarray:
        x_HE = 1.08
        # CIGM = self.CHII(z)
        CIGM = 3.0
        nh = self.n_HI(deltaV)*(u.Mpc**-3).to(u.cm**-3)
        Q_HII = xHII_Field
        # alpha_A = 4.2e-13 #cm**3/s
        alpha_B = 2.59e-13  # cm^3/s at 10^4 K
        differential_trans = self.dtdz(self.z)
        return -CIGM*x_HE*alpha_B*nh*Q_HII*(1+self.z)**3 * differential_trans
    
    ### Minihalo term
    def dndmeps(self,M,Mr,deltaV,z):
        deltaL = self.cosmo.delta_L(deltaV,self.z)
        sig1 = self.cosmo.sigma2_interp(M) - self.cosmo.sigma2_interp(Mr)
        del1 = self.cosmo.deltac(z) - deltaL
        return self.cosmo.rhom * (1 + deltaV) / M / np.sqrt(2 * np.pi) * abs(self.cosmo.dsigma2_dm_interp(M)) * del1 / (sig1**(3 / 2)) * np.exp(-del1 ** 2 / (2 * sig1))

    def dEPS_dz(self,M,Mv,deltaV,z):
        return (self.dndmeps(M,Mv,deltaV,z+0.001*z) - self.dndmeps(M,Mv,deltaV,z-0.001*z)) / (0.002*z)

    def Xim(self,m,z):
        A,B,C,D,E,F,G = 4.4 , 0.334, 0.023, 0.199, -0.042, 0,1
        M7 = m/1e7
        F0 = 1.0
        return 1+A*M7**(B+C*np.log10(M7))*(F+G*(1+z)/10)* F0**(D+E*np.log10(F0))
    
    def dNxi_dz(self,m,deltaR,Mv,z):
        m_H = (m_p.to(u.M_sun)).value #M_sun
        return self.Xim(m,z)*m*self.dEPS_dz(m,Mv,deltaR,z)/ m_H * self.cosmo.omegab / self.cosmo.omegam

    
    def Simpson_Nxi_dz(self,delta_R:np.ndarray,M_v:float,z:float) -> np.ndarray:
        try:
            Nxi_arr = np.load(f'.Nxi_dz_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{M_v:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy')
        except FileNotFoundError:
            os.makedirs(f'.Nxi_dz_Interp_init/NxiAtz{self.z:.2f}', exist_ok=True)
            M_min = self.M_J
            M_max = self.M_min
            Mh_interp = np.linspace(np.log10(M_min), np.log10(M_max), 1000)  # log10(M)
            deltaR_grid, m_grid = np.meshgrid(self.deltaV_interp, Mh_interp, indexing='ij')
            mh = self.dNxi_dz(m=10**m_grid, deltaR=deltaR_grid, Mv=M_v, z=z)
            mh[mh > 0] = 0
            integrand = mh * 10**m_grid * np.log(10)  # dN/dlogM -> dN/dM
            # Integrate using Simpson's rule
            Nxi_arr_save = simpson(integrand, x=Mh_interp, axis=1)
            np.save(f'.Nxi_dz_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{M_v:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy', Nxi_arr_save)
            Nxi_arr = np.load(f'.Nxi_dz_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{M_v:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy')
        Simpson_Interp = interp1d(self.deltaV_interp, Nxi_arr, kind='cubic')
        return Simpson_Interp(delta_R)
    
    ### Minihalo Initialization
    def Nxi_Pure(self,Mv,deltaR):
        m_H = (m_p.to(u.M_sun)).value #M_sun
        def Nxi_Pure_diff(m,Mv,deltaR):
            return self.Xim(m,self.z)*m*self.dndmeps(m,Mv,deltaR,self.z)/ m_H * self.cosmo.omegab / self.cosmo.omegam
        mslice = np.logspace(np.log10(self.M_J), np.log10(self.M_min), 12)
        ans = np.zeros_like(deltaR)
        for i in range(len(mslice)-1):
            ans += quad_vec(Nxi_Pure_diff, mslice[i], mslice[i+1],args=(Mv,deltaR), epsrel=1e-6)[0]
        return ans
    
    def Nxi_interp(self, Mv, deltaR):
        try:
            Nxi_arr = np.load(f'.Nxi_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{Mv:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy')
        except FileNotFoundError:
            os.makedirs(f'.Nxi_Interp_init/NxiAtz{self.z:.2f}', exist_ok=True)
            nxi_pure = self.Nxi_Pure(Mv, self.deltaV_interp)
            np.save(f'.Nxi_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{Mv:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy', nxi_pure)
            Nxi_arr = np.load(f'.Nxi_Interp_init/NxiAtz{self.z:.2f}/Nxi_arr_Mv_{Mv:.3f}at_z={self.z:.2f}_A{self.A2byA1}_k{self.kMpc_trans}_alpha{self.alpha}_beta{self.beta}.npy')
        Nxi_interp_Mv = interp1d(self.deltaV_interp, Nxi_arr, kind='cubic')
        return Nxi_interp_Mv(deltaR) 
    
    ### Minihaolo ST_over_PS
    def Nxi_ST(self,z):
        m_H = (m_p.to(u.M_sun)).value #M_sun
        Min = self.cosmo.M_Jeans(z)
        Max = self.cosmo.M_vir(0.61,1e4,z)
        def Nxi_ST_diff(m):
            return (self.Xim(m,z)*1 / m_H * self.cosmo.omegab / self.cosmo.omegam * m * self.cosmo.dndmst(m, z))
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
