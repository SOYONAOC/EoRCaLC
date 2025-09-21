from .from_power import MassFunctions
import numpy as np
from scipy.interpolate import interp1d
import astropy.units as u
from cosfunc import n_H,dtdz
from astropy.constants import m_p
from .special import qion_sb99


class Ion:

    def __init__(self,z,fesc=0.2,A2byA1=0.1,ktrans=200,alpha=2.0,beta=0.0):
        self.cosmo = MassFunctions(A2byA1=A2byA1,kMpc_trans=ktrans,alpha=alpha,beta=beta)
        self.z = z
        self.mH = (self.cosmo.mHu.to(u.M_sun)).value #M_sun
        self.ob0 = self.cosmo.omegab
        self.om0 = self.cosmo.omegam
        self.nH = self.cosmo.nH  #cm^-3
        self.fesc = fesc
        self.deltaV_interp = np.linspace(-0.999, 2, 100)  # delta_V
        self.M_min = self.cosmo.M_vir(0.61, 1e4, z)  # M_sun
        self.M_max = self.cosmo.M_vir(0.61, 2e8, z)
        self.M_J = self.cosmo.M_Jeans(z)  # M_sun
        self.qion = qion_sb99(z)  # ionizing photon per baryon in stars

    ### Source term
    def nion_interp(self,Mv:float,deltaV:np.ndarray)->np.ndarray:
        def diff(m,Mv,deltaR):
            return self.fstar(m)*m*self.cosmo.dndmeps(m,Mv,deltaR,self.z)/ self.mH * self.ob0 / self.om0
        return self._nion_trapez(diff,Mv,deltaV)*self.qion*self.fesc
    
    def nion_st(self,z):
        def diff(m):
            return self.fstar(m)*m*self.cosmo.dndmst(m,z)/ self.mH * self.ob0 / self.om0
        x = np.linspace(np.log10(self.M_min),np.log10(self.M_max),100000)
        y = diff(10**x)*10**x*np.log(10)
        return np.trapezoid(y,x)*self.qion*self.fesc
    
    def nion_ps(self,z):
        def diff(m):
            return self.fstar(m)*m*self.cosmo.dndmps(m,z)/ self.mH * self.ob0 / self.om0
        x = np.linspace(np.log10(self.M_min),np.log10(self.M_max),100000)
        y = diff(10**x)*10**x*np.log(10)
        return np.trapezoid(y,x)*self.qion*self.fesc
    
    ### Minihalo term
    def eps_dz(self,M,Mv,deltaV,z):
        return (self.cosmo.dndmeps(M,Mv,deltaV,z+0.001*z) - self.cosmo.dndmeps(M,Mv,deltaV,z-0.001*z)) / (0.002*z)
    
    def nxi_dz_interp(self, Mv, deltaV):
        def diff(m, Mv, deltaV):
            return self.xim(m,self.z) * m * self.eps_dz(m, Mv, deltaV, self.z) / self.mH * self.ob0 / self.om0
        return self._nxi_trapez(diff, Mv, deltaV,zero_cut=True)
        
    def nxi_interp(self, Mv, deltaV):
        def diff(m, Mv, deltaV):
            return self.xim(m,self.z) * m * self.cosmo.dndmeps(m, Mv, deltaV, self.z) / self.mH * self.ob0 / self.om0
        return self._nxi_trapez(diff, Mv, deltaV)
    
    ### Minihaolo ST_over_PS
    def nxi_st_ini(self,z):
        def diff(m):
            return (self.xim(m,z)*1 / self.mH * self.ob0 / self.om0 * m * self.cosmo.dndmst(m, z))
        x = np.linspace(np.log10(self.cosmo.M_Jeans(z)), np.log10(self.cosmo.M_vir(0.61,1e4,z)), 1000)
        y = diff(10**x)*10**x *np.log(10)
        return np.trapezoid(y,x)
    
    def st_dz(self,m,z):
        return (self.cosmo.dndmst(m,z+0.01) - self.cosmo.dndmst(m,z-0.01)) / 0.02

    def nxi_dz_st(self, z: np.ndarray):
        z = np.asarray(z)
        z_shape = z.shape
        z = z.ravel()  # 展平，便于广播，最后再 reshape

        # 固定大质量网格（覆盖所有可能 z 的质量范围）
        log_m_min_global = 3.0   # 1e3 M☉
        log_m_max_global = 10.0  # 1e10 M☉
        x = np.linspace(log_m_min_global, log_m_max_global, 1000)  # (1000,)
        m = 10**x  # (1000,)

        # 为每个 z 计算实际质量边界
        M_J_z = self.cosmo.M_Jeans(z)      # (N,)
        M_min_z = self.cosmo.M_vir(0.61, 1e4, z)  # (N,)

        # 扩展维度以便广播
        m_grid = m[:, None]      # (1000, 1)
        z_grid = z[None, :]      # (1, N)
        M_J_grid = M_J_z[None, :]    # (1, N)
        M_min_grid = M_min_z[None, :] # (1, N)

        # 计算被积函数（先不管质量范围）
        diff_val = self.xim(m_grid, z_grid) * m_grid * self.st_dz(m_grid, z_grid) / self.mH * self.ob0 / self.om0
        y = diff_val * m_grid * np.log(10)  # (1000, N)

        # 应用质量范围 mask：只在 M_J(z) <= m <= M_min(z) 内有效，其余设为 0
        mask = (m_grid >= M_J_grid) & (m_grid <= M_min_grid)  # (1000, N)
        y = np.where(mask, y, 0.0)

        # 只保留负贡献（可选，但保险）
        y = np.minimum(y, 0)

        # 沿质量轴积分 → (N,)
        result = np.trapezoid(y, x, axis=0)

        return result.reshape(z_shape)

    def nxi_st(self,z):
        zlin = np.arange(20.1, z-0.05, -0.1)
        step = self.nxi_dz_st(zlin)
        return abs(np.trapezoid(step, zlin)) + self.nxi_st_ini(20.1)


    ### IGM term
    def n_HI(self,deltaV):
        return n_H(deltaV).to(u.Mpc**-3).value  ### coming number density of hydrogen atoms in cm^-3
    
    def CHII(self,z):
        return 2.9*((1+z)/6)**-1.1

    def dnrec_dz_path(self,deltaV,xHII_Field:np.ndarray)->np.ndarray:
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
    def xim(self,m,z):
        A,B,C,D,E,F,G = 4.4 , 0.334, 0.023, 0.199, -0.042, 0,1
        M7 = m/1e7
        F0 = 1.0
        return 1+A*M7**(B+C*np.log10(M7))*(F+G*(1+z)/10)* F0**(D+E*np.log10(F0))

    def fstar(self, Mh): # Donnan 25 
        eps0 = 0.16
        Mc = 10**11.7
        beta = 0.9
        gamma = 0.65
        return 2*eps0 * ((Mh / Mc)**-beta + (Mh / Mc)**gamma)**-1


    #### method functions
    def _nxi_trapez(self, diff_func, Mv, deltaV,zero_cut=False):
        x = np.linspace(np.log10(self.M_J), np.log10(self.M_min), 1000)
        deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, x, indexing='ij')
        m_vals = 10**m_grid
        y = diff_func(m_vals, Mv, deltaV_grid) * m_vals * np.log(10)
        if zero_cut:
            y = np.minimum(y,0)
        integrand = np.trapezoid(y, x, axis=1)
        interp_func = interp1d(self.deltaV_interp, integrand, kind='cubic', bounds_error=False, fill_value=0)
        return interp_func(deltaV)
    
    def _nion_trapez(self, diff_func, Mv, deltaV):
        Mh_max = 0.9*min(self.M_max,Mv)
        x = np.linspace(np.log10(self.M_min), np.log10(Mh_max), 1000)
        deltaV_grid, m_grid = np.meshgrid(self.deltaV_interp, x, indexing='ij')
        m_vals = 10**m_grid
        y = diff_func(m_vals, Mv, deltaV_grid) * m_vals * np.log(10)
        integrand = np.trapezoid(y, x, axis=1)
        interp_func = interp1d(self.deltaV_interp, integrand, kind='cubic', bounds_error=False, fill_value=0)
        return interp_func(deltaV)