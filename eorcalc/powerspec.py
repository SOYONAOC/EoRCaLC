"""
From the matter power spectrum to sigma(M) , dsigma/dM and Massfunction.
"""
from dataclasses import dataclass
import os
import numpy as np
import camb
from scipy.interpolate import interp1d
from scipy.integrate import quad,quad_vec
from filelock import FileLock
import astropy.units as u
from scipy.integrate import simpson
import massfunc as mf 


@dataclass
class CosmologyParams:
    """Set Cosmological parameters."""
    h: float = 0.674
    omegam: float = 0.315

    def __post_init__(self):
        self.omegab = 0.0224 * self.h**-2
        self.omegalam = 1 - self.omegam
        self.rhocrit = 2.775366e11 * self.h**2 * u.Msun / u.Mpc**3 
        rhoc = self.rhocrit.value
        self.rhom = rhoc * self.omegam
        self.H0u = 100 * self.h * (u.km * u.s**-1 * u.Mpc**-1)
        self.mHu = 1.6726e-27 * u.kg                        # the mass of a hydrogen Unit: kg
        self.X = 0.75                                       # the mass fraction of hydrogen
        self.nHu = self.rhocrit.to(u.kg/u.cm**3) * self.omegab * self.X / self.mHu    # hydrogen density
        self.nH = self.nHu.value
        self.omegak = 0.0   
        self.omegar = 0.0

class SteepPowerCalculator:

    #here As, ns are from Planck constraints, you should not change them much
    #alpha, beta, A2byA1, kMpc_trans are model parameters; alpha in [-2, 2], beta in [-pi, pi], A2byA1 in (0, 1), kMpc_trans is where you want the boost of power to start.
    # To get the standard case without power boost, you can set A2byA1 = 1 or kMpc_trans = + infinity.
    def __init__(self, As=2.1e-9, ns = 0.965, A2byA1 = 0.1, kMpc_trans = 1200, alpha=2.0, beta=0.0):
        self.alpha = alpha
        self.beta = beta
        self.c_alpha = np.cosh(alpha)
        self.s_alpha = np.sinh(alpha)
        self.c_beta = np.cos(beta)
        self.s_beta = np.sin(beta)
        self.As = As
        self.ns = ns
        self.A2byA1 = A2byA1
        self.kMpc_trans = kMpc_trans
        self.kMpc_pivot = 0.05 ##where As is defined
        eps = (1.-ns)/6.  #\epsilon        
        self.boost = (self.c_alpha-self.c_beta*self.s_alpha)**2 + (self.s_beta*self.s_alpha)**2

        self.HbyMp = np.pi*np.sqrt(As*eps*8./self.boost)  #H/M_p

        self.A1byMp3 = np.sqrt(eps*18.)*self.HbyMp**2
        self.A2byMp3 = self.A1byMp3*A2byA1

        self.norm = (3.*self.HbyMp**3/(2.*np.pi*self.A2byMp3))**2
        self.coef  = 3.*(1. - self.A2byA1)
        self.c1 = (self.c_alpha-self.c_beta*self.s_alpha)
        self.c2 = self.s_beta*self.s_alpha
        self.c3 = (self.c_alpha+self.c_beta*self.s_alpha)
        self.mp_prepared = False
                


    def PrimordialPower(self, kMpc:list):
        k = kMpc/self.kMpc_trans
        ksq = k**2
        if(k < 0.03): 
            kjy = -1./3.+ ksq*( (2./35.)*ksq -2./15. )
            kjj = k*ksq/9.*(1.-ksq/5.)
        else:
            kjy = ((k**2-1.)*np.sin(2*k)+2*k*np.cos(2*k))/(2.*k**3)
            kjj = (k*np.cos(k)-np.sin(k))**2/k**3
        return self.norm*(kMpc/self.kMpc_pivot)**(self.ns-1.)*((self.c1*(1.+self.coef*kjy) + self.c2*self.coef*kjj)**2 + (self.c2*(1.+self.coef*kjy) + self.c3*self.coef*kjj)**2)

   # Transfer function from BBKS to primordial power
    def BBKS_trans(self, x:float):
        return  np.log(1.+0.171*x)**2 *x**(self.ns-2.) / np.sqrt( 1.+ x*(0.284+x*(1.3924+x*(0.0635212+x*0.057648))))

    def Prepare_MatterPower(self, H0 = 67.5, ombh2 = 0.022, omch2 = 0.12, kMpc_max=10., zmax = 100.):
        nz = 100 #number of steps to use for the radial/redshift integration
        kmax=min(kMpc_max, 30.)  #kmax to use
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
        pars.InitPower.set_params(As = self.As, ns=self.ns)
        results= camb.get_background(pars)
        self.mpi =  camb.get_matter_power_interpolator(pars, nonlinear=False, hubble_units=False, k_hunit=False, kmax=kmax, zmax=zmax)
        self.k_eq = 0.073*(omch2+ombh2)
        self.kmax = kmax
        self.norm_edge = self.BBKS_trans(self.kmax/self.k_eq)
        self.mp_prepared = True



    def MatterPower(self, kMpc:list, z): #kMpc can be an ordered array; z is the redshift
        if(not self.mp_prepared):
            self.Prepare_MatterPower()
        if(isinstance(kMpc, float)):
            if(kMpc <= self.kmax):
                pkraw = self.mpi.P(z, kMpc)
            else:
                pkraw = self.BBKS_trans(kMpc/self.k_eq)/self.norm_edge * self.mpi.P(z, self.kmax)
            return pkraw * (self.PrimordialPower(kMpc)/(self.As*(kMpc/self.kMpc_pivot)**(self.ns-1.)))
        else:
            nk = len(kMpc)
            Pk = np.empty(nk)
            for i in range(nk):
                if(kMpc[i] <= self.kmax):
                    pkraw = self.mpi.P(z, kMpc[i])
                else:
                    pkraw = self.BBKS_trans(kMpc[i]/self.k_eq)/self.norm_edge*self.mpi.P(z, self.kmax)
                Pk[i] = pkraw * (self.PrimordialPower(kMpc[i])/(self.As*(kMpc[i]/self.kMpc_pivot)**(self.ns-1.)))
            return Pk

class MassFunctions(CosmologyParams):
    def __init__(self, A2byA1=0.1, kMpc_trans=2e2, alpha=2.0, beta=0.0,h=0.674, omegam=0.315):
        super().__init__(h=h, omegam=omegam)
        self.cosmo = mf.SFRD(h=h,omegam=omegam)
        self.spc = SteepPowerCalculator(As = 2.1e-9, ns = 0.965, A2byA1 = A2byA1, kMpc_trans = kMpc_trans, alpha = alpha, beta = beta) 
        self.A2byA1, self.kMpc_trans, self.alpha, self.beta = A2byA1, kMpc_trans, alpha, beta
        self.dsigma2_dm_interpolation_completed = False
        self.pk_interpolation_completed = False
        self.sigma2_interpolation_completed = False
        self.m_interp_range = np.logspace(0.0,18.0,1000)
        self.ps_interp_range = np.logspace(-4.0,6.0,50000)
    

    def PowerSpectrum_Interp_Set(self):
        filename = f'.ps_init_out/Pk_interp_A{self.A2byA1}_K{self.kMpc_trans}_Alpha{self.alpha}_Beta{self.beta}.npy'
        lockfile = filename + '.lock'
        try:
            Pk_arr = np.load(filename)
        except FileNotFoundError:
            os.makedirs('.ps_init_out', exist_ok=True)
            with FileLock(lockfile):
                if not os.path.exists(filename):
                    self.spc.Prepare_MatterPower(H0=67.4, ombh2=0.022, omch2=0.12, kMpc_max=self.ps_interp_range[-1])
                    Pk = self.spc.MatterPower(self.ps_interp_range, z=0.0)
                    np.save(filename, Pk)
                Pk_arr = np.load(filename)
        self.pk0_interpolation = interp1d(np.log10(self.ps_interp_range), np.log10(Pk_arr), kind='cubic')
        self.pk_interpolation_completed = True

    def Sigma2_Interp_Set(self):
        filename = f'.ps_init_out/Sigma2_interp_A{self.A2byA1}_K{self.kMpc_trans}_Alpha{self.alpha}_Beta{self.beta}.npy'
        lockfile = filename + '.lock'
        try:
            sig_arr = np.load(filename)
        except FileNotFoundError:
            os.makedirs('.ps_init_out', exist_ok=True)
            with FileLock(lockfile):
                if not os.path.exists(filename):
                    sig = self.sigma2(self.m_interp_range)
                    np.save(filename, sig)
                sig_arr = np.load(filename)
        self.sigma2m_interpolation = interp1d(np.log10(self.m_interp_range), np.log10(sig_arr), kind='cubic')
        self.sigma2_interpolation_completed = True

    def Dsigma2dm_Interp_Set(self):
        filename = f'.ps_init_out/Dsigma2_dm_interp_A{self.A2byA1}_K{self.kMpc_trans}_Alpha{self.alpha}_Beta{self.beta}.npy'
        lockfile = filename + '.lock'
        try:
            dsig2dm_arr = np.load(filename)
        except FileNotFoundError:
            os.makedirs('.ps_init_out', exist_ok=True)
            with FileLock(lockfile):
                if not os.path.exists(filename):
                    dsig2dm = np.log10(-self.dsigma2_dm(self.m_interp_range))
                    np.save(filename, dsig2dm)
                dsig2dm_arr = np.load(filename)
        self.dsigma2_dm_interpolation = interp1d(np.log10(self.m_interp_range), dsig2dm_arr, kind='cubic')
        self.dsigma2_dm_interpolation_completed = True

    def ps_interp(self, kMpc):
        if not self.pk_interpolation_completed:
            self.PowerSpectrum_Interp_Set()
        return 10**self.pk0_interpolation(np.log10(kMpc))

    def sigma2_interp(self,M):
        if not self.sigma2_interpolation_completed:
            self.Sigma2_Interp_Set()
        return 10**self.sigma2m_interpolation(np.log10(M))
    
    def dsigma2_dm_interp(self,M):
        if not self.dsigma2_dm_interpolation_completed:
            self.Dsigma2dm_Interp_Set()
        return -10 ** self.dsigma2_dm_interpolation(np.log10(M))

    def dsigma2_dk(self,M,k):
        r=(3.0*M/(4.0*np.pi*self.rhom))**(1./3.) 
        x=r*k
        w=3.0*(np.sin(x)-x*np.cos(x))/x**3.0
        return 4.0*np.pi*k**2.0/(2.0*np.pi)**3.0 * self.ps_interp(k) * w*w

    def dsigma2_dlnk(self,lnk,M):
        k=np.exp(lnk)
        return k*self.dsigma2_dk(M,k)
    
    def sigma2(self,M):
        precision = 1e-6
        k_grid_1 = np.logspace(-4,0.9,20)
        k_grid_2 = np.logspace(1,5.3,200)
        k_grid = np.concatenate([k_grid_1, k_grid_2])
        int_ans = np.zeros_like(M)
        for i in range(len(k_grid)-1):
            int_ans += quad_vec(self.dsigma2_dlnk,np.log(k_grid[i]),np.log(k_grid[i+1]), args=(M,),epsrel=precision,limit=200)[0]
        return int_ans

    def sigma2_tarpz(self, M):
        k_grid_1 = np.logspace(-4, 0.999, 500)  # 增加点数提高精度
        k_grid_2 = np.logspace(1, 5.3, 2000)
        k_grid = np.concatenate([k_grid_1, k_grid_2])
        lnk_grid = np.log(k_grid)
        M_mesh, lnk_mesh = np.meshgrid(M, lnk_grid, indexing='ij')
        integrand_values = self.dsigma2_dlnk(lnk_mesh.T, M_mesh.T).T
        return np.trapezoid(integrand_values, x=lnk_grid, axis=1)

    def sigma2_simpson(self, M):
        k_grid_1 = np.logspace(-4, 0.999, 501)  
        k_grid_2 = np.logspace(1, 5.3, 2001)
        k_grid = np.concatenate([k_grid_1, k_grid_2])
        lnk_grid = np.log(k_grid)
        M_mesh, lnk_mesh = np.meshgrid(M, lnk_grid, indexing='ij')
        integrand_values = self.dsigma2_dlnk(lnk_mesh.T, M_mesh.T).T
        return simpson(integrand_values, x=lnk_grid, axis=1)

    def dsigma2_dlnk_dm(self,lnk,M):
        k=np.exp(lnk)
        r=(3.0*M/(4.0*np.pi*self.rhom))**(1./3.)
        x=r*k
        w=3.0*(np.sin(x)-x*np.cos(x))/(x)**3.0
        dw_dx=(9.0*x*np.cos(x)+np.sin(x)*(3.0*x**2.0-9.0))/(x)**4.0
        dw_dm=dw_dx*k/(4.0*np.pi*r**2.0*self.rhom)
        return 4.0*np.pi*k**3.0/(2.0*np.pi)**3.0*self.ps_interp(k)*2.0*w*dw_dm

    def dsigma2_dm(self,M):
        precision = 1e-6
        k_grid_1 = np.logspace(-4,0.999,20)

        k_grid_2 = np.logspace(1,5.3,200)
        k_grid = np.concatenate([k_grid_1, k_grid_2])
        int_ans = np.zeros_like(M)
        for i in range(len(k_grid)-1):
            int_ans += quad_vec(self.dsigma2_dlnk_dm,np.log(k_grid[i]),np.log(k_grid[i+1]), args=(M,),epsrel=precision,limit=200)[0]
        return int_ans

    def dsigma2_dm_simpson(self, M):
        k_grid_1 = np.logspace(-4, 0.999, 501)  # 奇数点数
        k_grid_2 = np.logspace(1, 5.3, 2001)
        k_grid = np.concatenate([k_grid_1, k_grid_2])
        lnk_grid = np.log(k_grid)
        
        # 向量化计算
        M_mesh, lnk_mesh = np.meshgrid(M, lnk_grid, indexing='ij')
        integrand_values = self.dsigma2_dlnk_dm(lnk_mesh.T, M_mesh.T).T
        
        return simpson(integrand_values, x=lnk_grid, axis=1)    

    def dsigma2_dm_trapz(self, M):
        k_grid_1 = np.logspace(-4, 0.999, 1000)  # 任意点数
        k_grid_2 = np.logspace(1, 5.3, 5000)
        k_grid = np.concatenate([k_grid_1, k_grid_2])
        lnk_grid = np.log(k_grid)
        
        # 向量化计算
        M_mesh, lnk_mesh = np.meshgrid(M, lnk_grid, indexing='ij')
        integrand_values = self.dsigma2_dlnk_dm(lnk_mesh.T, M_mesh.T).T
        
        return np.trapezoid(integrand_values, x=lnk_grid, axis=1)
    
    def deltac(self, z):
        return self.cosmo.deltac(z)
    
    def dndmps(self, m, z):
        sigm = np.sqrt(self.sigma2_interp(m))
        dsig_dm = abs(self.dsigma2_dm_interp(m)) / (2.0 * sigm)
        return np.sqrt(2.0 / np.pi) * self.rhom / m * self.deltac(z) / sigm**2 * dsig_dm * np.exp(-self.deltac(z) ** 2 / (2 * sigm ** 2))

    def dndmst(self,M,z):
        Ast=0.353
        ast2=0.73
        Pst=0.175
    #   ST parameters from Jenkins et al. 2001,0005260
        sigma = np.sqrt(self.sigma2_interp(M))
        nu=self.deltac(z)/sigma
        nup=np.sqrt(ast2)*nu
        dsigmadm=self.dsigma2_dm_interp(M)/(2*sigma)
        return -self.rhom/M*(dsigmadm/sigma)*(2*Ast)*(1+1.0/(nup)**(2*Pst))*(nup**2/(2*np.pi))**(1./2.)*np.exp(-nup**2/2)

    def fcoll_st(self, z,Mmin, Mmax=1e16,n=2000):
        logm_max = np.log10(Mmax)
        logm_min = np.log10(Mmin)
        def diff(m,z):
            return self.dndmst(m,z)*m
        x = np.linspace(logm_min, logm_max, n)
        y = diff(10**x, z)*10**x*np.log(10)
        integral = np.trapezoid(y, x)        
        return integral / self.rhom

    def delta_L(self, deltar,z):
        return (1.68647 - 1.35 / (1 + deltar) ** (2 / 3) - 1.12431 / (1 + deltar) ** (1 / 2) + 0.78785 / (1 + deltar) ** (0.58661)) / self.cosmo.Dz(z)

    def dndmeps(self, M, Mr, deltar, z):
        delta_L = self.delta_L(deltar,z)
        sig1 = self.sigma2_interp(M) - self.sigma2_interp(Mr)
        del1 = self.cosmo.deltac(z) - delta_L
        return self.rhom * (1 + deltar) / M / np.sqrt(2 * np.pi) * abs(self.dsigma2_dm_interp(M)) * del1 / sig1 ** (3 / 2) * np.exp(-del1 ** 2 / (2 * sig1))
    
    def M_Jeans(self, z):
        return 5.73e3*(self.omegam*self.h**2/0.15)**(-1/2) * (self.omegab*self.h**2/0.0224)**(-3/5) * ((1+z)/10)**(3/2)

    def Delta_cc(self,z):
        d=self.cosmo.omegam_z(z)-1.0
        return 18*np.pi**2+82.0*d-39.0*d**2

    def M_vir(self,mu,Tvir,z):
        a1=(self.cosmo.omegam_z(z)*self.Delta_cc(z)/(18*np.pi**2))**(-1.0/3.0)
        a2=a1*(mu/0.6)**(-1.0)*((1.0+z)/10)**(-1.0)/1.98e4*Tvir
        return a2**(3.0/2.0)*1e8/self.h