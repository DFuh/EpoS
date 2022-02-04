'''
storage auxilliary file:
- class: STRG (storage)
- class: HFE_model

# TODO:
    - hardcoded parameters !

'''
import numpy as np
import scipy.integrate as scpi
import scipy.optimize as scpo
import scipy.interpolate as scpin

# import HFE_model_v0 as hfe
class STRG():
    '''Object representing a cavern storage.
    Inputs : height [m], Volume [m^3], initial pressure [Pa], initial Temperature [K]'''

    # Cavern infromation out of excelfile
    # a = os.path.dirname(__file__)

    def __init__(self,simu_obj): #height=300,V=500e3,p0=40e5,T0=317.15):
        # ci = add_cavern_information(read_data("Information_cavern.xlsx"))

        # general specifications
        # Reference : Hydrogen Science & Engineering, Chapter: Thermodynamics of Pressurized Gas Storage(p.601-628), Tietze et al., 2016
        # salt temperature [K]
        self.R = simu_obj.pec.R
        self.M_H2 = simu_obj.pec.M_H2
        self.T_salt = simu_obj.pstrg.T_amb # ci['T_salt']
        self.par = simu_obj.pstrg

        self.p_max = self.par.pressure_max
        self.p_min = self.par.pressure_min
        # heat convection coefficient between storage and and salt surface [W/(m^2K)]
        self.alpha_storage_salt = simu_obj.pstrg.alpha_storage # 133    # ci['alpha_storage_salt']
        # heat conduction coefficient of the salt [W/(mK)]
        self.lambda_salt = simu_obj.pstrg.lambda_wall # 5.5           # ci['lambda_salt']
        # thikness of salt [m]
        self.d_salt = simu_obj.pstrg.d_wall #ä2                  # ci['d_salt']
        # --------------------------------------------------------------------
        #initialize EOS model
        # self.EoS = hfe.HFE_EOS()
        self.EoS = HFE_EOS()
        # initialize Volume and geometry of cavern assuming a cylinder shape
        self.V = simu_obj.pstrg.capacity_volume
        self.L = simu_obj.pstrg.height
        if not simu_obj.pstrg.diameter_inner:
            self.d_i = np.sqrt(4*self.V/(self.L*np.pi))
        self.d_o = self.d_i+2*self.d_salt
        self.cap_m = simu_obj.pstrg.capacity_mass
        if not self.cap_m:
            self.cap_m = self.clc_cap_cavern(p_min=self.par.pressure_min,
                                     p_max=self.par.pressure_max,
                                     T=self.par.TN,
                                     V_cvrn=self.V)
            simu_obj.prms['parameters_strg']['strg_0']['capacity_mass']['value'] = self.cap_m

        print('STRG: cap_m=', self.cap_m)
        # initialize Temperature and pressure and coresponding mass, density and internal energy if given
        # if not(p0==None) and not(T0==None):
        self.T = self.par.T0
        self.p = self.par.p0
        self.rho = self.EoS.density(self.T,self.p)
        self.m = self.rho*self.V
        self.u = self.EoS.internal_energy(self.T,self.rho)

        self.m_min = self.V* self.EoS.density(self.par.T0,self.p_min)
        self.m_max = self.V* self.EoS.density(self.par.T0,self.p_max)
        print('Strorage capacity (from EOS): ', self.m_max-self.m_min)
        print('Strorage capacity (simple clc): ', self.cap_m)
        print(f'Strorage m0 (@ p = {self.p}): {self.m}')
        # self.m_max = self.cap_m + self.m_min

    def Heat(self):
        '''function returning Heat flow [W] between cavern storage and surrounding salt'''
        R_cavern = 1/(np.pi*self.L)*(1/(self.alpha_storage_salt*self.d_i)+np.log(self.d_o/self.d_i)/(2*self.lambda_salt))
        return 1/R_cavern*(self.T_salt-self.T)

    def Heat2(self,T_store):
        '''function returning Heat flow [W] between cavern storage and surrounding salt taking storage temperature [K] as input'''
        R_cavern = 1/(np.pi*self.L)*(1/(self.alpha_storage_salt*self.d_i)+np.log(self.d_o/self.d_i)/(2*self.lambda_salt))
        return 1/R_cavern*(self.T_salt-T_store)

    def set_initial_conditions(self,p0,T0):
        '''function to set initial conditions for T0 and p0 to reset values from __init__ method and recompute rho and m'''
        self.T = T0
        self.p = p0
        self.rho = self.EoS.density(self.T,self.p)
        self.m = self.rho*self.V
        self.u = self.EoS.internal_energy(self.T,self.rho)
        return


    def clc_state_ostp(self, T_in, t_span, dm_in, dm_out, heatloss_wall=True):
        '''
        clc storage state for one timestep
        '''
        # T_in = ? # T of H2 entering system (T of gas-treatment)

        results = scpi.solve_ivp(self.ODE,t_span,
                                [self.u,self.rho],
                                args=[T_in, dm_in, dm_out, heatloss_wall],
                                method='RK45', max_step=np.inf,
                                t_eval=t_span, dense_output=True)

        '''
        following 8 lines adopted from Henry's original cavern.py'
        '''
        # rho = results.y[1,:]    #kg/m^3
        # u = results.y[0,:]      #J/kg
        # T = np.array([HFE_EOS().temperature(x,y) for x,y in zip(rho,u)])    #K
        # p = np.array([HFE_EOS().pressure(x,y) for x,y in zip (T,rho)])    #Pa
        # t = results.t           #h
        # m = rho*self.V*1e-3     #ton
        # U = u*m*1e-9            #TJ
        # Q = np.array([self.Heat2(x)*1e-6 for x in T]) #MW
        rho = results.y[1,-1]    # // in kg/m^3
        u = results.y[0,-1]      # // in J/kg
        T = self.EoS.temperature(rho,u) # in K
        p = self.EoS.pressure(T,rho)    #Pa
        t = results.t           # h
        m = rho*self.V          # // in kg
        U = u*m *1e-3                 # // in kJ
        Q = self.Heat2(T)*1e-3 # // in kW
        return (rho, u, T, p, t, m, U, Q)


    def clc_state_mstp(self, T_in, p_in, t_span, dm_in, dm_out, max_step=np.inf,
                        heatloss_wall=True):
        '''
        clc for multiple time steps at once
            dm_in: massflow to storage
            dm_out: massflow from storage

        '''
        step_in = 10
        m_dot = np.vstack((np.arange(0,len(dm_in)*step_in,step=step_in),
                           dm_in,dm_out)).T

        results = scpi.solve_ivp(self.ODE_mstp,t_span,
                                [self.u,self.rho],
                                args=[T_in, m_dot, heatloss_wall],
                                method='RK45', max_step=max_step,
                                t_eval=m_dot[:,0], dense_output=True)

        rho = results.y[1,:]    #kg/m^3
        u = results.y[0,:]      #J/kg
        T = np.array([self.EoS.temperature(x,y) for x,y in zip(rho,u)])    #K
        p = np.array([self.EoS.pressure(x,y) for x,y in zip (T,rho)])    #Pa
        t = results.t           #h (?)
        m = rho*self.V          # kg    ##*1e-3     #ton
        U = u*m                 # J       ##*1e-9            #TJ
        Q = np.array([self.Heat2(x)*1e-3 for x in T]) #kW
        P_cmp = np.array([self.EoS.isothermal(T_in,[p_in,x])[0] for x in p])*dm_in *1e-3  #kW (isothermal returns spec. work in J/kg)

        m_dot_in = scpin.interp1d(m_dot[:,0],m_dot[:,1],kind='nearest',fill_value=(m_dot[0,1],m_dot[-1,1]),bounds_error=False)
        m_dot_out = scpin.interp1d(m_dot[:,0],m_dot[:,2],kind='nearest',fill_value=(m_dot[0,2],m_dot[-1,2]),bounds_error=False)

        m_dot_in_act = np.zeros(len(m_dot))
        for i in range(len(t)):
            m_dot_in_act[i] = m_dot_in(t[i])
        return (rho, u, T, p, t, m, U, Q, P_cmp), m_dot_in_act


    # differntial equation describing filling/unfilling process
    def ODE(self,t,Y, T_in, dm_in, dm_out, c):
        '''system of ODE's representing mass and energy balance of cavern storage'''
        # Reference : Hydrogen Science & Engineering, Chapter: Thermodynamics of Pressurized Gas Storage(p.601-628), Tietze et al., 2016

        #assign variables
        self.u = Y[0]
        self.rho = Y[1]

        #interpolate massflow for timestep
        # print('m_dot[:,0],m_dot[:,1],m_dot[0,1],m_dot[-1,1]:', m_dot[:,0],m_dot[:,1],m_dot[0,1],m_dot[-1,1])

        # Dismissed interpolation of massflow, since values of
        # massflows are considered constant during timestep

        # m_dot_in = interp1d(m_dot[:,0],m_dot[:,1],kind='nearest',fill_value=(m_dot[0,1],m_dot[-1,1]),bounds_error=False)
        # m_dot_out = interp1d(m_dot[:,0],m_dot[:,2],kind='nearest',fill_value=(m_dot[0,1],m_dot[-1,1]),bounds_error=False)


        # calculate state variables from HFE_EOS functions
        rho_in = self.EoS.density(T_in,self.p)
        self.m = self.rho*self.V
        self.T = self.EoS.temperature(self.rho,self.u)
        self.p = self.EoS.pressure(self.T,self.rho)

        # energy and massbalance in terms of derivative of specific internal energy and density
        dudt = (((self.EoS.enthalpy(T_in,rho_in) * dm_in
                  # (((self.EoS.enthalpy(T_in,rho_in) * dm_in(t)
                  -self.EoS.enthalpy(self.T,self.rho) * dm_out)
                 # -self.EoS.enthalpy(self.T,self.rho) * dm_out(t))
                 -self.u * (dm_in - dm_out)
                 # -self.u * (dm_in(t) - dm_out(t))
                # +float(c)*self.Heat()*3600)/self.m  )
                + float(c)*self.Heat()) / self.m  )

        drhodt = (dm_in - dm_out) / self.V
        # drhodt = (dm_in(t) - dm_out(t)) / self.V

        return [dudt,drhodt]

    # differntial equation describing filling/unfilling process
    def ODE_mstp(self,t,Y,T_in,m_dot,c):
        '''system of ODE's representing mass and energy balance of cavern storage'''
        # Reference : Hydrogen Science & Engineering, Chapter: Thermodynamics of Pressurized Gas Storage(p.601-628), Tietze et al., 2016

        #assign variables
        self.u = Y[0]
        self.rho = Y[1]

        #interpolate massflow for timestep
        # print('m_dot[:,0],m_dot[:,1],m_dot[0,1],m_dot[-1,1]:', m_dot[:,0],m_dot[:,1],m_dot[0,1],m_dot[-1,1])
        m_dot_in = scpin.interp1d(m_dot[:,0],m_dot[:,1],kind='nearest',fill_value=(m_dot[0,1],m_dot[-1,1]),bounds_error=False)
        m_dot_out = scpin.interp1d(m_dot[:,0],m_dot[:,2],kind='nearest',fill_value=(m_dot[0,2],m_dot[-1,2]),bounds_error=False)


        # calculate state variables from HFE_EOS functions
        rho_in = self.EoS.density(T_in,self.p)
        self.m = self.rho*self.V
        self.T = self.EoS.temperature(self.rho,self.u)
        self.p = self.EoS.pressure(self.T,self.rho)

        # energy and massbalance in terms of derivative of specific internal energy and density
        dudt = (((self.EoS.enthalpy(T_in,rho_in)*m_dot_in(t)
                  -self.EoS.enthalpy(self.T,self.rho)*m_dot_out(t))
                 -self.u*(m_dot_in(t)-m_dot_out(t))
                 +float(c)*self.Heat())/self.m  )
        drhodt = (m_dot_in(t)-m_dot_out(t))/self.V

        return [dudt,drhodt]

# ==============================================================================
# ==============================================================================

    def clc_cap_cavern(self, p_min=None, p_max=None, T=None, V_cvrn=None):
        '''
        calc. gravimetric capacity of cavern/storage
        ## adopted from powersyn_v55
        '''
        # Crotogino: Huntor CAES -> V=310000m³ // p_min =43 bar, p_max = 70 bar

        #=========================================================================
        #CAUTION: WesPe -> 4,8 kt H2 @ 600000 m³ cvrn (probably 60 ...190 bar)
        #=========================================================================
        # if not p_max:
        #     p_max = 70 *1e5     # // in Pa          z_H2 (313K) -> 1.0409
        # if not p_min:
        #     p_min = 43 *1e5     # // in Pa              z_H2 (313K) -> 1.0249
        # if not vol:
        #    V_cvrn = 310000 # // in m³
        #else:
        #    V_cvrn = vol
        # R = 8.314       # // in J/K mol
        # M_H2 = 2.01588 * 1e-3 # // in kg/mol
        # T = 313
        z_min = H2_compressibility_factor(p_min,T) #H2_compressibility_factor(self,p_in,T):
        z_max = H2_compressibility_factor(p_max,T) #H2_compressibility_factor(self,p_in,T):
        #m_cap_cvrn = M_H2 * ((p_max - p_min)/ (T*z)) * (V_cvrn)/R
        m_cap_cvrn = self.M_H2 * (p_max/ (T*z_max) - p_min/(z_min*T)) * (V_cvrn)/self.R
        # print('m_cap_cvrn: ', m_cap_cvrn)
        return m_cap_cvrn

def H2_compressibility_factor(p_in,T):
        '''
        calc. compressibility factor of H2
        Based on Zheng 2016: Standardized equation for hydrogen gas
        compressibility factor for fuel consumption applications
        -----
        0.1 to 100 MPa
        270 to 500 K
        '''
        # coeff = (i,j)
        p = p_in/1e6
        # fitting coefficients based on table 1 in Zheng2016

        coeff = np.array([
                        [1.00018,        -0.0022546,     0.01053,       -0.013205],
                        [-0.00067291,    0.028051,       -0.024126,     -0.0058663],
                        [0.000010817,    -0.00012653,    0.00019788,    0.00085677],
                        [-1.4368*1e-7,   1.2171*1e-6,    7.7563*1e-7,   -1.7418*1e-5],
                        [1.2441*1e-9,    -8.965*1e-9,    -1.6711*1e-8,  1.4697*1e-07],
                        [-4.4709*1e-12,  3.0271*1e-11,   6.3329*1e-11,  -4.6974*1e-10]
                        ])
        # calc based on eq. 1 in Zheng2016
        z=0
        for i in range(6): #i -> 1 -6
            for j in range(4):
                z += coeff[i,j] * p**(i) *(100/T)**(j)

        return z
# ==============================================================================
# ==============================================================================
# ==============================================================================

def fill_strg(m0, m_max, m_dot_diff, time_step):

    dm_resid = m_dot_diff * time_step # kg/s*s -> kg

    m_act = np.zeros(len(dm_resid))
    dm_grid = np.zeros(len(dm_resid))
    dm_cvrn = np.zeros(len(dm_resid))
    m_act[0] = m0


    '''
    dm <0 -> purchase of H2         # dm -> absolute mass in kg
    dm >0 -> feed to cavern/grid

    dm_grid >0 feed to grid
    dm_grid <0 purchase from grid
    '''
    #print('m_grid: ', m_grid)
    # print('m0: ', m_act[0])
    # print('m_max: ', m_max)
    for i,dm in enumerate(dm_resid):
        if i >0:

            #    print('fill_cvrn: dm', dm)
            if (m_act[i-1] + dm) > m_max:
                #print('res: ', m_max-(m_act[i-1] + dm))
                dm_grid[i] = (m_act[i-1] + dm)-m_max
                dm_cvrn[i] = dm-dm_grid[i]
                m_act[i]   = m_max #m_act[i-1]

            elif m_act[i-1] + dm <0:
                dm_grid[i] = (m_act[i-1]+dm)
                dm_cvrn[i] = dm-dm_grid[i]
                m_act[i]   = m_act[i-1]
            else:
                dm_grid[i] = 0
                dm_cvrn[i] = dm-dm_grid[i]
                m_act[i]   = m_act[i-1]+dm
            # if i < 5:
                # print('fill_cvrn: m[i-1] dm dt m[i]', m_act[i-1], dm, time_step, m_act[i])
    return m_act[:], dm_grid[:]/time_step, dm_cvrn[:]/time_step
# ==============================================================================
# ==============================================================================
# ==============================================================================

class HFE_EOS:
    '''class for the reduced Helmholtz free energy Equation of state model to compute real gas thermodynamic properties'''
    # Thermodynamic constants
    # Reference https://webbook.nist.gov/cgi/cbook.cgi?ID=1333740
    R = 8.314e3 #J/(kmolK)
    T_crit = 33.18 #K
    p_crit = 13e5 #bar
    M = 2.01588 #kg/kmol
    rho_crit =15.4*M #kg/m^3

    # Parameters for Real Gas Component of Reduced Helmholtz Free Energy for hydrogen
    # Reference: A thermodynamic analysis of refueling of a hydrogen tank , Yang, 2009

    a_k = np.array([-1.4579856475, 1.888076782, 1.616, -0.4117, -0.792, 0.758, 1.217])
    b_k = np.array([0, 0, -16.0205159149, -22.6580178006, -60.0090511389, -74.9434303817, -206.9392065168])
    N_i = np.array([-6.93643, 0.01, 2.1101, 4.52059, 0.732564, -1.34086,
                    0.130985,-0.777414, 0.351944,-0.0211716, 0.0226312,
                    0.032187, -0.0231752, 0.0557346])
    t_i = np.array([0.6844, 1, 0.989, 0.489, 0.803, 1.1444, 1.409, 1.754,
                    1.311, 4.187, 5.646, 0.791, 7.249, 2.986])
    d_i = np.array([1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1])
    p_i = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1])
    phi_i = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, -1.685, -0.489, -0.103, -2.506, -1.607])
    beta_i = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, -0.171, -0.2245, -0.1304, -0.2785, -0.3967])
    gamma_i = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7164, 1.3444, 1.4517, 0.7204, 1.5445])
    D_i = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1.506, 0.156, 1.736, 0.67, 1.662])

    #Helmholtz reduced energy functions and their derivatives
    #Reference: A thermodynamic analysis of refueling of a hydrogen tank , Yang, 2009
    def alpha_0(self,T,rho):
        '''function returning alpha_0, taking Temperature[K] and density [kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return (np.log(delta)+1.5*np.log(tau)+self.a_k[0]+self.a_k[1]*tau+
                np.sum(self.a_k[2:] * np.log(1-np.exp(self.b_k[2:]*tau))))

    def alpha_r(self,T,rho):
        '''function returning alpha_r, taking Temperature[K] and density [kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return (np.sum(self.N_i[:7]*delta**self.d_i[:7]*tau**self.t_i[:7])+
                np.sum(self.N_i[7:9]*delta**self.d_i[7:9]*tau**self.t_i[7:9]*np.exp(-delta**self.p_i[7:9]))+
                np.sum(self.N_i[9:]*delta**self.d_i[9:]*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)))

    def alpha_0_diff_tau(self,T,rho):
        '''function returning first partial derivative (wrt:tau) of alpha_0, taking Temperature[K] and density [kg/m^3] as input'''
        tau = self.T_crit/T
        return 1.5/tau+self.a_k[1]-np.sum(self.a_k[2:]*self.b_k[2:]*np.exp(self.b_k[2:]*tau)/(1-np.exp(self.b_k[2:]*tau)))

    def alpha_r_diff_tau(self,T,rho):
        '''function returning first partial derivative (wrt:tau) of alpha_r, taking Temperature[K] and density [kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return (np.sum(self.N_i[:7]*delta**self.d_i[:7]*self.t_i[:7]*tau**(self.t_i[:7]-1))+
                np.sum(self.N_i[7:9]*delta**self.d_i[7:9]*self.t_i[7:9]*tau**(self.t_i[7:9]-1)*np.exp(-delta**self.p_i[7:9]))+
                np.sum(self.N_i[9:]*delta**self.d_i[9:]*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*
                       (self.t_i[9:]/tau+2*self.beta_i[9:]*(tau-self.gamma_i[9:]))))

    def alpha_r_diff_delta(self,T,rho):
        '''function returning first partial derivative (wrt:delta) of alpha_r, taking Temperature[K] and density [kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return (np.sum(self.N_i[:7]*self.d_i[:7]*delta**(self.d_i[:7]-1)*tau**self.t_i[:7])+
                np.sum(self.N_i[7:9]*tau**self.t_i[7:9]*np.exp(-delta**self.p_i[7:9])*delta**self.d_i[7:9]*(self.d_i[7:9]/delta-self.p_i[7:9]*delta**self.p_i[7:9]))+
                np.sum(self.N_i[9:]*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*delta**self.d_i[9:]*
                       (self.d_i[9:]/delta+2*self.phi_i[9:]*(delta-self.D_i[9:]))))

    def alpha_0_diff_tau_tau(self,T,rho):
        '''function returning second partial derivative (wrt:tau,tau) of alpha_0, taking Temperature[K] and density [kg/m^3] as input'''
        tau = self.T_crit/T
        return -1.5/tau**2-np.sum(self.a_k[2:]*(self.b_k[2:]**2*np.exp(self.b_k[2:]*tau)/((1-np.exp(self.b_k[2:]*tau))**2)))

    def alpha_r_diff_tau_tau(self,T,rho):
        '''function returning second partial derivative (wrt:tau,tau) of alpha_r, taking Temperature[K] and density [kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return (np.sum(self.N_i[:7]*delta**self.d_i[:7]*self.t_i[:7]*(self.t_i[:7]-1)*tau**(self.t_i[:7]-2))+
                np.sum(self.N_i[7:9]*delta**self.d_i[7:9]*self.t_i[7:9]*(self.t_i[7:9]-1)*tau**(self.t_i[7:9]-2)*np.exp(-delta**self.p_i[7:9]))+
                np.sum(self.N_i[9:]*delta**self.d_i[9:]*self.t_i[9:]*tau**(self.t_i[9:]-1)*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*
                       (self.t_i[9:]/tau+2*self.beta_i[9:]*(tau-self.gamma_i[9:])))+
                np.sum(self.N_i[9:]*delta**self.d_i[9:]*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*
                       (2*self.beta_i[9:]*(tau-self.gamma_i[9:]))*(self.t_i[9:]/tau+2*self.beta_i[9:]*(tau-self.gamma_i[9:])))+
                np.sum(self.N_i[9:]*delta**self.d_i[9:]*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*
                (-self.t_i[9:]/tau**2+2*self.beta_i[9:])))

    def alpha_r_diff_delta_tau(self,T,rho):
        '''function returning second partial derivative (wrt:tau,delta) of alpha_r, taking Temperature[K] and density [kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return (np.sum(self.N_i[:7]*self.d_i[:7]*delta**(self.d_i[:7]-1)*self.t_i[:7]*tau**(self.t_i[:7]-1))+
                np.sum(self.N_i[7:9]*self.d_i[7:9]*delta**(self.d_i[7:9]-1)*self.t_i[7:9]*tau**(self.t_i[7:9]-1)*np.exp(-delta**self.p_i[7:9]))+
                np.sum(self.N_i[7:9]*delta**self.d_i[7:9]*self.t_i[7:9]*tau**(self.t_i[7:9]-1)*np.exp(-delta**self.p_i[7:9])*(-self.p_i[7:9]*delta**(self.p_i[7:9]-1)))+
                np.sum(self.N_i[9:]*self.d_i[9:]*delta**(self.d_i[9:]-1)*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*
                       (self.t_i[9:]/tau+2*self.beta_i[9:]*(tau-self.gamma_i[9:])))+
                np.sum(self.N_i[9:]*delta**self.d_i[9:]*tau**self.t_i[9:]*np.exp(self.phi_i[9:]*(delta-self.D_i[9:])**2+self.beta_i[9:]*(tau-self.gamma_i[9:])**2)*
                       (2*self.phi_i[9:]*(delta-self.D_i[9:]))*(self.t_i[9:]/tau+2*self.beta_i[9:]*(tau-self.gamma_i[9:]))))

    def alpha_r_diff_delta_delta(self,T,rho):
        pass

    #functions to determine thermodynamic state variables
    #Reference Hydrogen Science & Engineering, Chapter: Thermodynamics of Pressurized Gas Storage(p.601-628), Tietze et al., 2016
    def pressure(self,T,rho):
        '''function returning pressure using Temperature[K] and density[kg/m^3] as input'''
        delta = rho/self.rho_crit
        return rho*self.R*T/self.M*(1+delta*self.alpha_r_diff_delta(T,rho))

    def Z_factor(self,T,rho):
        '''function returning compressability factor using Temperature[K] and density[kg/m^3] as input'''
        delta = rho/self.rho_crit
        return 1+delta*self.alpha_r_diff_delta(T,rho)

    def enthalpy(self,T,rho):
        ''''function returning enthalpy [J/kg] using Temperature[K] and density[kg/m^3] as input'''
        tau, delta = self.T_crit/T, rho/self.rho_crit
        return self.R*T/self.M*(tau*(self.alpha_0_diff_tau(T,rho)+self.alpha_r_diff_tau(T,rho))+delta*self.alpha_r_diff_delta(T,rho)+1)

    def internal_energy(self,T,rho):
        '''function returning internal energy [J/kg] taking temperatue [K] and density [kg/m^3] as input'''
        tau = self.T_crit/T
        return self.R*T/self.M*tau*(self.alpha_0_diff_tau(T,rho)+self.alpha_r_diff_tau(T,rho))

    def entropy(self,T,rho):
        '''function returning entropy [J/(kgK)] using Temperature [K] and density [kg/m^3] as input'''
        tau = self.T_crit/T
        return self.R/self.M*(tau*(self.alpha_0_diff_tau(T,rho)+self.alpha_r_diff_tau(T,rho))-self.alpha_0(T,rho)-self.alpha_r(T,rho))

    def density(self,T,p):
        '''function returning density [kg/m^3] using Temperatur [K] and pressure [Pa] as input'''
        #define expression for root finding solver
        rho_calc = lambda rho:rho*self.R*T/self.M*(1+rho/self.rho_crit*self.alpha_r_diff_delta(T,rho))-p
        # solve for density
        return scpo.fsolve(rho_calc,self.rho_crit/10)[0]

    def temperature(self,rho,u):
        '''function returning temperature [K] using density [kg/m^3] and internal energy [J/kg] as input'''
        #define expression for root finding solver
        T_calc = lambda T: self.R*T/self.M*self.T_crit/T*(self.alpha_0_diff_tau(T,rho)+self.alpha_r_diff_tau(T,rho))-u
        # solve for temperature
        return scpo.fsolve(T_calc,self.T_crit*10)[0]

    def c_v(self,T,rho):
        '''function returning heat capacity [J/(kgK)]for constant volume using temperature [K] and densitiy [kg/m^3] as input'''
        tau = self.T_crit/T
        return self.R/self.M*(-tau**2*(self.alpha_0_diff_tau_tau(T,rho)+self.alpha_r_diff_tau_tau(T,rho)))

    #functions to determine thermodynamic process variables for compression and expansion processes
    #Reference Hydrogen Science & Engineering, Chapter: Thermodynamics of Pressurized Gas Storage(p.601-628), Tietze et al., 2016

    def isothermal(self,T,p_bounds):
        '''function returning specific work [J/kg] and specific waste heat [J/kg]
        from caloric equation of state using enthalpy and entropy for an isothermal compression/expansion,
        taking Temperature [K] and in- and outlet pressures [Pa] as list like object as input'''
        #get densities for inlet and outlet state
        rho = [self.density(T,p) for p in p_bounds]
        q_spec = T*(self.entropy(T,rho[1])-self.entropy(T,rho[0]))
        w_spec = self.enthalpy(T,rho[1])-self.enthalpy(T,rho[0])-q_spec
        return (w_spec,q_spec)

    def isentropic(self,T_in,p_bounds):
        '''function returning specific work [J/kg] and outlet temperature [K]
        from caloric equation of state using enthalpy and entropy for an isentropic compression/expansion,
        taking inlet Temperature [K] and in- and outlet pressures [Pa] as list like object as input'''
        # compute inlet conditions
        rho_in = self.density(T_in,p_bounds[0])
        h_in = self.enthalpy(T_in,rho_in)
        s = self.entropy(T_in,rho_in)

        #define iterative expression for the outlet temperature assuming constant entropy
        def T_calc(T_out,s,p):
            rho = self.density(T_out,p)
            tau = self.T_crit/T_out
            return self.R/self.M*(tau*(self.alpha_0_diff_tau(T_out,rho)+self.alpha_r_diff_tau(T_out,rho))-self.alpha_0(T_out,rho)-self.alpha_r(T_out,rho))-s

        #solve for T and compute outlet conditions
        T_out = scpo.fsolve(T_calc,400,args=(s,p_bounds[1]))[0]
        h_out = self.enthalpy(T_out,self.density(T_out,p_bounds[1]))
        w_spec = h_out - h_in
        return (w_spec,T_out)


    def isenthalpic(self,T_in,p_bounds):
        '''function returning Outlet temperature [K] from caloric equation of state using enthalpy and entropy for an isenthalpic compression/expansion,
        taking inlet Temperature [K] and in- and outlet pressures [Pa] as list like object as input'''
        # compute inlet conditions
        rho_in = self.density(T_in,p_bounds[0])
        h = self.enthalpy(T_in,rho_in)

        # define iterative expression for outlet Temperature assuming constant enthalpy
        def T_calc(T_out,h,p):
            rho = self.density(T_out,p)
            tau, delta = self.T_crit/T_out,rho/self.rho_crit
            return self.R*T_out/self.M*(tau*(self.alpha_0_diff_tau(T_out,rho)+self.alpha_r_diff_tau(T_out,rho))+delta*self.alpha_r_diff_delta(T_out,rho)+1)-h

        # solve for T
        return scpo.fsolve(T_calc,300,args=(h,p_bounds[1]))[0]
