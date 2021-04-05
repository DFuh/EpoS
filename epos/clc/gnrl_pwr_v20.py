'''
general power calculations
'''
import numpy as np
import scipy.optimize as scpt

### calc optimal operation point

def plnt_ctrl():
    '''
    Control important operational functions/ values

    - on/off

    - massflow ely

    - load restrictions

    -
    '''

    return

def cntrl_pow_clc(obj, pec, T, i, p, pp, P_in, P_prev, u_prev, dt):
    '''
    Find optimal operation point of Stack based on given/ available power
    using scipy.optimize.minimize()
    if P_avail in critical range, calc limiting current density (based on maximum voltage)

    Parameters
    ----------
    obj: object
        Simu object
    pec: Namedtuple
        Parameters for electrochemical components
    T: float
        Temperature of stack // in K
    i: float
        Initial value for current density // in A/m²
    p: Namedtuple
        Electrode pressures // in Pa
    pp: tuple
        Partialpressure of species // in Pa
    P_in: float
        Available Power for one Stack // in kW
    P_prev: float
        Previous value of Stack power (single) // in kW
    u_prev: float
        Previous value of cell voltage // in V
    dt: float
        timeincrement of calculation step // in s ???

    Returns
    -------
    P: float
        Actual stack-power // in kW
    P_rect: float
        Actual rectifier power // in kW
    '''


    i_lim = [(0,obj.pcll.current_density_max)]
    P_avail = P_in / obj.pplnt.number_of_stacks_act

    #TODO: how to skip steps? (in case dPdt << dPdtmax)
    # -> if no dudt given -> static value
    ### if P_availabel in critical range (close to nominal power), calc i_max based on u_max


    # Maximum gradient of voltage
    if hasattr(obj.av, 'dudt_p'):
        umax = u_prev+obj.av.dudt_p*dt
    else:
        umax = obj.pcll.voltage_max
    if hasattr(obj.av, 'dudt_n'):
        umin = u_prev-obj.av.dudt_n*dt
    else:
        umin = None

    ### Power limit

    # Maximum stack power |||
    '''
    may partly be moved to ini auxvals, calculation only once necessary !
    '''
    if ((obj.pop.power_fraction_max >1)
        & (obj.av.t_ol < obj.pop.overload_time_max)
        & obj.av.t_nom >= obj.pop.overload_recovery_time_min):

        P_N = obj.pplnt.power_of_stack_nominal * obj.pop.power_fraction_max
        obj.av.t_ol += dt
        obj.av.t_nom = 0
    else:
        P_N = obj.pplnt.power_of_stack_nominal
        obj.av.t_ol = 0
        obj.av.t_nom += dt

    # Maximum gradients of power
    # Redundant code !
    '''
    -> P_N has to be P_N_rectifier !
    '''
    P_N_rect = obj.bop.power_rectifier_nominal
    if hasattr(obj.av, 'dPdt_p') and (P_in/P_prev >1):
        print('obj.av.dPdt_p: ', obj.av.dPdt_p)
        eta_rect = efficiency_rectifier(obj, (P_prev + obj.av.dPdt_p*dt) / P_N_rect)
        P_avail_max = (P_prev + obj.av.dPdt_p*dt) * (1+(1-eta_rect))
        if P_avail_max < P_in:
            P_avail = P_avail_max

    if hasattr(obj.av, 'dPdt_n') and (P_in/P_prev <1):
        print('obj.av.dPdt_n: ',obj.av.dPdt_n)
        eta_rect = efficiency_rectifier(obj, (P_prev - obj.av.dPdt_n*dt) / P_N_rect)
        P_avail_min = (P_prev - obj.av.dPdt_n*dt)* (1+(1-eta_rect))
        if P_avail_min > P_in:
            P_avail = P_avail_min

    ### Minimum stack power

    print('P_avail: ', P_avail)
    ### Calc limits for current density based on voltage and power limits
    ppout0 = op_opt(obj, pec, T, i, i_lim, p, pp, u_mx=umax)
    if umin:
        ppout1 = op_opt(obj, pec, T, i, i_lim, p, pp, u_mx=umin)
    else:
        ppout1 = [[0,],]

    i_lim = [(ppout1[0][0], ppout0[0][0])]
    print('i_lim: ', i_lim)


    ### Calc optimal power of stack
    pout = op_opt(obj, pec, T, i, i_lim, p, pp, P_avail)

    i_o = pout[0]
    u_o = pout[-1][-1]
    P_rect = pout[-1][3]
    #P = i_o * u_o * obj.pplnt.number_of_cells_in_plant_act * obj.pcll.active_cell_area/1e3
    P = i_o * u_o * obj.pcll.active_cell_area/1e3 * obj.pplnt.number_of_cells_in_stack_act

    return P, P_rect, i_o, u_o


def objective_popt(i, obj, pec, T, p, pp, P, P_N):
    '''
    objective function for popt within operational control
    '''

    pol     = obj.clc_m.plr.voltage_cell(obj, pec, T, i, p, pp=pp)
    # NO LONGER TRUE: returns i in A/m² ,U_cell in V, P in W/m² /// polarc ehem. polar4

    # TODO: break-condition (for low values of u_max) )in order to improve performance?
    u = pol[-1]
    #P_EL = (u *i* obj.pplnt.number_of_cells_in_plant_act * obj.pcll.active_cell_area)/1000 # // in kW
    A_cell = obj.pcll.active_cell_area
    n_clls = obj.pplnt.number_of_cells_in_stack_act
    P_EL = (u *i* n_clls * A_cell) / 1e3 # // in kW
    #P_N = obj.pplnt.power_of_plant_nominal
    eta_rect = efficiency_rectifier(obj, P_EL / P_N)
    P_rect = P_EL * (1-eta_rect)
    #P_diff  = P - (pol[1] * pv.N * pv.A_cell) # edit: 2019-06-13
    P_diff = P - P_EL*(1 + (1-eta_rect))
    objective_popt.out = (P, P_EL, P_N, P_rect, eta_rect, P_diff, u)

    return abs(P_diff)

def objective_popt_bsc(i, obj, pec, T, p, pp, P, A_cell):
    '''
    objective function for popt in pwr_vls calculation
    '''

    pol     = obj.clc_m.plr.voltage_cell(obj, pec, T, i, p, pp=pp, )

    u = pol[-1]
    #P_EL = (u *i* obj.pplnt.number_of_cells_in_plant_act * obj.pcll.active_cell_area)/1000 # // in kW
    P_EL = i * pol[-1] *1e-3# // in kW

    P_diff = P - P_EL
    objective_popt_bsc.out = (P, P_EL, u)

    return abs(P_diff*1e3)


def objective_uopt(i, obj, pec, T, p, pp, u_tar, P_N):
    '''
    objective function for iotp within operational control and pwr_val-calculation
    '''

    pol     = obj.clc_m.plr.voltage_cell(obj, pec, T, i, p, pp=pp) #returns (U_ca, U_an, U_cell in V, /// polarc ehem. polar4
    print(f'(uopt) i: {i}, pol (uopt): {pol}')
    u_diff  = u_tar - pol[-1]
    objective_uopt.out = pol

    return abs(u_diff*1e3)

def objective_uA_opt(i, obj, pec, T, p, pp, u_tar, P_N):
    '''
    objective function for iotp within operational control and pwr_val-calculation
    '''

    pol     = obj.clc_m.plr.voltage_cell(obj, pec, T, i, p, pp=pp) #returns (U_ca, U_an, U_cell in V, /// polarc ehem. polar4

    u_diff  = u_tar - pol[-1]
    objective_uA_opt.out = pol

    return abs(u_diff*1e3)


def op_opt(obj, pec, T_in, i, i_max, p, pp_in, P_in=None, u_mx=None):
    '''
    optimization for operational control values
    '''
    # initial value for current density
    if i < 0.5:
        i +=1
    x0 = [i]         #,0] # initial value for i_in

    # TODO: what, if u > u_max?
    # TODO: implement maximum cell_voltage -> shortcut, if reached
    # TODO: check tolerance
    # bounds for optimization
    #bnds = i_max #[(0 , i_max)]

    #if not A_cell:
    #A_cell = None #obj.pcll.active_cell_area
    P_N = None # ????
    #if not n_clls:
    #    n_clls = obj.number_of_cells_in_plant_act #???


    if u_mx is not None:
        tar_val = u_mx # // in V
        obj_fun = objective_uopt
        P_N = None
        #umx = None
    else:
        tar_val = P_in # // in kW
        P_N = obj.bop.power_rectifier_nominal #pplnt.power_of_stack_nominal
        obj_fun= objective_popt

    sol = scpt.minimize (obj_fun, x0,
                            args=(obj, pec, T_in, p, pp_in, tar_val, P_N),
                            method='SLSQP',bounds=i_max, tol=1e-2)#,constraints=cons)
    #print('Sol: ', sol)
    print(f'Sol: ({obj_fun.__name__}): x0: {x0} x: {sol.x} success: {sol.success} out:{obj_fun.out}')
    return sol.x ,sol.success, obj_fun.out


def bsc_opt(obj, pec, T_in, i, i_max, p, pp_in, P_in=None, u_mx=None, A_cell=None):
        '''
        optimization for basic calculation of plant/ cell properties (pwr_vls)
        '''
        # initial value for current density
        if i < 0.5:
            i +=1
        x0 = [i]         #,0] # initial value for i_in

        #if not A_cell:
        #A_cell = obj.pcll.active_cell_area
        P_N = None # ????

        if u_mx is not None:
            tar_val = u_mx # // in V
            obj_fun = objective_uopt
            args_in = (obj, pec, T_in, p, pp_in, tar_val, None)
            #umx = None
        else:
            tar_val = P_in # // in W/m²
            obj_fun = objective_popt_bsc
            args_in = (obj, pec, T_in, p, pp_in, tar_val, A_cell)

        sol = scpt.minimize (obj_fun,x0, args= args_in,
                                method='SLSQP',bounds=i_max, tol=1e-2)

        print('Sol: ',sol.x ,sol.success, obj_fun.out)
        return sol.x ,sol.success, obj_fun.out

### calc and limit power gradient

def grad_pwr(P_new, P_old, dtime, ):
    '''
    calc power and limit power gradient
    '''

    # dtime in s!
    P_N         = pv.P_N
    pct_time    = av.pct_t      # time count for operation at maximum power

    dPdt_max    = pv.dp_max * P_N      # set maximum power gradient +++REDUNDANT+++
    #P_min       = pv.fr_pmin * P_N     # set minimum power level +++RED+++

    # check, if maximum over-load-time is reached
    if pct_time > pv.t_pmax: #edit:2019-06-13
        P_max   = P_N   # max power is nominal power
    else:
        P_max   = pv.fr_pmax * P_N # maximal power is max overload power

    dPr     = P_old / P_new    # calc fraction of old/ new power value

    #dPdt    = (P_new-P_old) / dtime  # possible power gradient from input

    #if np.isnan(dPr) or (P_new < 0): # check low values of Pnew /// cannot happen
    #    P_out   = 0

    #elif dPr ==1:
    if dPr ==1: # New power value equal to old
        P_out = P_new

    elif dPr > 1: # no limit of neg. power gradient (?)
        P_out = P_new

    elif (dPr <1) or (P_old == 0):
        dPdt_act = (P_max - P_old) / dtime # possible power gradient
        dPdt_new = (P_new - P_old) / dtime # powr gradient due to new input power

        if dPdt_act <= dPdt_new <= dPdt_max: # gradient within limit but out of power-range
            #if dPdt_act < dPdt_max:#
            dPdt = dPdt_act
            #else:
                #dPdt = dPdt_max
        elif dPdt_new <= dPdt_act <= dPdt_max: # gradient within all limits
            #if dPdt_new < dPdt_max:
            dPdt = dPdt_new
            #else:
            #    dPdt = dPdt_max
        else:                                   # gradient out of all limits
            dPdt = dPdt_max

        P_out = P_old + (dPdt*dtime)

    return P_out, pct_time


### ===========================================================================
            # Power electronics
### ===========================================================================
def efficiency_cmp(obj, x_in, fun, fitvals, lolim=[0,0], hilim=[1,1] ):

    if x_in < lolim[0]:
        y = lb[1]
    elif x_in > hilim[0]:
        y = hilim[1]
    else:
        y = fun(x_in, fitvals)
    return y


def efficiency_rectifier(obj, x_in):
    '''
    rectifier efficiency
    '''
    #TOD: a,b,c = obj.fv.eff_rectifier_fit_vals
    #a, b, c = 0.0578373, 0.0417102, 0.811693
    a,b,c = obj.bop.fitvals_efficieny_rectifier
    lolim = obj.bop.xlim_efficieny_rectifier
    hilim = obj.bop.ylim_efficieny_rectifier
    #xmin= 0.05
    #ymin = 0.88

    #eff_mot fitvals:
    #a, b, c = 0.25236, 0.0437787, 0.194493
    #xmin = 0.01
    #ymin = 0.2

    if x_in < lolim[0]:
        y = lolim[1]
    elif x_in > hilim[0]:
        y = hilim[1]
    else:
        y = eff_efun(x_in, a,b,c)
    return y

def efficiency_emotor(obj, x_in):
    '''
    rectifier efficiency
    '''
    #TOD: a,b,c = obj.fv.eff_rectifier_fit_vals
    #a, b, c = 0.0578373, 0.0417102, 0.811693
    a,b,c = obj.bop.fitvals_efficieny_motor
    lolim = obj.bop.xlim_efficieny_motor
    hilim = obj.bop.ylim_efficieny_motor

    #eff_mot fitvals:
    #a, b, c = 0.25236, 0.0437787, 0.194493
    #xmin = 0.01
    #ymin = 0.2

    if x_in < lolim[0]:
        y = lolim[1]
    elif x_in > hilim[0]:
        y = hilim[1]
    else:
        y = eff_efun(x_in, a,b,c)
    return y

#original aus Rodriguez bzw. Tjarks abgelesen
def eff_efun(x,a,b,c):
    ''' fitted e-funct. for power electronics/ motor efficiency '''
    out = a * np.exp(1-(b/ x)) + c
    return out

### ===========================================================================
            # Auxilliary Power (BoP)
### ===========================================================================

def clc_pwr_bop(obj, m_ely, m_coolant, n_H2, P_heat):
    '''
    Calc power of balance of plant
    --------------
    P_pumps (Coolant, feedwater/electrolyte, )
        - coolant:
            T_in = 15 °C
            T_out = ?
            m_c = PID-output
            dp = dp_hex
        - feedwater/electrolyte:
            T_in = T_out_hex
            m_ely = PID-output ? (PEM -> Anode, AEL -> both compartments)
            dp = dp_cell_comaprtment

    P_gasdryer

    Returns
    --------
    P_aux (on plant level)
    '''
    P_di = 0 # Power consumption of DI-Unit ???

    Vdot_coolant = (m_coolant / obj.av.rho_H2O) #* obj.plnt.number_of_stacks_act
    Vdot_ely = (m_ely / obj.av.rho_ely)#* obj.plnt.number_of_stacks_act
    delta_p_coolant = obj.bop.dp_coolant_cycle
    delta_p_ely = obj.bop.dp_ely_cycle
    #print(f'V_ely: {Vdot_ely}|| V_clnt: {Vdot_coolant}')
    #print(f'eta_opt: {obj.bop.eta_opt_pump}|| V_clnt: {obj.bop.power_pump_coolant_nominal, obj.bop.power_pump_ely_nominal}')
    #print(f'dp_ely: {delta_p_ely}|| dp_clnt: {delta_p_coolant}')
    #print(f'Vely_nom: {obj.bop.volumetricflow_ely_nominal} || V_clnt_nom: {obj.bop.volumetricflow_coolant_nominal}')
    V0 = obj.bop.volumetricflow_coolant_nominal * obj.pplnt.number_of_cells_in_stack_act
    P_pc = pwr_pmp(obj, Vdot_coolant, delta_p_coolant, V0,
                    obj.bop.eta_opt_pump, obj.bop.power_pump_coolant_nominal) # Coolant pump
    P_pely = pwr_pmp(obj, Vdot_ely, delta_p_ely, V0,
                    obj.bop.eta_opt_pump, obj.bop.power_pump_ely_nominal) # Electrolyte/ feedwater Pump

    P_gt = 0#pwr_gasdryer(obj, obj.pec, n_H2)#*obj.plnt.number_of_stacks_act)

    print('P_aux_cmpnts: ', P_di, P_gt, P_pc, P_pely)

    return (P_di + P_gt + P_pc + P_pely + P_heat)


def pwr_pmp(obj, v_dot, delta_p, v0, eta_opt, P_N_pmp):
    '''
    calc power of electric pump
    v_dot: float
        volumetric flowrate in m³/s
    delta_p: float
        pressure difference (losses) of fluid-system in Pa
    '''
    if (P_N_pmp >0) & (v_dot >0):
        P_clc = v_dot * delta_p
        #if np.isinf(P_clc):
        #print('P_clc: (v_dot, dp): ', v_dot, delta_p)
        print('eff_pmp: ', eff_pmp(v_dot, v0, eta_opt))
        P_act = P_clc / eff_pmp(v_dot, v0, eta_opt)
        #if np.isinf(P_clc):
        #print('P_act: (v_dot, v0, eta_opt): ', v_dot, v0, eta_opt)
        P_e = P_act / efficiency_emotor(obj, P_act/P_N_pmp)
        #if np.isinf(P_clc):
        #print('P_e: (P_act, P_N_pmp): ', P_act, P_N_pmp)
    else:
        P_e = 0
    print(f'pwr_pmp: P_clc: {P_clc}, P_act={P_act}, P_e={P_e}')
    return P_e*1e-3 #// in kW

def eff_pmp(Q_P, Q_Popt, eta_Popt):
    '''
    Efficiency of fluid-pump based on operating point

    eta_P(Q_P)
    Diss SchuetzoldS.45
    '''
    #(Q_Popt, eta_Popt) = prs
    #print(f'eta_pmp: Q_P={Q_P}, Q_Popt={Q_Popt}, eta_Popt={eta_Popt}')
    eta = eta_Popt * ( (2*Q_P / Q_Popt) - (Q_P**2 / Q_Popt**2)  )
    if eta< 0.1:
        eta=0.1
    #print('eta_ ', eta)
    return eta

def pwr_gasdryer(obj, pec, n_H2):        # Tjarks
    '''
    gas treatment -> drying via tsa

    '''
    X_tar = getattr(obj.pop,'purity_target_H2', 0.0005)#None)
    #if not X_tar:
    #    X_tar = 0.0005

    ## ->>> n_H2 in mol/s !
    #if Xtar:
    #X_tar       = 0.0005 #? # target output purity (max. water content)
    X_in        = obj.bop.tsa_Xin #0.08      # fixed humidity of hydrogen after condenser // in molH2O/molH2
    X_out2      = obj.bop.tsa_Xdes #0.9       # fixed humidity of hydrogen after desorbtion // in molH2O/molH2
    H_ads       = obj.bop.tsa_Hads #48.6*10**(3)      # adsorption enthapy // in J/mol
    deltaT      = obj.bop.tsa_dT #40        # Temp diff heater from Tjarks: 40K
    n_H2in2     = n_H2 * X_in / (X_out2 * X_in)
    n_H2in1     = (n_H2 + n_H2in2)
    n_H2Oin1    = n_H2in1 * X_in
    n_H2out     = n_H2in1 - n_H2in2
    n_H2Oout    = X_tar * n_H2out
    delta_H2O   = n_H2Oin1 - n_H2Oout
    n_H2in2     = 1 / X_out2 * (delta_H2O + n_H2Oin1)

    P_ht        = pec.cp_H2 * pec.M_H2 * n_H2in2 * deltaT # J/kgK * kg/mol * mol/s * K = J/s
    Q_des       = H_ads * delta_H2O     # in J/s
    P_gt        = P_ht + Q_des
    #else:
    #    P_gt = 0
    #print('n_H2:',n_H2,'P_gt:',P_gt)
    return P_gt*1e-3 # in kW
