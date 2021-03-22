'''
calculation: power consumption and control
    -> optimal power/ operation point
        -> suitable power gradient
        ->
    -> power consumption BoP
        -> power consumption pumps (feedwater, (coolant?, compressor?)
        -> power characteristics (efficiencies): rectifier, el-motor, fan?
        -> power consumption gas treatment (dryer)

'''
print(__name__ + ' imported...')


### calc optimal operation point
def objective_potp(T, i, P, p, cnc):
        '''
        objective function for potp
        '''
        pol     = m_plr.polar_clc(T, i, p, cnc, aux3in) #returns i in A/m² ,U_cell in V, P in W/m² /// polarc ehem. polar4

        P_diff  = P - (pol[1] * pv.N * pv.A_cell) # edit: 2019-06-13
        objective.pout = pol

        return abs(P_diff)

def op_opt(T, i, ):
    # TODO: what, if u > u_max?

    # initial value for current density
    if i < 0.5:
        i +=1
    x0 = [i]         #,0] # initial value for i_in

    # bounds for optimization
    bnds = [0 , pv.i_max]
    sol = minimize (objective,x0,args=(P_in, T_in, p_in, cnc_in, pv,av),method='SLSQP',bounds=bnds)#,constraints=cons)

    return sol.x ,sol.success, objective.pout


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


### Balance of plant
def pwr_BoP(n_H2, m_c_in, M_ely_in, p, P_heat):
    ''' power consumption of plant/ periphery'''


    V0 = pv.V0_ely
    P_ely   = pow_pump_V(Vely, V0, p_C, pv0, pv)
    P_pmp   = pow_pump_m(m_c_in, m_c0, p_A, pv0,pv,av)    # power input feed/cooling pump     // in W

    P_gt    = pow_gas_trt(n_H2, pv0) # power input gas treatment         // in W

    m_c0    = pv.m_cmax

    #eff_pe  = ip_rect()     # rectifier efficiency
    #P_pe    = P_st * (1-eff_pe) # power input/loss power electronics
#print('P_gt:',P_gt,'P_pmp:',P_pmp, 'P_ely:',P_ely)
    P_bop = P_gt + P_pmp + P_ely # // in W

    return P_bop

def pwr_pmp(m_flow, p_diff):
    '''
    calc power of electric pump

    '''
    P_clc = dV * dp
    P_act = P_clc / eff_pmp(P_clc, P0)
    return

def pwr_gasdryer(obj, pec, n_H2):        # Tjarks
    '''
    gas treatment -> drying via tsa

    '''
    Xtar = getattr(obj.pop,'purity_target_H2', None)
    ## ->>> n_H2 in mol/s !
    if Xtar:
        #X_tar       = 0.0005 #? # target output purity (max. water content)
        X_in        = pec.tsa_Xin #0.08      # fixed humidity of hydrogen after condenser // in molH2O/molH2
        X_out2      = pec.tsa_Xdes #0.9       # fixed humidity of hydrogen after desorbtion // in molH2O/molH2
        H_ads       = pec.tsa_Hads #48.6*10**(3)      # adsorption enthapy // in J/mol
        deltaT      = pec.tsa_dT #40        # Temp diff heater from Tjarks: 40K
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
    else:
        P_gt = 0
    #print('n_H2:',n_H2,'P_gt:',P_gt)
    return P_gt # in W
