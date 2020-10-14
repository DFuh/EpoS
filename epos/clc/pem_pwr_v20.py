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
import scipy.optimize as scpt

print(__name__ + ' imported...')


### calc optimal operation point
def objective_popt(i, obj, pec, P, T, p, pp, ifu, ini):
        '''
        objective function for potp
        '''
        pol     = m_plr.polar_clc(obj, pec, T, i, pp=pp) #returns i in A/m² ,U_cell in V, P in W/m² /// polarc ehem. polar4

        P_diff  = P - (pol[1] * pv.N * pv.A_cell) # edit: 2019-06-13
        objective_popt.out = pol

        return abs(P_diff)

def objective_iopt(i, obj, pec, u_tar, T, p, pp, ifu, ini):
        '''
        objective function for potp
        '''
        if ifu:
            pol     = ifu(obj, pec, T, i, p, pp=pp, ini=ini) #returns
        else:
            pol     = m_plr.polar_clc(T, i, p, pp=pp) #returns i in A/m² ,U_cell in V, P in W/m² /// polarc ehem. polar4

        u_diff  = u_tar - pol[0]
        objective_iopt.out = pol

        return abs(u_diff)


def op_opt(obj, pec, T_in, i, i_max, p_in, pp_in, P_in=None, u_mx=None, ifun=None, ini=False):
    # initial value for current density
    if i < 0.5:
        i +=1
    x0 = [i]         #,0] # initial value for i_in

    # TODO: what, if u > u_max?
    # TODO: implement maximum cell_voltage -> shortcut, if reached
    # bounds for optimization
    bnds = [(0 , i_max)]

    if u_mx is not None:
        tar_val = u_mx
        obj_fun = objective_iopt
    else:
        tar_val = P_in
        obj_fun= objective_popt

    sol = scpt.minimize (obj_fun,x0,args=(obj, pec, tar_val, T_in, p_in, pp_in, ifun, ini),method='SLSQP',bounds=bnds)#,constraints=cons)

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
def BoP_power(n_H2, m_c_in, Vely, p_in, ):
    ''' power consumption of plant/ periphery'''

    ### gas treatment (drying)
    P_gt    = pow_gas_trt(n_H2, ) # power input gas treatment         // in W

    m_c0    = pv.m_cmax
    P_pmp   = pow_pump_m(m_c_in, m_c0, p_A, )    # power input feed/cooling pump     // in W
    #eff_pe  = ip_rect()     # rectifier efficiency
    #P_pe    = P_st * (1-eff_pe) # power input/loss power electronics

    P_bop = P_gt + P_pmp + P_ely # // in W

    return P_bop
