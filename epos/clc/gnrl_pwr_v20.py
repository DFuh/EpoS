'''
general power calculations
'''
import numpy as np
import scipy.optimize as scpt

### calc optimal operation point

def cntrl_pow_clc(obj, pec, T, i, p, pp, P_avail):
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
    P_avail: float
        Available Power // in kW

    Returns
    -------
    P: float
        Actual stack-power // in kW
    P_rect: float
        Actual rectifier power // in kW
    '''

    ### if P_availabel in critical range (close to nominal power), calc i_max based on u_max
    i_lim = obj.pcll.current_density_max

    if P_avail / obj.pplnt.power_of_plant_max >= 0.8:
        ppout = op_opt(obj, pec, T, i, i_lim, p, pp, u_mx=obj.pcll.voltage_max)
        i_lim = ppout[0]
        print('i_lim: ', i_lim)
    ### if P_availabel not critical, calc i_opt directly

    pout = op_opt(obj, pec, T, i, i_lim, p, pp, P_avail)

    i_o = pout[0]
    u_o = pout[-1][-1]
    P_rect = pout[-1][3]
    P = i_o * u_o * obj.pplnt.number_of_cells_in_plant_act * obj.pcll.active_cell_area/1e3
    return P, P_rect, i_o, u_o


def objective_popt(i, obj, pec, P, T, p, pp, ifu, ini):
        '''
        objective function for popt
        '''

        pol     = obj.clc_m.plr.voltage_cell(obj, pec, T, i, p, pp=pp)
        # NO LONGER TRUE: returns i in A/m² ,U_cell in V, P in W/m² /// polarc ehem. polar4

        # TODO: break-condition (for low values of u_max) )in order to improve performance?
        #if pol[-1] > umx:
    #        u = umx
        #else:
        u = pol[-1]
        P_EL = (u *i* obj.pplnt.number_of_cells_in_plant_act * obj.pcll.active_cell_area)/1000 # // in kW
        P_N = obj.pplnt.power_of_plant_nominal
        eta_rect = efficiency_rectifier(obj, pec, P_EL / P_N)
        P_rect = P_EL * (1-eta_rect)
        #P_diff  = P - (pol[1] * pv.N * pv.A_cell) # edit: 2019-06-13
        P_diff = P - P_EL*(1 + (1-eta_rect))
        objective_popt.out = (P, P_EL, P_N, P_rect, eta_rect, P_diff, u)

        return abs(P_diff)


def objective_iopt(i, obj, pec, u_tar, T, p, pp, ifu, ini):
        '''
        objective function for potp
        '''
        if ifu:
            pol     = ifu(obj, pec, T, i, p, pp=pp, ini=ini) #returns
        else:
            pol     = obj.clc_m.plr.voltage_cell(obj, pec, T, i, p, pp=pp) #returns (U_ca, U_an, U_cell in V, /// polarc ehem. polar4

        #print('u_tar: ', u_tar)
        #print('u: ', pol[-1])
        u_diff  = u_tar - pol[-1]
        objective_iopt.out = pol

        return abs(u_diff*1e3)


def op_opt(obj, pec, T_in, i, i_max, p, pp_in, P_in=None, u_mx=None, ifun=None, ini=False):
    # initial value for current density
    if i < 0.5:
        i +=1
    x0 = [i]         #,0] # initial value for i_in

    # TODO: what, if u > u_max?
    # TODO: implement maximum cell_voltage -> shortcut, if reached
    # TODO: check tolerance
    # bounds for optimization
    bnds = [(0 , i_max)]

    if u_mx is not None:
        tar_val = u_mx
        obj_fun = objective_iopt
        #umx = None
    else:
        tar_val = P_in
        obj_fun= objective_popt
        #if obj.pcll.voltage_max:
        #    umx = obj.pcll.voltage_max
        #else:
        #    umx=1e3
    #print('Bounds: ', bnds)
    sol = scpt.minimize (obj_fun,x0,
                            args=(obj, pec, tar_val, T_in, p, pp_in, ifun, ini),
                            method='SLSQP',bounds=bnds, tol=1e-2)#,constraints=cons)
    #print('Sol: ', sol)
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

def efficiency_rectifier(obj, pec, x_in):
    '''
    rectifier efficiency
    '''
    #TOD: a,b,c = obj.fv.eff_rectifier_fit_vals
    a, b, c = 0.0578373, 0.0417102, 0.811693

    if x_in < 0.05:
        yo = 0.88
    else:
        yo = eff_efun(x_in,a,b,c)
    return yo

#original aus Rodriguez bzw. Tjarks abgelesen
def eff_efun(x,a,b,c):
    ''' fitted e-funct. for power electronics/ motor efficiency '''
    out = a * np.exp(1-(b/ x)) + c
    return out
