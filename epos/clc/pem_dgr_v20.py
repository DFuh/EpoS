'''
PEM
calculation: voltage increase trough degradation
'''
print(__name__ + ' imported...')

import numpy as np

def voltage_increase_lin(obj, ):

    ### absolute voltage increase
    dU_dgr_abs = obj.av.t_op/3600 * obj.pec.fctr_vlr # h * V/h

    ### incremental voltage increase
    dU_dgr_incr = obj.av.t_diff/3600 * obj.pec.fctr_vlr # h * V/h

    return dU_dgr_incr, dU_dgr_abs


def voltage_increase_lfun(obj, pec, T, i):
    '''
    ...at this state just a wrapper for clc_vlr_tot()
    '''
    (dU_abs,
    dU_act,
    dU_rev,
    dU_irr) = clc_vlr_tot(obj, pec, T ,i,
                            obj.av.dU_dgr_abs, obj.av.t_diff, obj.av.t_uni)

    return dU_act, dU_abs


def clc_vlr_tot(obj, pec, T ,i, dU_abs_pre, dt_in, t_uninterrupt, ):
    '''
    adopted from  EpoS_apply_dgr_fun_v01 (in mod_3/dgr )

    Parameters
    ----------
    obj : class
        instance of simulation
    pec : namedtuple
        container for tec-specific parameters
    T : float
        cell (Stack) Temperature // in K
    i : float
        current density // in A/m2
    dU_abs_pre : float
        previous value of voltage increase // in V
    dt_in : float
        timedelta // in s
    t_uninterrupt : float
        period of time, where cell voltage lay above treshold // in s

    Returns
    -------
    dU_abs : float
        DESCRIPTION.
    dU_act : float
        DESCRIPTION.
    dU_rev : float
        DESCRIPTION.
    dU_irr : float
        DESCRIPTION.

    '''
    # a = obj.pec.lim_hi_vlr_rev
    # b = obj.pec.tfrc_vlr_rev # 1/600
    vlr_rev = efun_incr(t_uninterrupt,
                            obj.pec.lim_hi_vlr_rev,
                            obj.pec.tfrc_vlr_rev)

    vlr_irr = 1 # f(i_ref)

    # dU_abs_pre = 0
    # t_uninterrupt = 0 # time increment of uninterrupted operation

    f_i = 1 # Factor accounting for influence of curent density on vlr
    f_T = fctr_vlr_T(T, obj.pec.vlr_min_Tfun_vlr, obj.pec.vlr_max_Tfun_vlr,
                        obj.pec.T_min_Tfun_vlr, obj.pec.T_max_Tfun_vlr)
    # f_T = fctr_vlr_T(T) # Factor accounting for influence of temperature on vlr

    f_i_ref = 1 # Factor accounting for impact of vlr at different current densities
    dt = dt_in/3600 # Time increment // in h

    dU_rev = vlr_rev * dt * f_T * f_i               # Reversible voltage increase in current timestep
    dU_irr = vlr_irr * dt * f_T * f_i               # Irreversible voltage increase in current timestep
    dU_act = (dU_abs_pre + dU_rev + dU_irr) * f_i_ref    # Total voltage increase in curent timestep
    dU_abs = dU_abs_pre + dU_irr                        # Absolute voltage increase at current state of plant

    return dU_abs, dU_act, dU_rev, dU_irr

def fctr_vlr_T(T, vlr_min, vlr_max, T_min, T_max):
    # vlr_min = 0
    # vlr_max = 50
    # T_min = 300
    # T_max = 353
    f = (vlr_max-vlr_min)/(T_max-T_min) * (T - T_min) + vlr_min
    f = f if T> T_min else vlr_min
    f = f if T< T_max else vlr_max
    return f/vlr_max

def efun_incr(x, a,b):
    '''
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    '''
    a - max output value
    b - 'gain', how fast, does y reach a?
        -> possible approach: -b = -4...-9/t0 (ln(0.01) = -4.605 || ln(0.0001) = -9.2103 )
    '''
    return a * (1-np.exp(-b*x))
