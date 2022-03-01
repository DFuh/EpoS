'''
PEM
calculation: voltage increase trough degradation
'''
print(__name__ + ' imported...')

import numpy as np

def voltage_increase_lin(obj, pec, T, i):

    ### absolute voltage increase
    dU_dgr_act = obj.av.t_op/3600 * obj.pec.fctr_vlr # h * V/h

    ### incremental voltage increase
    # dU_dgr_incr = obj.av.t_diff/3600 * obj.pec.fctr_vlr # h * V/h

    return dU_dgr_act, 0 #dU_dgr_abs


def voltage_increase_lfun(obj, pec, T, i):
    '''
    ...at this state just a wrapper for clc_vlr_tot()
    '''
    (dU_abs,
    dU_act,
    dU_rev,
    dU_irr) = clc_vlr_tot(obj, pec, T ,i,
                            obj.av.dU_dgr_abs, obj.av.t_diff, obj.av.t_uni)
    # print('Res (dgr_abs_in volt incr lfun: ', obj.av.dU_dgr_abs)
    # print('Res (dgr_abs_out volt incr lfun: ', dU_abs)
    return dU_act, dU_abs

################################################################################

def clc_vlr_tot(obj, pec, T ,i, dU_abs_pre, dt_in, t_uninterrupt, print_vals=False):
    '''


    Parameters
    ----------
    T : TYPE
        DESCRIPTION.
    i : TYPE
        DESCRIPTION.
    dU_abs_pre : TYPE
        DESCRIPTION.
    dt_in : TYPE
        DESCRIPTION.
    t_uninterrupt : TYPE
        DESCRIPTION.

    Returns
    -------
    dU_abs : TYPE
        DESCRIPTION.
    dU_act : TYPE
        DESCRIPTION.
    dU_rev : TYPE
        DESCRIPTION.
    dU_irr : TYPE
        DESCRIPTION.

    '''
    i = i/1e4
    T = T-273
    # a = 1
    # b = 1/600
    #vlr_rev = efun_incr(t_uninterrupt, a,b)*100

    # vlr_irr = 1 # f(i_ref)

    # dU_abs_pre = 0
    # t_uninterrupt = 0 # time increment of uninterrupted operation

    vlr_i =  vlr_lin_i(pec,i) # Factor accounting for influence of curent density on vlr

    # t_uni in sec; since tfrc_vlr_rev in 1/s
    vlr_rev = efun_incr(t_uninterrupt*3600, vlr_i*pec.frc_vlr_rev, pec.tfrc_vlr_rev)

    vlr_irr = vlr_i * (1-pec.frc_vlr_rev) # 0.4

    f_T = fctr_vlr_T(pec,T) # Factor accounting for influence of temperature on vlr

    f_i_ref = fctr_i_ref(pec,i) # Factor accounting for impact of vlr at different current densities
    dt = dt_in/3600 # Time increment // in h


    dU_rev = vlr_rev * dt * f_T #* f_i_ref               # Reversible voltage increase in current timestep
    dU_irr = vlr_irr * dt * f_T #* f_i_ref               # Irreversible voltage increase in current timestep
    dU_act = (dU_abs_pre + dU_rev + dU_irr) * f_i_ref    # Total voltage increase in curent timestep
    dU_abs = dU_abs_pre + dU_irr                        # Absolute voltage increase at current state of plant

    if print_vals:
        print('i: ', i)
        print('dt: ', dt)
        print('vlr_i: ', vlr_i)
        print('vlr_rev: ', vlr_rev)
        print('vlr_irr: ', vlr_irr)
        print('f_T: ', f_T)
        print('f_i_ref: ', f_i_ref)

    return dU_abs, dU_act, dU_rev, dU_irr

def fctr_i_ref(pec,i):
    '''


    Parameters
    ----------
    i : TYPE
        DESCRIPTION.
        in A/cm2

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''

    # slope =  -9.666182469382544
    # intercept =  35.76183499231406
    # slope = pec.slope_vlr_i_ref
    # intercept = pec.intercpt_vlr_i_ref
    fctr_vlr= (pec.slope_vlr_i_ref*i +pec.intercpt_vlr_i_ref)/pec.nrmdiv_vlr_i_ref
    if True:
        fctr_vlr = fctr_vlr if fctr_vlr >pec.lolim_vlr_i_ref else pec.lolim_vlr_i_ref
    if False:
        fctr_vlr = fctr_vlr if fctr_vlr <pec.hilim_vlr_i_ref else pec.hilim_vlr_i_ref
    return fctr_vlr # (intercept + slope*i)/intercept


def vlr_lin_i(pec,i):
    # lolim_vlr = 4 # Lower limit of vlr @ i = 0 // minimal value from Buttler
    # slope =  7.724910394265227
    # intercept =  66.78333333333335

    vlr= (pec.slope_vlr_i*i +pec.intercpt_vlr_i)/pec.nrmdiv_vlr_i
    if True:
        vlr = vlr if vlr >pec.lolim_vlr_i else pec.lolim_vlr_i
    if False:
        vlr = vlr if vlr <pec.hilim_vlr_i else pec.hilim_vlr_i
    return vlr if i > 0 else pec.lolim_vlr_i


def fctr_vlr_T(pec,T):#, vlr_min, vlr_max, T_min, T_max):

    #print('T: ', T)
    # vlr_min = 0
    # vlr_max = 50
    # T_min = 300
    # T_max = 353
    # f = (vlr_max-vlr_min)/(T_max-T_min) * (T - T_min) + vlr_min
    # f = f if T> T_min else vlr_min
    # f = f if T< T_max else vlr_max
    fvlr= (pec.slope_vlr_T*T +pec.intercpt_vlr_T)/pec.nrmdiv_vlr_T
    # print('fvlr: ', fvlr)
    if True:
        fvlr = fvlr if fvlr >pec.lolim_vlr_T else pec.lolim_vlr_T
    if False:
        fvlr = fvlr if fvlr <pec.hilim_vlr_T else pec.hilim_vlr_T
    return fvlr


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
