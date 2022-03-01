'''
AEL
calculation: degaradation/ ageing effects
'''
# import numpy as np

print(__name__ + ' imported...')

def voltage_increase_lin(obj, pec, T, i):

    ### absolute voltage increase
    dU_dgr_act = obj.av.t_op/3600 * obj.pec.fctr_vlr # h * V/h

    ### incremental voltage increase
    # dU_dgr_incr = obj.av.t_diff/3600 * obj.pec.fctr_vlr # h * V/h

    return dU_dgr_act, 0 #dU_dgr_abs
