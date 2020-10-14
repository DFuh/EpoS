'''
calculation: polarisation characteristics
'''
import numpy as np


print(__name__ + ' imported...')



def partial_pressure(T, p_in):
    ''' partial pressure of product gases dependend on water-vapor-pressure'''
    #T in K
    p_ca, p_an = p_in


    #p_H2O
    pp_H2O    = 1e5 *(6.1078 * 1e-3 * np.exp( 17.2694 * (T - 273.15) / (T - 34.85) )) # espinoza-lopez // in Pa
    #p_O2_an
    #if i > 0: # +++++ edit: 2019-10-28
    #av.plr_ppr[0][1] = p_A - p_H2O
    pp_H2_ca = p_ca - pp_H2O
        #p_H2_ca
    #av.plr_ppr[0][2] = p_C - p_H2O
    pp_O2_an = p_an - pp_H2O
    #else:
    #av.plr_ppr[0][1] = 0#p_A - p_H2O
        #p_H2_ca
    #av.plr_ppr[0][2] = 0#p_C - p_H2O

    #av.plr_ppr[1] = p_H2O
    return pp_H2_ca, pp_O2_an, pp_H2O #av.plr_ppr[0][1:3]
