'''
general auxilliary calculations
'''
import numpy as np

def clc_auxvals(obj, T):
    '''
    Calc auxilliary values
    '''

    # obj.av.rho_H2O      = clc_rho_H2O(T)
    clc_rho_H2O(obj, T)
    # obj.av.rho_ely = clc_rho_KOH(obj, T, obj.pec.w_KOH)
    clc_rho_KOH(obj, T, obj.pec.w_KOH)
    return


# ==============================================================================
def clc_rho_KOH(obj,T, w_KOH):
    ''' calc density of aqueous KOH-solution
        valid: 0.01 ... 200 °C /// w = 0 ... 0.5 *100 wt% KOH'''
    #w_KOH = 0.3 # mass fraction potassium hydroxide // in 1

    T_K0 = 273.15
    theta = T-T_K0
    rL_0 = 1001.53053
    rL_1 = -0.08343
    rL_2 = -0.00401
    rL_3 = 5.51232 *1e-6
    rL_4 = -8.20994*1e-10

    rL = rL_0,rL_1,rL_2,rL_3,rL_4
    rho_int  = 0
    for j in range(5):
        rho_int += rL[j]*theta**j
    obj.av.rho_ely = rho_int * np.exp(0.86*w_KOH)
    return # rho_L_out # in kg/m³

def clc_rho_H2O(obj, T):
    obj.av.rho_H2O= 999.972 - 7*10**(-3)*(T-273.15-20)**2 # Source=? factors "20" vs "4" ??? (both not valid)
    return
### miscellaneous

def division(n, d):
    return n / d if d else 0
