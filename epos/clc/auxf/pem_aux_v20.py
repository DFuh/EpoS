'''
general auxilliary calculations
'''


import numpy as np

def clc_auxvals(obj, T):
    '''
    Calc auxilliary values
    '''
    clc_lambda_mem(obj, T)
    clc_c_hpl(obj, T)
    clc_rho_H2O(obj, T)
    clc_pp_H2O(obj, obj.pec, T)
    return


# ==============================================================================

def clc_lambda_mem(obj, T):
    '''
    Calc lambda of Membrane for respective Temnperature (in K)
    based on Ito 2011 (lin fit on Fig 1); valid for T<=100°C
    lin fit
    -> Yoshitake        | m,b = [  0.11272062, -19.34588299] /// for °C: [ 0.11272062, 11.44375396]
    -> Hinatsu (N-Form) | m,b = [  0.13319059, -26.77932753] /// for °C: [ 0.13319059, 9.60168093]
    '''

    ### Yoshitake
    m,b = [  0.11272062, -19.34588299]
    obj.av.lambda_mem = T*m + b
    return


def clc_c_hpl(obj,T):
    ### concentration of H+-ions in Naf
    ### see Sethuraman2008 eq. 17...19
    #T = 353
    #lambd   = lam(T) # sethuraman eq 15/16
    #lambd   = 14
    EW      = 1100 # // in g/mol resp.: g/ equiv.
    #c_hp    = (1980 + 32.4 * lambd) / ( (1 + (0.0648 * lambd)) * EW) # // in mol/kg ???
    c_hp    = ((1.98 + 0.0324 * obj.av.lambda_mem) / (1 + (0.0648 * obj.av.lambda_mem))) * 1/EW # // in mol/cm³ ???
    obj.av.c_hpl = c_hp*1e3*1e3 #+++ edit 202106: factor: 1e3 -> change mol/l to mol/m³
    #return c_hp*1000#*100e3/1000 # in mol/kg +++changed to mol/m³ (see above)
    return

def clc_rho_H2O(obj, T):
    obj.av.rho_H2O= obj.av.rho_ely = 999.972 - 7*10**(-3)*(T-273.15-20)**2 # Source=? factors "20" vs "4" ??? (both not valid)
    return

def clc_pp_H2O(obj, pec, T):
    '''
    Partial pressure of water vapor,
     # espinoza-lopez // in Pa
    '''
    obj.av.pp_H2O = 1e5 *(6.1078 * 1e-3 * np.exp( 17.2694 * (T - 273.15) / (T - 34.85) ))
    return 

### miscellaneous

def division(n, d):
    return n / d if d else 0
