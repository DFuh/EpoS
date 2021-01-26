'''
calculation: flow conditions in PEMEL-cell
literature:
 - Marangio_2009: Theoretical Model and experimental analysis of high pressure PEM water electrolyser for hydrogen production
 - Abdin 2015 (same approach, but different calc. of binary diffusion coefficients)

'''

import numpy as np
#from collections import namedtuple

print(__name__ + ' imported...')

k_L_H2_ver = 2 # How to calc. k_L_H2 || 0 -> const. | 1 -> lin. | 2 -> sqrt
diff_coeff_H2_ver = 0 # How to calc. D_H2 || 0 -> Wise | 1 -> Bernardi | 2 -> Chandesris |

def matbal(i, p_an, p_ca):
    '''
    basic calculations for overall material balance
    '''

    ### pressures



    ###
    return


def clc_flws_auxpars(obj, T):
    '''
    calculate all necessary auxilliary parameters for flw-balance
    '''
    obj.av.d_mem = obj.pec.d0_mem
    obj.av.D_eff_H2, obj.av.D_eff_O2 = clc_diffusion_coefficient(obj, obj.pec, T)
    ### H2O
    obj.av.rho_H2O      = 999.972 - 7*10**(-3)*(T-273.15-20)**2 # Source=? factors "20" vs "4" ??? (both not valid)
    obj.av.viscos_H2O   = 1 / (0.1 * T**2 - 34.335 * T + 2472) # Dynamic viscosity of water // in Pa s | wikipedia
    ###

    return

##### clc Diff-Coefficient
''' output of dico-functions H2 , O2'''

def clc_diffusion_coefficient(obj, pec, T):
    ### Diffusion coefficient D_eff(T)
    #Diff_corr = 0 #(eps / tau)
    #obj.av.D_eff_H2 = 1 #None
    ### Diffusion ???
    dif_fun         = dico_wise, dico_Bernardi, dico_Chandesris
    #eps             = 0.37 # Trinke 2017
    #tau             = 1.5 # Trinke 2017
    D_H2, D_O2 = dif_fun[pec.clc_diffusion](pec,T) # *av.corr_Di * (eps/tau)
    D_eff_H2 = D_H2 * (pec.epsilon_ca / pec.tau_ca)
    D_eff_O2 = D_O2 * (pec.epsilon_an / pec.tau_an)
    return D_eff_H2, D_eff_O2

def dico_wise(pec,T):
    ### experiments: 10...60°C
    B = np.array([4.9,4.2]) *1e-2        # // in cm²/s error= [+-0.3 / +-0.2]
    calJ = 4.1858
    dE = np.array([3960,4390])*calJ # // in calJ*cal/mol error= [+-40 / +-20]
    D_l = 1e-4 * B * np.exp(-dE/(pec.R*T))
    return D_l #H2, O2

def dico_Bernardi(pec,T):
    B = np.array([4.1,3.1]) *1e-3        # // in cm²/s
    C = np.array([2602,2768])       # // in K (?)
    D_l = 1e-4*B * np.exp(-C/(T))
    return D_l #H2, O2

def dico_Chandesris(pec,T):
    B = np.array([1.23,4.1]) *1e-6        # // in cm²/s
    C = np.array([2602,18380/pec.R])       # // in K (?)
    D_l = B * np.exp(-C/(T))
    return D_l #H2, O2

def clc_K_L_H2(i, ver=None):
    '''
    Trinke 2017 (fig 8)
    ------
    Returns
    k_L in m/s
    '''
    i_cm2 = i/1e4 # Convert current density to A/cm²
    ### square root
    if ver == 2:
        #i_cm2 = i/1e4 # Convert current density to A/cm²
        k_L = np.sqrt(1.2 * 1e-5 * i_cm2)
    ### linear
    elif ver == 1:
        m = 1.3*1e-3 # m/s cm²/A
        b = 1.9*1e-3 # m/s
        k_L = m*i_cm2 +b
    ### constant
    else:
        k_L = 3*1e-3 # m/s
    return k_L # in m/s


def diff_coeff_corr_H2(ver=None):

    if ver == 1: # Marangio 2009  eq. 30
        #eps = 0.3       # Porosity of electrodes
        #eps_p = 0.11    # percolation treshold
        #alpha = 0.785   # Empirical coefficient
        dc = pec.epsilon_ca * ( (pec.epsilon_ca - pec.epsilon_perc_ca) /
                                (1 - pec.epsilon_perc_ca) )**pec.alpha_diff_corr
    else:
        dc = pec.epsilon_ca / pec.tau_ca
    return dc


def absolute_pressure():
    '''
    calc. absolute presure at electrodes
    '''

    return


def clc_mlr_flws_prod(obj, pec, i, A_cell):
    '''
    area-specific molar flows from chemical reactions in cell
    Marangio_2009 eq.:
    '''
    N_H2_prd_ca = i * A_cell/ (2*pec.F)     # Hydrogen production // in mol/(s)
    N_O2_prd_an = i * A_cell/ (4*pec.F)     # Oxygen production // in mol/(s)
    N_H2O_cns = i * A_cell/ (2*pec.F)       # Water consumption // in mol/(s)

    return N_H2_prd_ca, N_O2_prd_an, N_H2O_cns


def clc_mlr_flw_balance(i, n_H2O_ch_in_an):

    return



def clc_hydrogen_permeation(obj, pec, i, c_H2_henry=0):
    '''
    permeation of hydrogen through membrane
    - diffusion (fick)
    - differential pressure (darcy)
    - electro osmotic drag (mainly affects Oxygen)
    (ion moving direction in PEM: an -> ca)

    '''
    ### k_L(i)
    obj.av.k_L_H2 = clc_K_L_H2(i, k_L_H2_ver)

    #n_H2_prm = n_H2_prm_drc + n_H2_prm_fck
    # // in mol/m²s | Trinke2017_H2-perm eq.11
    n_H2_prm = obj.av.D_eff_H2 * (( (i / (2*pec.F)) + (obj.av.k_L_H2 * c_H2_henry)) /
                                        (obj.av.k_L_H2 * obj.av.d_mem + obj.av.D_eff_H2) )

    return n_H2_prm

def clc_oxygen_permeation(obj, ):
    '''
    permeation of oxygen through membrane
    - diffusion (fick)
    - differential pressure (darcy)
    - electro osmotic drag (mainly affects Oxygen)
    (ion moving direction in PEM: an -> ca)

    '''
    n_O2_prm = 0#n_O2_prm_drg #
    return n_O2_prm

def clc_mlr_frc(numer, denom=[]):
    if denom:
        sum = 0
        for denomi in denom:
            sum += denomi
        frc = numer / (numer + sum)
    else:
        frc = None
    return frc

def clc_mlr_frc_ch():
    '''
    calc. molar fractions of species in outlet-channels of cell
    Marangio_2009 eq.: 37
    '''
    ### TODO: consider permeation?

    x_H2_ch_ca = n_H2_ch_ca / (n_H2_ch_ca + n_O2_ch_ca + n_H2O_ch_ca)
    x_H2_ch_an = n_H2_ch_an / (n_H2_ch_an + n_O2_ch_an + n_H2O_ch_an)

    x_O2_ch_ca = n_O2_ch_ca / (n_H2_ch_ca + n_O2_ch_ca + n_H2O_ch_ca)
    x_O2_ch_an = n_O2_ch_an / (n_H2_ch_an + n_O2_ch_an + n_H2O_ch_an)

    x_H2O_ch_ca = n_H2O_ch_ca / (n_H2_ch_ca + n_O2_ch_ca + n_H2O_ch_ca)
    x_H2O_ch_an = n_H2O_ch_an / (n_H2_ch_an + n_O2_ch_an + n_H2O_ch_an)

    return

def clc_conc_px(pressure, mlrfrc):
    '''
    clc concentration based on pressure (in compartment)
    and molar fraction

    '''

    return (pressure * mlrfrc)/1 #???


def conc_ch(obj, pec, T, i):
    '''
    Concentration of species in channels
    Marangio 2009 eq. 38
    '''
    #c_H2_ch_ca = p_ca * x_H2_ch_ca / (pec.R * T)
    #c_O2_ch_an = p_an * x_O2_ch_an / (pec.R * T)

    #c_H2O_ch = obj.av.rho_H2O / pec.M_H2O
    return


def conc_mem(obj, pec, T, i):
    '''
    Concentration of species at membrane
    Marangio 2009 eq. 39
    '''


    c_H2_mem_ca = c_H2_ch_ca * (d_e_ca/D_eff_H2) * n_H2
    c_O2_mem_an = c_O2_ch_an * (d_e_an/D_eff_O2) * n_O2

    # TODO: check D_eff for H2O
    # liquid water?
    c_H2O_mem_ca = c_H2_ch_ca * (d_e_ca/D_eff_H2) * n_H2O_ca # Check !
    c_H2O_mem_an = c_O2_ch_an * (d_e_an/D_eff_O2) * n_H2O_an

    return


def mlr_frc_mem():
    '''
    calc. molar fractions of species at membrane
    Marangio_2009 eq.: 40
    '''
    ### TODO: consider permeation?

    #x_H2_mem_ca = + n_O2_chca + n_H2O_chca)
    #x_H2_mem_an = + n_O2_chan + n_H2O_chan)

    #x_O2_mem_ca = + n_O2_chca + n_H2O_chca)
    #x_O2_mem_an = + n_O2_chan + n_H2O_chan)

    #x_H2O_mem_ca =  + n_O2_chca + n_H2O_chca)
    #x_H2O_mem_an =  + n_O2_chan + n_H2O_chan)

    return


def partial_pressure(pec, T, p_in):
    ''' partial pressure of product gases dependend on water-vapor-pressure'''

    p_ca, p_an = p_in

    pp_H2_mem_ca = p_ca * x_H2_mem_ca
    pp_H2_mem_an = p_an * x_H2_mem_an

    pp_O2_mem_ca = p_ca * x_O2_mem_ca
    pp_O2_mem_an = p_an * x_O2_mem_an

    pp_H2O_


    return pp_H2_ca, pp_O2_an, pp_H2O #av.plr_ppr[0][1:3]


def clc_crssflw_membrane_H2O(obj, pec, i, T, A_cell, c_H2O_m_ca=0, c_H2O_m_an=0):
    '''
    water crossflow through membrane, due to:
    - electro osmotiv drag
    - diffusion (fick)
    - differential pressure (darcy)
    '''

    '''
    see Panchenko 2018 for gas/liq behaviour: O2-bubbles causing local dry-out of mem
    '''

    ### electro osmotic drag
    n_H2O_eo = pec.n_drag * (i * A_cell) / pec.F # an > ca

    ### diffusion
    n_H2O_dd = 0#(A_cell * pec.D_wm) / obj.av.d_mem * (c_H2O_m_ca - c_H2O_m_an)

    ### differential pressure ## ca -> an
    ''' delta-p missing for darcy ?! check Abdin!
    '''
    n_H2O_pe = - pec.k_drc * (A_cell * obj.av.rho_H2O) / (obj.av.viscos_H2O * pec.M_H2O)
    # fix values: (?)
    # rho_H2O = 1000 kg/m³
    # viscos_H2O = 1.1*1e-3 Pa s

    return n_H2O_eo + n_H2O_dd + n_H2O_pe

def clc_crssflw_membrane_H2O_marangio(obj, pec, i, T, A_cell,):
    '''
    Marangio eq. 36

    -> significant influence of temperature on diffusion, if implemented !

    yet to be checked:
    -> D_eff_i !
    '''
    n_H2O_eo = pec.n_drag * (i * A_cell) / pec.F # an > ca
    n_H2O_pe = - pec.k_drc * (A_cell * obj.av.rho_H2O) / (obj.av.viscos_H2O * pec.M_H2O)


    rho_H2O_ca = obj.av.rho_H2O # no temperature deviations from an to ca considered yet
    rho_H2O_an = obj.av.rho_H2O
    #rho_H2O_ca = 999.972 - 7*10**(-3)*(T-273.15-20)**2
    #rho_H2O_an = 999.972 - 7*10**(-3)*(T-273.15-20 -1)**2

    n_H2O_cns = (i*A_cell) / (2*pec.F)

    n_H2O_t_a = (n_H2O_eo - n_H2O_pe + (A_cell * obj.pec.D_wm)/ obj.av.d_mem
                 *((rho_H2O_ca - rho_H2O_an) / obj.pec.M_H2O)
                    + ((obj.pec.d_el_an * n_H2O_cns) / (obj.av.D_eff_O2 * A_cell) )
                 )
    n_H2O_t_b = 1 - ((obj.pec.D_wm / obj.av.d_mem)
                        * ((obj.pec.d_el_ca / obj.av.D_eff_H2)
                            + (obj.pec.d_el_an / obj.av.D_eff_O2)))
    n_H2O_t = n_H2O_t_a / n_H2O_t_b
    return n_H2O_t


def clc_crssflw_membrane_H2O_chandesris(obj, pec, i, T, A_cell,):
    '''
    water crossflow through membrane, due to:
    - electro osmotiv drag
    - diffusion (fick)
    - differential pressure (darcy)

    derived experimentally by Chandesris (2015)

    '''
    #def water_bal(i ,pv0, pv, av):
    ## i-input in A/m²
    '''water balance/drag |Chandesris '''
    ### TODO: validate!!!
    if i >0:
        N_bsc = i/(2*pec.F)                                      # // in mol/(m² s)
        #q_c_H2O = (M_H2O * i) / (z * F * rho_H2O)           # water consumption // in m³/(m² s)
        q_c_H2O = N_bsc                                      # water consumption // in mol/(m² s)
        #Q_t_H2O = (- 0.332 * np.log10(i) +5.59) * Q_c_H2O   # water perm through membrane // in m³/(m² s)
        #n_H2O_t = Q_t_H2O*(rho_H2O/M_H2O)                   # // in mol/(m² s)
        #n_H2O_t = (- 0.332 * np.log10(i) +5.59) *N_bsc      # // in mol/(m² s)
        N_H2O_t = (- 0.332 * np.log10(i) +5.59) *q_c_H2O     # // in mol/(m² s)
        #rho_H2O = obj.av.rho_H2O
        #v_H2O_t = (- 0.332 * np.log10(i) +5.59) *q_c_H2O* (pec.M_H2O / rho_H2O)  #* av.corr_v_H2O   # // in m³/(m² s) bzw. m/ s
        #av.v_H2O = v_H2O_t
        #print('v_H2O_t/q_c_H2O ( in water-bal)',v_H2O_t/q_c_H2O)
    else:
        q_c_H2O, v_H2O_t, N_H2O_t = 0,0,0
    n_H2O = q_c_H2O * A_cell # water consumption // in mol/s

    return N_H2O_t * A_cell # // in mol/s
