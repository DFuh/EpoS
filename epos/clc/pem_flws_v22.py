'''
calculation: flow conditions in PEMEL-cell
literature:
 - Marangio_2009: Theoretical Model and experimental analysis of high pressure PEM water electrolyser for hydrogen production
 - Abdin 2015 (same approach, but different calc. of binary diffusion coefficients)

'''

import numpy as np


print(__name__ + ' imported...')

k_L_H2_ver = 2 # How to calc. k_L_H2 || 0 -> const. | 1 -> lin. | 2 -> sqrt
diff_coeff_H2_ver = 0 # How to calc. D_H2 || 0 -> Wise | 1 -> Bernardi | 2 -> Chandesris |

def clc_flws_auxvals(obj, T, i):
    '''
    calculate all necessary auxilliary values for flw-balance
    '''
    ### D_eff(T)
    Diff_corr = (eps / tau)
    obj.av.D_eff_H2 =

    ### k_L(i)

    obj.av.k_L_H2 = clc_k_L(i, k_L_H2_ver)

    ###

    return

def clc_K_L_H2(i, ver=None):
    '''
    Trinke 2017 (fig 8)
    ------
    Returns
    k_L in m/s
    '''
    i_cm2 = i/1e4 # Convert current density to A/cm²
    if ver == 2:
        i_cm2 = i/1e4 # Convert current density to A/cm²
        k_L = np.sqrt(1.2 * 1e-5 * i_cm2)
    elif ver == 1:
        m = 1.3*1e-3 # m/s cm²/A
        b = 1.9*1e-3 # m/s
        k_L = m*i_cm2 +b
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

def mlr_flws_prod(i):
    '''
    area-specific molar flows from chemical reactions in cell
    Marangio_2009 eq.:
    '''
    N_H2_prd_ca = i / (2*F)     # Hydrogen production // in mol/(s m²)
    N_O2_prd_an = i / (4*F)     # Oxygen production // in mol/(s m²)
    N_H2O_cns = i / (2*F)       # Water consumption // in mol/(s m²)

    return

def mlr_flw_balance():
    '''
    flow balance in cell

    // in mol/s
    '''
    n_H2O_an        = n_H2O_m + n_H2O_cns
    n_H2O_an_out    = n_H2O_in - n_H2O_cns
    n_H2O_ca_out    = n_H2O_m

    return

def gas_permeation():
    '''
    gas permeation through membrane
    - diffusion (fick)
    - differential pressure (darcy)
    - electro osmotic drag (mainly affects Oxygen)
    (ion moving direction in PEM: an -> ca)

    '''

    #n_H2_prm = n_H2_prm_drc + n_H2_prm_fck
    n_H2_prm = obj.av.Deff_H2 * ( (i / (2*F)) + (obj.av.k_L * c_H2_henry)) / (obj.av.k_L * av.d_mem + Deff_H2) # // in mol/m²s | Trinke2017_H2-perm eq.11
    n_O2_prm = n_O2_prm_drg #
    return


def mlr_frc_ch():
    '''
    calc. molar fractions of species in outlet-channels of cell
    Marangio_2009 eq.: 37
    '''
    ### TODO: consider permeation?

    x_H2_ch_ca = n_H2_chca / (n_H2_chca + n_O2_chca + n_H2O_chca)
    x_H2_ch_an = n_H2_chan / (n_H2_chan + n_O2_chan + n_H2O_chan)

    x_O2_ch_ca = n_O2_chca / (n_H2_chca + n_O2_chca + n_H2O_chca)
    x_O2_ch_an = n_O2_chan / (n_H2_chan + n_O2_chan + n_H2O_chan)

    x_H2O_ch_ca = n_H2O_chca / (n_H2_chca + n_O2_chca + n_H2O_chca)
    x_H2O_ch_an = n_H2O_chan / (n_H2_chan + n_O2_chan + n_H2O_chan)

    return


def conc_ch(obj, pec, T, i):
    '''
    Concentration of species in channels
    Marangio 2009 eq. 38
    '''
    c_H2_ch_ca = p_ca * x_H2_ch_ca / (pec.R * T)
    c_O2_ch_an = p_an * x_O2_ch_an / (pec.R * T)

    c_H2O_ch = obj.av.rho_H2O / pec.M_H2O
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

    x_H2_mem_ca = + n_O2_chca + n_H2O_chca)
    x_H2_mem_an = + n_O2_chan + n_H2O_chan)

    x_O2_mem_ca = + n_O2_chca + n_H2O_chca)
    x_O2_mem_an = + n_O2_chan + n_H2O_chan)

    x_H2O_mem_ca =  + n_O2_chca + n_H2O_chca)
    x_H2O_mem_an =  + n_O2_chan + n_H2O_chan)

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


def crssflw_membrane_H2O(obj, pec):
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
    n_H2O_eo = pec.n_d * (i * pec.A_cell) / pec.F

    ### diffusion
    n_H2O_dd = (pec.A_cell * pec.D_wm) / obj.av.d_mem * (c_H2O_m_ca - c_H2O_m_an)

    ### differential pressure
    n_H2O_pe = pec.k_drc * (pec.A_cell * obj.av.rho_H2O) / (obj.av.viscos_H2O * M_H2O)
    # fix values: (?)
    # rho_H2O = 1000 kg/m³
    # viscos_H2O = 1.1*1e-3 Pa s
    return
