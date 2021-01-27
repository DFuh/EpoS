'''
AEL
calculation: polarisation characteristics
'''
import numpy as np

print(__name__ + ' imported...')



def testf(x, pec):
    print(f'x:{x} * par: {pec.b}:' , x*pec.b)
    return

def pwr_cell(u, i, A_cell=None, N_cells=None):
    '''
    calc (specific) power of single or multiple cell(s)
    '''
    if not A_cell:
        p = u * i
    else:
        p = u * i * A_cell

    if N_cells:
        return p * N_cells
    else:
        return p


def voltage_cell(obj, pec,
                appl_cv_rev=True, appl_ov_act=True, appl_ov_cnc=False,
                appl_ov_ohm=True):
    '''
    calculate cell voltage
    - Gibbs
    - Nernst
    - Activation overvoltage
    - dV due to Ohmic losses
    '''
    ### Reversible cell voltage
    ((dE_rev_ca, dE_rev_an), dE_rev, U_tn)      = cv_rev(obj, pec, T, pp)

    ### Activation overpotential
    # TODO: p_sat_KOH required
    U_act_ca, U_act_an = ov_act(obj,pec, T, pp)

    ### Concentration overpotential
    U_conc_ca, U_conc_an = ov_cnc(apply_funct=False)

    ### Additional Voltage due to Ohmic losses
    U_ohm_ca, U_ohm_an = None


    U_ca = U_rev_ca + U_act_ca + U_ohm_ca #dG_ca / (z*F) # cathodic halfcell potential
    U_an = U_rev_an + U_act_an + U_ohm_an # Anodic halfcell potential

    #U_rev, U_tn = None

    U_cell = U_ca + U_an

    return


def cv_rev(obj, pec, T, pp_H2_ca, pp_O2_an):
    '''
    calc reversible cell voltage
    '''
    '''
    @ standard pressure and temp.
    Schalenbach 2013 eq. 6 -11 (!)
    '''

    ((dG_ca, dG_an), (dH_ca, dH_an)) = clc_gibbs_free_energy(obj, pec, T)

    dE_ca = dG_ca / (2 * pec.F)
    dE_an = dG_an / (2 * pec.F)

    ### thermoneutral voltage
    U_tn = (dH_an + dH_ca) /(2 * pec.F)

    ### Nernst voltage
    '''
    deviations from temp., pressure (conc.)
    '''
    dE_N_ca = pec.R * T / (2 * pec.F) * np.log(pp_H2_ca /pec.p0_ref )
    dE_N_an = pec.R * T / (2 * pec.F) * np.log((pp_O2_an /pec.p0_ref)**(1/2))

    dE_rev_ca = dE0_ca + dE_N_ca
    #print(f'dE_N_ca: {dE_N_ca}')
    dE_rev_an = dE0_an + dE_N_an
    #print(f'dE_N_an: {dE_N_an}')
    dE_rev = dE_rev_ca - dE_rev_an
    #print(f'dE_rev: {dE_rev}')
    return ((dE_rev_ca, dE_rev_an),dE_rev,U_tn)


def ov_act(obj, pec, T, spaceholder,pp):
    '''
    calculate activation overvoltage
    '''

    pp_H2_ca, pp_O2_an, p_sat_KOH = pp[:2]


    i01   = i0_A + i0_C
    if i_in >0:

        ### exchange current density
        i0_an = pec.gamma_an * np.exp( (-pec.dG_c_an / pec.R) * ( (1 / T) - (1 / pec.T_ref_i0_an) )) * pec.i0_ref_an  # Abdin eq. 26
        i0_ca = pec.gamma_ca * np.exp( (-pec.dG_c_ca / pec.R) * ( (1 / T) - (1 / pec.T_ref_i0_ca) )) * pec.i0_ref_ca


        ### bubble coverage  theta
        # See: Olivier eq. 23 ff.
        #theta      = adv.f_thet(T,i,p_H2O) ### in aux

        # TODO: double check pressure vals used here
        theta_an = (-97.25 + 182 * (T / pec.T_ref_bc_an) - 84 * ( T / pv.T_ref_bc_an)**2) * ((i / pec.i_lim_bc)**0.3)*(p_an/(p_an - p_sat_KOH))
        theta_ca = (-97.25 + 182 * (T / pec.T_ref_bc_ca) - 84 * ( T / pv.T_ref_bc_ca)**2) * ((i / pec.i_lim_bc)**0.3)*(p_ca/(p_ca - p_sat_KOH))
        bcf_an      =  (1 - theta_an)    # caution: use different theta for an/ ca if pressure is different!
        bcf_ca      =  (1 - theta_ca)

        #alph_A      = 0.0675 + 0.00095 * T
        #alph_C      = 0.1175 + 0.00095 * T # HAmmoudi eq. 30

        ### ###
        dU_act_an   = ( (pec.R * T) / (pec.alpha_an * pec.F)) * np.log( i / (i0_an * bcf_an)) #Abdin eq. 29

        dU_act_ca   = ( (pec.R * T) / (pec.alpha_ca * pec.F)) * np.log( i / (i0_ca * bcf_ca)) #Abdin eq. 30


    else:
        dU_act_an   = 0
        dU_act_ca   = 0
    #dU_act     = dU_act_an + dU_act_ca
    return (dU_act_an, dU_act_ca)

def ov_cnc(apply_funct=True):

    '''
    ### Concentration Overpotential dU_con
    conc_O2_el_an = conc_O2_ch_an + ((d_an * n_O2_an) / (D_eff_an))
    conc_O2_ref =
    conc_H2_el_ca = conc_O2_ch_ca + ((d_ca * n_O2_an) / (D_eff_ca))
    conc_H2_ref =
    dU_con     = R*T/(4*F) * np.log(mol_con(T,i)[1]/mol_con(T,i)[0]) -  R*T/(2*F) * np.log(mol_con(T,i)[3]/mol_con(T,i)[2])
    dU_con     = -R*T/(2*F) * np.log((1 - i*R*T*dlt_CA/(2*F*D_eff(T)[1]*part_press(T,i)[2]))/(1+i*R*T*dlt_CA/(2*F*D_eff(T)[1]*part_press(T,i)[0])))    # NI_ SOL model

    '''
    if apply_funct:
        ret = None, None
    else:
        ret = 0,0
    return ret


def ov_ohm(T,i,):
    '''
    calc. additional cell-voltage-component due to ohmic losses
    '''
    ### Ohmic overpotential dU_ohm
    # TODO: check, if complete -> Hammoudi, Abdin
    # R_electrode |Henao 2013 eq. 13,14

    # calc conductivity of current collector
    sigma_ni = (6000000 - 279650*T + 532*T**2 - 0.38057*T**3) # Conductivity of nickel current collector // in S/cm
    sgm_crr = 1e2
    R_AN = 1 / (sgm_ni * sgm_crr) * (pv.dlt_AN / pv.srf_AN) # in 1/S
    R_CA = 1 / sgm_ni * (pv.dlt_CA / pv.srf_CA)

    ### calc. resistance of separator

    ### calc. resistance of electrolyte
    '''
    value of rho_el for eq. 60 in Abding unknown ???
    '''

    ### calc. resistance of electrodes


    # resistance electrolyte | Henao2013
    '''
    c_KOH     = 1.26054 * 1000 * pv.w_KOH/(100 * pv.M_KOH)  # molar concentration KOH at 30w% KOH solution   mol/m³
    sgm_KOH_free = (-2.04*c_KOH - 0.0028*c_KOH**2 + 0.005332*c_KOH*T + 207.2 * c_KOH/T + 0.001043*c_KOH**3 - 0.0000003*c_KOH**2*T**2) # // in S/m Olivier eq. 15

    R_el_free = 1 / sgm_KOH_free * (pv.l_AN_S / pv.srf_AN + pv.l_CA_S / pv.srf_CA) # // in 1/S
    R_el_b    = R_el_free * ((1/((1 - 2 * av.theta / 3)**(3/2))) -1)
    '''
    #R_mem = (0.060 + 80*np.exp(T/50))/(10000*srf_SE*10**4)            #//henao
    #R_mem = pv.rho_SE * pv.tau_SE**2 * pv.dlt_SE / (pv.ome_SE * pv.eps_SE * pv.srf_SE)             #//abdin
    #sigma_ely = None # ?

    ### Calculations adopted from Haug
    '''
    check: phi !!! (Haug below eq. 58)
    '''
    R_elec_an = (1/res_elec_an) * (1/(1-epsilon_an)**(3/2)) * (d_an / A_e)
    R_elec_ca = (1/res_elec_ca) * (1/(1-epsilon_ca)**(3/2)) * (d_ca / A_e)
    R_ely_an = (1/res_ely_ref_an) * ( (1 / (1 + kappa_ely*(T - T_ref)))
                * ( ((l_s_an - beta_an ) / A_e ) + (1/(1 - phi_an)**(3/2)) * beta_an/A_e  ) )
    R_ely_ca = (1/res_ely_ref_ca) * ( (1 / (1 + kappa_ely*(T - T_ref)))
                * ( ((l_s_ca - beta_ca ) / A_e ) + (1/(1 - phi_ca)**(3/2)) * beta_ca/A_e  ) )
    kappa_KOH = clc_KOH_conductivity(T, w_KOH) # based on Gilliam 2007
    R_sep = (1/kappa_KOH) * (tau_sep**2 * d_sep) / (omega_sep * epsilon_sep * A_sep)

    dU_ohm = i * (R_elec_an + R_elec_ca + R_ely_an + R_ely_ca + R_sep)
    #dU_ohm = i * pv.srf_CA * (R_AN + R_CA + R_el_free + R_el_b + R_mem + av.R_sep) # // in V

    return(dU_ohm)


### auxilliary calculations
################################################################################

def clc_KOH_conductivity(T, w_KOH):
    '''
    based on Gilliam 2007
    '''
    rho_KOH = clc_rho_KOH(T, w_KOH) # // in kg/m³ = g/L
    #m_KOH = (wt%*rho_KOH) / (100 * M_w)
    # w_KOH # in 1 (not % !)
    M_w = 56.10564 # Molar mass of potassium hydroxide // in g/mol
    m_KOH = (w_KOH * rho_KOH) / (M_w) # molarity KOH // in mol/L
    coeff = (-2.041, -0.0028, 0.005332, 207.2, 0.001043, -3*1e-7) # Coefficeints
    kappa_KOH = (coeff[0]*m_KOH + coeff[1] * m_KOH**2
                + coeff[2]* (m_KOH * T) + coeff[3] * (m_KOH / T)
                + coeff[4] * m_KOH**3 + coeff[5] * (m_KOH**2 * T**2))
    return kappa_KOH # S/cm

def clc_rho_KOH(T, w_KOH):
    ''' adopted from Haug (?) ->CHECK! '''
    ''' calc density of aqueous KOH-solution'''
    '''valid: 0.01 ... 200 °C /// w = 0 ... 0.5 *100 wt% KOH'''
    #w_KOH = 0.3 # mass fraction potassium hydroxide // in 1

    T_K0 = 273.15
    theta = T-T_K0
    rL = (1001.53053, -0.08343, -0.00401, 5.51232 *1e-6, -8.20994*1e-10)

    rho_int  = 0
    for j in range(5):
        rho_int += rL[j]*theta**j
    rho_L_out = rho_int * np.exp(0.86*w_KOH)
    return rho_L_out # in kg/m³


def clc_gibbs_free_energy(obj, pec, T):
    '''
    orig version in mod_3
    '''
    dH_ca = 0
    dG_ca = 0
    dH_an = (((1 * pec.H_H2) + ((1/2) * pec.H_O2))-(1 * pec.H_H2Ol))
    dG_an = ( dH_an -( ((1  *pec.S_H2) + ((1/2) * pec.S_O2) - (1 * pec.S_H2Ol))*T) )  #Carmo eq. 9 // Marangio eq. 6
    #dGa = -(-237.19*1e3)
    #U_rev = dGa / (2 * F)
    #U_tn = H_H2O /(2 * F) #- ???
    return ((dG_ca, dG_an), (dH_ca, dH_an))# // in
