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


def voltage_cell(obj, pec, T, i, p, pp=None, ini=False,
                appl_cv_rev=True, appl_ov_act=True, appl_ov_cnc=False,
                appl_ov_ohm=True):
    '''
    calculate cell voltage
    - Gibbs
    - Nernst
    - Activation overvoltage
    - dV due to Ohmic losses
    '''
    p_ca, p_an = p.cathode, p.anode

    if pp is None:
        pp = obj.clc_m.flws.partial_pressure(obj, pec, T, p)

    #if not A_cell:
    #A_cell = obj.pcll.active_cell_area

    ### Reversible cell voltage
    ((dE_rev_ca, dE_rev_an), dE_rev, U_tn)      = cv_rev(obj, pec, T, pp)

    ### Activation overpotential
    # TODO: p_sat_KOH required
    U_act_ca, U_act_an = ov_act(obj,pec, T, i, p, pp)

    ### Concentration overpotential
    U_conc_ca, U_conc_an = ov_cnc(apply_funct=False)

    ### Additional Voltage due to Ohmic losses
    U_ohm_ca, U_ohm_an, U_ohm_sep = ov_ohm(obj,pec, T, i, pp)


    U_ca = dE_rev_ca + U_act_ca + U_ohm_ca #dG_ca / (z*F) # cathodic halfcell potential
    U_an = dE_rev_an + U_act_an + U_ohm_an # Anodic halfcell potential

    #U_rev, U_tn = None

    U_cell = U_ca + U_an + U_ohm_sep
    print('U_rev={0} |U_act_An={1} |U_act_ca={2} |U_ohm_sep={3}'.format(dE_rev, U_act_an, U_act_ca, U_ohm_sep))
    return (U_ca, U_an, U_cell)


def cv_rev(obj, pec, T, pp):
    '''
    calc reversible cell voltage
    '''
    '''
    @ standard pressure and temp.
    Schalenbach 2013 eq. 6 -11 (!)
    '''
    #pp_H2_ca, pp_O2_an, pp_H2O = pp
    pp_H2_ca, pp_O2_an = pp[:2]

    ((dG_ca, dG_an), (dH_ca, dH_an)) = clc_gibbs_free_energy(obj, pec, T)

    dE_ca = dG_ca / (2 * pec.F)
    dE_an = dG_an / (2 * pec.F)

    ### thermoneutral voltage
    U_tn = (dH_an + dH_ca) /(2 * pec.F)

    ### Nernst voltage
    '''
    deviations from temp., pressure (conc.)
    '''
    if pp_H2_ca>0:
        dE_N_ca = pec.R * T / (2 * pec.F) * np.log(pp_H2_ca /pec.p0_ref )
    else:
        dE_N_ca = 0
    if pp_O2_an >0:
        dE_N_an = pec.R * T / (2 * pec.F) * np.log((pp_O2_an /pec.p0_ref)**(1/2))
    else:
        dE_N_an = 0
    dE_rev_ca = dE_ca + dE_N_ca
    #print(f'dE_N_ca: {dE_N_ca}')
    dE_rev_an = dE_an + dE_N_an
    #print(f'dE_N_an: {dE_N_an}')
    dE_rev = dE_rev_ca - dE_rev_an
    #print(f'dE_rev: {dE_rev}')
    return ((dE_rev_ca, dE_rev_an),dE_rev,U_tn)


def ov_act(obj, pec, T, i, p, pp):
    '''
    calculate activation overvoltage
    '''

    pp_H2_ca, pp_O2_an, p_sat_KOH = pp#[:2]

    clc_bubble_cvrg(obj, pec, T, i, p, pp)
    #print('i: ', i)
    #i01   = i0_A + i0_C
    if i >0:

        ### exchange current density
        i0_an = (pec.gamma_an *
                np.exp( (-pec.Ae_an / pec.R) * ( (1 / T) - (1 / pec.T_ref_i0_an) ))
                * pec.i0_ref_an ) # Abdin eq. 26
        i0_ca = (pec.gamma_ca *
                np.exp( (-pec.Ae_ca / pec.R) * ( (1 / T) - (1 / pec.T_ref_i0_ca) ))
                * pec.i0_ref_ca )


        ### bubble coverage  theta
        # See: Olivier eq. 23 ff.
        #theta      = adv.f_thet(T,i,p_H2O) ### in aux
        #print(' ------ |    i = ', i)
        #print('i0_an = ', i0_an)
        #print('i0_ca = ', i0_ca)
        # ??? if not hasattr(obj.av, 'theta_an'):


        bcf_an      =  (1 - obj.av.theta_an)    # caution: use different theta for an/ ca if pressure is different!
        bcf_ca      =  (1 - obj.av.theta_ca)
        #print('bcf_an = ', bcf_an)
        #print('bcf_ca = ', bcf_ca)

        #alph_A      = 0.0675 + 0.00095 * T
        #alph_C      = 0.1175 + 0.00095 * T # HAmmoudi eq. 30

        #print('i0 an/ca: ', i0_an, i0_ca)
        #print('bcf an/ca: ', bcf_an, bcf_ca)
        ### ###
        dU_act_an   = ( (pec.R * T) / (pec.alpha_an * pec.F)) * np.log( i / (i0_an * bcf_an)) #Abdin eq. 29

        dU_act_ca   = ( (pec.R * T) / (pec.alpha_ca * pec.F)) * np.log( i / (i0_ca * bcf_ca)) #Abdin eq. 30
        #print('dU_act_an = ', dU_act_an)
        #print('dU_act_ca = ', dU_act_ca)

    else:
        dU_act_an   = 0
        dU_act_ca   = 0
    #dU_act     = dU_act_an + dU_act_ca
    return (dU_act_ca, dU_act_an)

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
        dU_cnc_ca, dU_cnc_an = None, None
    else:
        dU_cnc_ca, dU_cnc_an = 0,0
    return (dU_cnc_ca, dU_cnc_an)

def ov_ohm(obj, pec, T, i, pp):
    '''
    Calc. ohmic conbtributions to cell voltage
    Adopted from Abdin 2017 and Henao 2013
    -> normalized to cell area (removed A_cell, A_electrode in denominator)

    -> See also Haug 2017 on gas voidage, phi, etc... !
    '''
    ### electrodes
    # calc conductivity of current collector |Henao 2013
    sgm_ni = (6000000 - 279650*T + 532*T**2 - 0.38057*T**3) *1e2 # Conductivity of nickel current collector // in S/m

    # Ohmic resistance of electrodes
    f_geom_an = 1/ (1 - pec.epsilon_an)**(3/2) # Geometry factor anode
    f_geom_ca = 1/ (1 - pec.epsilon_ca)**(3/2) # Geometry factor anode
    R_lctr_an = (1 / (sgm_ni)) * f_geom_an * (pec.d_lctr_an / pec.srf_lctr_an) # in 1/S
    R_lctr_ca = (1 / sgm_ni) * f_geom_ca * (pec.d_lctr_ca / pec.srf_lctr_ca)

    ### Ohmic resistance of electrolyte
    # bubble free electrolyte
    kappa_KOH = clc_conductivity_KOH(obj, pec, T, pec.w_KOH) *1e2 # // in S/m | Gilliam 2007
    # TODO: what, if zero-gap -> width of bubble zone??
    R_ely_free_an = (1/kappa_KOH) * (pec.l_anlctr_sep * (1 - pec.f_lbz_an ))#/ (pec.A_lctr_an) # without 1/A_electrode -> Ohm*m²
    R_ely_free_ca = (1/kappa_KOH) * (pec.l_calctr_sep * (1 - pec.f_lbz_ca ))#/ (pec.A_lctr_ca)

    # TODO: check order of calc. -> theta !
    f_geom_bc_an = 1/ (1 - obj.av.theta_an)**(3/2) # Geometry factor anode
    f_geom_bc_ca = 1/ (1 - obj.av.theta_ca)**(3/2) # Geometry factor anode
    R_ely_bc_an = (1/kappa_KOH) * f_geom_bc_an * pec.l_anlctr_sep * pec.f_lbz_an #/ (pec.A_lctr_an) # without 1/A_electrode -> Ohm*m²
    R_ely_bc_ca = (1/kappa_KOH) * f_geom_bc_ca * pec.l_calctr_sep * pec.f_lbz_ca #/ (pec.A_lctr_ca)

    #print('R_ely_free_an',R_ely_free_an)
    #print('R_ely_bc_an', R_ely_bc_an)
    R_ely_an = R_ely_free_an + R_ely_bc_an
    R_ely_ca = R_ely_free_ca + R_ely_bc_ca
    ### Ohmic resistance of separator
    R_sep = ((1/kappa_KOH) * (pec.tau_sep**2 * pec.d_sep) /
            (pec.omega_sep * pec.epsilon_sep)) #* pec.A_sep) -> Ohm*m²

    dU_ohm_an   = i * (R_lctr_an + R_ely_an )
    dU_ohm_ca   = i * (R_lctr_ca + R_ely_ca )
    dU_ohm_sep  = i * R_sep
    print('dU_ohm_ca={0}| dU_ohm_an={1}| dU_ohm_sep={2}'.format(dU_ohm_ca, dU_ohm_an, dU_ohm_sep))
    return (dU_ohm_ca, dU_ohm_an, dU_ohm_sep)


### auxilliary calculations
################################################################################

def clc_conductivity_KOH(obj, pec, T, w_KOH):
    '''
    Calc conductivity of aqous KOH solution
    based on Gilliam 2007

    Parameters
    ----------
    obj: Simulation instances
    pec: named tuple, containing electrochemical parameters
    T: Temperature in K
    w_KOH: weight percent KOH in 1 (31wt% -> w_KOH=0.31)

    Returns
    -------
    kappa_KOH: float
        unit: S/cm
    '''
    if not hasattr(obj, 'rho_KOH'):
        obj.av.rho_KOH = clc_rho_KOH(T, w_KOH) # // in kg/m³ = g/L
    #m_KOH = (wt%*rho_KOH) / (100 * M_w)
    # w_KOH # in 1 (not % !)
    #M_w = 56.10564 # Molar mass of potassium hydroxide // in g/mol
    m_KOH = (w_KOH * obj.av.rho_KOH) / (pec.M_KOH*1e3) # molarity KOH // in mol/m³
    coeff = (-2.041, -0.0028, 0.005332, 207.2, 0.001043, -3*1e-7) # Coefficeints
    kappa_KOH = (coeff[0]*m_KOH + coeff[1] * m_KOH**2
                + coeff[2]* (m_KOH * T) + coeff[3] * (m_KOH / T)
                + coeff[4] * m_KOH**3 + coeff[5] * (m_KOH**2 * T**2))
    return kappa_KOH # S/cm


def clc_rho_KOH(T, w_KOH):
    '''
    =======> DUPLICATE <====== (clc_auxvals in ael_aux_vXX)

    Calculate density of aqousPotassium hydroxide solution
    adopted from Haug
    valid: 0.01 ... 200 °C /// w = 0 ... 0.5 *100 wt% KOH

    Parameters
    ----------
    T: Temperature in K
    w_KOH: wt% KOH in 1 #w_KOH = 0.3 # mass fraction potassium hydroxide // in 1

    Returns
    --------
    rho_L_out: float
        unit: kg/m³

    '''


    T_K0 = 273.15
    theta = T-T_K0
    rL = (1001.53053, -0.08343, -0.00401, 5.51232 *1e-6, -8.20994*1e-10)

    rho_int  = 0
    for j in range(5):
        rho_int += rL[j]*theta**j
    rho_L_out = rho_int * np.exp(0.86*w_KOH)
    return rho_L_out # in kg/m³


def clc_bubble_cvrg(obj, pec, T, i, p, pp):
    '''
    Calc. Bubble coverage of electrodes based on Abdin 2017 eq. 31

    '''
    #TODO: check order of calculation (prevent redundant calc. !)
    #p_sat_KOH = xflws.clc_pp_H2O()
    p_sat_KOH = pp[-1]

    # TODO: double check pressure vals used here

    obj.av.theta_an = ((-97.25 + 182 * (T / pec.T_ref_bc_an) - 84 * ( T / pec.T_ref_bc_an)**2)
                * ((i / pec.i_lim_bc)**0.3)*(p.anode/(p.anode - p_sat_KOH)))

    obj.av.theta_ca = ((-97.25 + 182 * (T / pec.T_ref_bc_ca) - 84 * ( T / pec.T_ref_bc_ca)**2)
                * ((i / pec.i_lim_bc)**0.3)*(p.cathode/(p.cathode - p_sat_KOH)))
    return

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
