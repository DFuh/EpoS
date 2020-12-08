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


def voltage_cell(obj, pec,):
    '''
    calculate cell voltage
    - Gibbs
    - Nernst
    - Activation overvoltage
    - dV due to Ohmic losses
    '''
    ### Reversible cell voltage
    U_rev_ca = None
    U_rev_an = None

    ### Activation overpotential
    U_act_ca, U_act_an = None

    ### Concentration overpotential
    U_conc_ca, U_conc_an = 0,0

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

    dE_ca = pec.dG_ca / (2 * pec.F)
    dE_an = pec.dG_an / (2 * pec.F)

    ### thermoneutral voltage
    U_tn = dH_an/(2*pec.F)

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


def ov_act():
    '''
    calculate activation overvoltage
    '''



    i01   = i0_A + i0_C
    if i_in >0:

        ### exchange current density
        i0_an = pv.gma_M_an * np.exp( (-pv.dG_c_an / gpv.R) * ( (1 / T) - (1 / pv.T_ref_an) )) * pv.i0_ref_an  # Abdin eq. 26
        i0_ca = pv.gma_M_ca * np.exp( (-pv.dG_c_ca / gpv.R) * ( (1 / T) - (1 / pv.T_ref_ca) )) * pv.i0_ref_ca


        ### bubble coverage  theta
        # See: Olivier eq. 23 ff.
        #theta      = adv.f_thet(T,i,p_H2O) ### in aux
        theta = (-97.25 + 182 * (T / pv.T_0_b) - 84 * ( T / pv.T_0_b) **2) * ((i / pv.i_lim)**0.3)*(p/(p-p_H2O))
        bf          =  1 #(1-theta)    # caution: use different theta for an/ ca if pressure is different!


        #alph_A      = 0.0675 + 0.00095 * T
        #alph_C      = 0.1175 + 0.00095 * T # HAmmoudi eq. 30

        ### ###
        dU_act_A   = ( (gpv.R * T) / (pv.alph_A * gpv.F)) * np.log( i / (i0_An * bf))

        dU_act_A_b = 0#( (R*T) / (pv.alph_A * F) ) * np.log( 1 / (1-theta) )

        dU_act_C   = ( (gpv.R * T) / (pv.alph_C * gpv.F)) * np.log( i / (i0_Ca * bf))

        dU_act_C_b = 0#( (R*T) / (pv.alph_C * F) ) * np.log( 1 / (1-theta))

        dU_act     = dU_act_A + dU_act_C
        dU_act_b   = dU_act_A_b + dU_act_C_b
    else:
        dU_act     = 0
        dU_act_b   = 0
    return


def ov_ohm(T,i,):


    '''
    ### Concentration Overpotential dU_con

    dU_con     = R*T/(4*F) * np.log(mol_con(T,i)[1]/mol_con(T,i)[0]) -  R*T/(2*F) * np.log(mol_con(T,i)[3]/mol_con(T,i)[2])
    dU_con     = -R*T/(2*F) * np.log((1 - i*R*T*dlt_CA/(2*F*D_eff(T)[1]*part_press(T,i)[2]))/(1+i*R*T*dlt_CA/(2*F*D_eff(T)[1]*part_press(T,i)[0])))    # NI_ SOL model

    '''

    ### Ohmic overpotential dU_ohm
    # TODO: check, if complete -> Hammoudi, Abdin
    # R_electrode |Henao 2013 eq. 13,14

    sgm_ni = (6000000 - 279650*T + 532*T**2 - 0.38057*T**3) # // in S/cm
    sgm_crr = 1e2
    R_AN = 1 / (sgm_ni * sgm_crr) * (pv.dlt_AN / pv.srf_AN) # in 1/S
    R_CA = 1 / sgm_ni * (pv.dlt_CA / pv.srf_CA)

    # R_electrolyte | Henao2013

    c_KOH     = 1.26054 * 1000 * pv.w_KOH/(100 * pv.M_KOH)  # molar concentration KOH at 30w% KOH solution   mol/m³
    sgm_KOH_free = (-2.04*c_KOH - 0.0028*c_KOH**2 + 0.005332*c_KOH*T + 207.2 * c_KOH/T + 0.001043*c_KOH**3 - 0.0000003*c_KOH**2*T**2) # // in S/m Olivier eq. 15

    R_el_free = 1 / sgm_KOH_free * (pv.l_AN_S / pv.srf_AN + pv.l_CA_S / pv.srf_CA) # // in 1/S
    R_el_b    = R_el_free * ((1/((1 - 2 * av.theta / 3)**(3/2))) -1)

    #R_mem = (0.060 + 80*np.exp(T/50))/(10000*srf_SE*10**4)            #//henao
    R_mem = pv.rho_SE * pv.tau_SE**2 * pv.dlt_SE / (pv.ome_SE * pv.eps_SE * pv.srf_SE)             #//abdin


    dU_ohm = i * pv.srf_CA * (R_AN + R_CA + R_el_free + R_el_b + R_mem + av.R_sep) # // in V


    return(dU_ohm)


### auxilliary calculations
################################################################################

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
