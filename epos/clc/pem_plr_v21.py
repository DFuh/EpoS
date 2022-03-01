'''
PEM
calculation: polarisation characteristics

act_overpot: (changed chandesris eq. to marangio (v20 -> v21)
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


def voltage_cell(obj, pec, T,i,p, pp=None, ini=False): #, A_cell=None):
    '''
    calculate cell voltage
    - Gibbs
    - Nernst
    - Activation overvoltage (changed chandesris eq. to marangio (v20 -> v21))
    - dV due to Ohmic losses
    '''
    p_ca, p_an = p.cathode, p.anode

    if pp is None:
        pp = obj.clc_m.flws.partial_pressure(obj, pec, T, p)

    #if not A_cell:
    #    A_cell = obj.pcll.active_cell_area
    ### Reversible cell voltage
    ((dE_rev_ca, dE_rev_an), dE_rev, U_tn)      = cv_rev(obj, pec, T, pp)
    #U_rev_ca = U_rev_[0]
    #U_rev_an = U_rev_[1]

    ### Activation overpotential
    U_act_ca, U_act_an = ov_act(obj, pec, T, i, ini=ini)

    ### Concentration overpotential
    U_conc_ca, U_conc_an = 0,0#ov_conc(obj, pec, T, i,apply_funct=True)
    # print('U_conc_ca, U_conc_an: ',U_conc_ca, U_conc_an)
    ### Additional Voltage due to Ohmic losses
    U_ohm = ov_ohm(obj, pec, T, i, ini=ini)


    U_ca = U_act_ca + U_conc_ca # dG_ca / (z*F) # cathodic halfcell potential
    U_an = U_act_an + U_conc_an # Anodic halfcell potential

    #U_rev, U_tn = None
    U_cell_0 = dE_rev + U_ca + U_an + U_ohm

    ### update t_uninterrupted
    if  U_cell_0 < obj.pec.treshold_voltage_interrupt:
        obj.av.t_uni = 0
    else:
        obj.av.t_uni += obj.av.t_diff

    # print(f'dE_rev: {dE_rev}, U_ca: {U_ca}, U_an: {U_an}, U_ohm: {U_ohm} ')

    (obj.av.dU_dgr_act,
    obj.av.dU_dgr_abs) = obj.clc_m.dgr.voltage_increase(obj, pec, T, i)
    # print('U_dgr_incr, U_dgr_abs: ',obj.av.dU_dgr_act, obj.av.dU_dgr_abs)

    U_cell = U_cell_0 + 0#obj.av.dU_dgr_act * obj.av.enable_dgr
    #U_dgr_abs = 0
    print('U_cell, dU_abs (plr) = ', U_cell, obj.av.dU_dgr_abs)

    #print(f'----> U_ca: {U_ca}  // U_an: {U_an}     // U_ohm: {U_ohm}   ///U_cell: {U_cell}')
    return (U_ca, U_an, obj.av.dU_dgr_act, U_cell)


def cv_rev(obj, pec, T, pp):
    '''
    calc reversible cell voltage
    '''
    '''
    @ standard pressure and temp.
    Schalenbach 2013 eq. 6 -11 (!)
    '''
    pp_H2_ca, pp_O2_an = pp[:2]
    #print('pp:', pp)
    ((dG_ca, dG_an), (dH_ca, dH_an)) = clc_gibbs_free_energy(obj, pec, T)
    dE0_ca = -dG_ca / (2 * pec.F)
    dE0_an = -dG_an / (2 * pec.F)
    #print(f'dE0_ca: {dE0_ca}    // dE0_an: {dE0_an}')

    U_tn = dH_an/(2*pec.F)
    ### Nernst voltage
    '''
    deviations from temp., pressure (conc.)
    '''
    if pp_H2_ca >0:
        dE_N_ca = pec.R * T / (2 * pec.F) * 0#np.log(pp_H2_ca /pec.p0_ref )
    else:
        dE_N_ca = 0
    if pp_O2_an >0:
        dE_N_an = pec.R * T / (2 * pec.F) * 0#np.log((pp_O2_an /pec.p0_ref)**(1/2))
    else:
        dE_N_an = 0

    dE_rev_ca = dE0_ca + dE_N_ca
    #print(f'dE_N_ca: {dE_N_ca}')
    dE_rev_an = dE0_an + dE_N_an
    #print(f'dE_N_an: {dE_N_an}')
    dE_rev = dE_rev_ca - dE_rev_an
    #print(f'dE_rev: {dE_rev}')
    return ((dE_rev_ca, dE_rev_an),dE_rev,U_tn)


def ov_act(obj, pec, T, i, ini=False, apply_funct=True):
    '''
    calculate activation overvoltage
    '''
    # TODO: implement dRct (charge transfer resistance (degr))

    ### exchange current density
    i0_ca     = 2 * pec.F * pec.k0_ca * T * np.exp( (- pec.Ae_ca / (pec.R*T) )) # Chandesris2014 // in A/m²
    i0_an     = 2 * pec.F * pec.k0_an * T * np.exp( (- pec.Ae_an / (pec.R*T) )) # Chandesris2014 // in A/m²
    #print('i0_an: ', i0_an)
    #print('i0_ca: ', i0_ca)
    if (i >0) & apply_funct:
        dU_act_ca  = (pec.R * T / (pec.alpha_ca * pec.F ) ) * np.arcsinh( i / ( 2*i0_ca)) # // MArangio eq. 12
        #dU_act_ca  = (pec.R * T / (pec.alpha_ca * pec.F *2) ) * np.log( i / ( i0_ca * pec.rugos_ca * 1) )
        dU_act_an  = (pec.R * T / (pec.alpha_an * pec.F ) ) * np.arcsinh( i / ( 2*i0_an)) # // MArangio eq. 11
        #if not ini:
        #    dU_act_an   = (pec.R * T / (pec.alpha_an * pec.F * 2) ) * np.log( i / ( i0_an * pec.rugos_an * obj.av.dRct) ) # arsinh ???
        #else:
        #    dU_act_an   = (pec.R * T / (pec.alpha_an * pec.F * 2) ) * np.log( i / ( i0_an * pec.rugos_an) ) # arsinh ???
    #print('i:',i,'dU_act_An:',dU_Act_an)
    else:
        dU_act_ca  = 0
        dU_act_an  = 0
    return dU_act_ca, dU_act_an

def ov_conc(obj, pec, T, i,apply_funct=True):

    if apply_funct:
        if obj.av.cnc_O2_mem_an !=0:
            dU_cnc_an = (pec.R*T)/(4*pec.F)*np.log(obj.av.cnc_O2_mem_an / pec.cnc_O2_mem_an_ref)
        else:
            dU_cnc_an = 0
        if obj.av.cnc_H2_mem_ca != 0:
            dU_cnc_ca = (pec.R*T)/(2*pec.F)*np.log(obj.av.cnc_H2_mem_ca / pec.cnc_H2_mem_ca_ref )

        else:
            dU_cnc_ca = 0
        # print(f'(plr,cnc) pec.cnc_H2_mem_ca_ref, obj.av.cnc_H2_mem_ca (@ i = {i})-> ', pec.cnc_H2_mem_ca_ref, obj.av.cnc_H2_mem_ca)
        ret = dU_cnc_ca, dU_cnc_an
    else:
        ret = 0,0
    return ret


def ov_ohm(obj, pec, T, i, ini=False):
    '''
    see: Marangio_2009
    - eq. 60 vs. 61 (simplified)
    '''
    #print(f'i: {i}')
    ### Resistance of membrane

    #if not ini:
    #    obj.av.lambda_mem = clc_lambda_mem() # See: Ito et al2011 ! (eq.3-6 and table 1)
    #    sigma_mem   = obj.av.corr_dgr *  ((0.005139 * obj.av.lambda_mem) - 0.00326 ) * np.exp( 1268 * ( (1 / 303) - (1 / T) )) # ionic conductivity of membrane // in S/m | Springer1991: 1/(ohm*cm) , Olivier2017, Chandesris, Tjarks
    #    R_mem_c     = ( (obj.av.d_mem / (sigma_mem)) )
    #if True: #else:
        # --> clc lambda --> |Yigit 2016 eq. 14 // Medina2010
        #lambda_mem = 16
    sigma_mem   = ((0.005139 * obj.av.lambda_mem) - 0.00326 ) * np.exp( 1268 * ( (1 / 303) - (1 / T) )) *1e2# ionic conductivity of membrane // in S/cm *1e2 | Springer1991: 1/(ohm*cm) , Olivier2017, Chandesris, Tjarks
    R_mem_c     = ( (pec.d0_mem / (sigma_mem)) )
    ### Resistance of current collector
    R_cc_an      = (pec.d_cc_an  / pec.sigma_cc_an)
    R_cc_ca     = (pec.d_cc_ca / pec.sigma_cc_ca )                               # current collector resistance | in (ohm * m²)
    #print(f'-R_mem: {R_mem_c}   // R_cc_ca: {R_cc_ca}   // R_cc_an: {R_cc_an} ')
    du_ohm = (R_mem_c + R_cc_an + R_cc_ca ) * i# |in V (ohm*m² * A/m² ////// old:A/cm² * 10000cm²/m²)
    return du_ohm

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

    # print(f'dG_an({T}) = {dG_an}')

    #dGa = -(-237.19*1e3)
    #U_rev = dGa / (2 * F)
    #U_tn = H_H2O /(2 * F) #- ???
    return ((dG_ca, dG_an), (dH_ca, dH_an))# // in
