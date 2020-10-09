'''
PEM
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
    dE_ca = pec.dG_ca / (2 * pec.F)
    dE_an = pec.dG_an / (2 * pec.F)

    ### Nernst voltage
    '''
    deviations from temp., pressure (conc.)
    '''
    dE_N_ca = pec.R * T / (2 * pec.F) * np.log(pp_H2_ca /pec.p0_ref )
    dE_N_an = pec.R * T / (2 * pec.F) * np.log((pp_O2_an /pec.p0_ref)**(1/2))

    dE0_ca = None
    dE0_an = None
    return


def ov_act(obj, pec, T, i):
    '''
    calculate activation overvoltage
    '''
    # TODO: implement dRct (charge transfer resistance (degr))
    ### exchange current density
    i0_ca     = 2 * pec.F * pec.k0_ca * T * np.exp( (- pec.Ae_ca / (pec.R*T) )) # Chandesris2014 // in A/m²
    i0_an     = 2 * pec.F * pec.k0_an * T * np.exp( (- pec.Ae_an / (pec.R*T) )) # Chandesris2014 // in A/m²

    if i >0:
        dU_act_ca  = 0#(R * T / (alph_cat * F *z) ) * np.log( i / ( i0_cat * rug_c * corrf * 1) )

        dU_act_an   = (pec.R * T / (pec.alpha_an * pec.F * 2) ) * np.log( i / ( i0_an * pec.rugos_an * obj.av.dRct) ) # arsinh ???
    #print('i:',i,'dU_act_An:',dU_Act_an)
    else:
        U_act_ca  = 0
        U_act_an  = 0
    return dU_act_ca, dU_act_an


def ov_ohm(obj, pec, T, i):

    ### Resistance of membrane
    obj.lambda_mem = clc_lambda_mem() # See: Ito et al2011 ! (eq.3-6 and table 1)
    sigma_mem   = obj.av.corr_dgr *  ((0.005139 * obj.av.lambda_mem) - 0.00326 ) * np.exp( 1268 * ( (1 / 303) - (1 / T) )) # ionic conductivity of membrane // in S/m | Springer1991: 1/(ohm*cm) , Olivier2017, Chandesris, Tjarks
    R_mem_c     = ( (obj.av.d_mem / (sigma_mem)) )
    ### Resistance of current collector
    R_cc_na      = (pv.d_Acc  / pv.sig_Acc)
    R_cc_ca     = (pv.d_Ccc / pv.sig_Ccc )                               # current collector resistance | in (ohm * m²)

    du_ohm = (R_mem_c + R_cc_an + R_cc_ca ) * i  # |in V (ohm*m² * A/m² ////// old:A/cm² * 10000cm²/m²)
    return du_ohm
