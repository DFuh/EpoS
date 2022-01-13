'''
auxilliary calculations
flws
AEL
'''

#TODO: check current density unit in inputs/ methods
#TODO: "self" vs. "obj."

#TODO: check partial pressures !
#TODO: redundant rho-calc. (KOH) in clc_k_Li !

import numpy as np
import math


def clc_flws_auxpars(obj, T):
    '''
    calculate necessary auxilliary parameters for flw-balance
    '''
    #obj.av.d_mem = obj.pec.d0_mem
    #obj.av.D_eff_H2, obj.av.D_eff_O2 = clc_diffusion_coefficient(obj, obj.pec, T)
    ### H2O
    # obj.av.rho_H2O      = clc_rho_H2O(T)
    # obj.av.rho_ely = clc_rho_KOH(obj, T, obj.pec.w_KOH)
    obj.av.viscos_H2O   = 1 / (0.1 * T**2 - 34.335 * T + 2472) # Dynamic viscosity of water // in Pa s | wikipedia
    obj.av.viscos_KOH = viscos_KOH(obj, T)
    obj.av.pp_H2O = clc_pp_H2O(obj, obj.pec, T, )
    ###

    return

def clc_cp_H2O(obj, T):
    '''
    adopted from Petipas 2013, eq. 10
    and check:
    https://www.engineeringtoolbox.com/molarity-molality-mole_fraction-weight_percent-grams_solute_liter_solution-conversion-formula-calculator-d_1998.html

    '''
    return 92,67*1,1e1*T + 1.78e-4*T**2

def clc_cp_KOH():
    '''
    see Hnedkovsky2017
    '''
    return

def testcalc(x,y):
    return x+y*2

def clc_Deff():
    '''
    calc. effective binary diffusion cioefficient
    '''

    return


def clc_conc():
    '''
    calc. concentrations in respective interfaces in compartments
    -> separator/electrolyte
    -> electrolde/electrolyte

    materialbalance: concentrations in ...?
    '''

    return

'''
================================================================================
 all following methods adopted from:
 /home/dafu_res/0_modeling/py_scripts_v101/mod_3/flws/
 parameter_sensibility_bsc_matbal_v1002.py
'''

def matbal_preclc(self, T,i, ):
    '''
    input unit of current density: A/m²
    expected unit in functions: A/cm²
    see line 86 !
    '''
    ### comp. HAug 2017 table 5 !!! -> implement
    ''' EDIT: V constant (no dependence on i)'''

    i = i*1e-4 # TODO: remove factorization !

    # clc delta-p : "results from the concave curvature of the spherically shaped gas bubbles" (Haug)
    d_b_an, d_b_ca  = self.av.d_b
    dp_an           = 4 * self.av.gamma / d_b_an
    dp_ca           = 4 * self.av.gamma / d_b_ca

    ### pressures of gaseous phase (???) in respective compartment
    self.av.pp_an = self.p.anode + dp_an
    self.av.pp_ca = self.p.cathode + dp_ca

    # clc volume of gaseous species in resp. compartment
    eps_an, eps_ca = self.av.epsilon
    self.av.V_gas_an =  self.av.Vhcell * (eps_an *(self.p.anode/ (self.p.anode + dp_an)  ))
    self.av.V_gas_ca =  self.av.Vhcell * (eps_ca *(self.p.cathode/ (self.p.cathode + dp_ca) ))

    # clc generation of species (molar flows)
    # n_gn_H2 = (i * 1e4 * self.pcll.active_cell_area) / (2*self.pec.F) # H2 generation; cathode // in mol/s
    # n_gn_O2 = (i * 1e4 * self.pcll.active_cell_area) / (4*self.pec.F) # O2 generation; anode // in mol/s
    n_gn_H2 = (i * 1e4 * self.pcll.active_cell_area_matbal) / (2*self.pec.F) # H2 generation; cathode // in mol/s
    n_gn_O2 = (i * 1e4 * self.pcll.active_cell_area_matbal) / (4*self.pec.F) # O2 generation; anode // in mol/s
    self.av.n_gn = n_gn_H2, n_gn_O2

    return


def matbal(self,y_in, t, T, i, n_in, prm_dff=False, prm_drc=False): # i-input in A/cm²

    c_H2_an, c_H2_ca, c_O2_an, c_O2_ca = y_in#[:4] #c_in
    n_H2_an, n_H2_ca, n_O2_an, n_O2_ca = n_in#[4:] #n_in

    ''' calculation of material balance, based on HAUG2017;
    edit: permeation??? '''
    #print('--matbal -> kLi: ', self.k_Li)
    #print('--matbal -> n_in: ', n_in)
    kL_H2_an        = self.av.k_Li[0] # mass transfer-coeff.; anode // in m/s
    kL_H2_ca        = self.av.k_Li[1] # mass transfer-coeff.; cathode // in m/s
    kL_O2_an        = self.av.k_Li[2] # mass transfer-coeff.; anode // in m/s
    kL_O2_ca        = self.av.k_Li[3] # mass transfer-coeff.; cathode // in m/s

    # H2 mixer outlet concentration // in mol/m³
    c_mix_H2        = ((self.av.VL_an * c_H2_an) +
                        (self.av.VL_ca * c_H2_ca)) / (2*self.av.VL)
    # O2 mixer outlet concentration // in mol/m³
    c_mix_O2        = ((self.av.VL_an * c_O2_an) +
                        (self.av.VL_ca * c_O2_ca)) / (2*self.av.VL)

    ### partial pressures
    pp_H2_an, pp_H2_ca, pp_O2_an, pp_O2_ca= clc_partialpressures(self, n_in)
    #print('--matbal -> pps: ', pp_H2_an, pp_H2_ca, pp_O2_an, pp_O2_ca)
    ### clc equilibrium concentrations
    S_H2, S_O2 = self.av.S_ik # // in mol/ m3 Pa
    c_eq_H2_an      = pp_H2_an * S_H2 # equilibrium concentration; anode
    c_eq_H2_ca      = pp_H2_ca * S_H2 # equilibrium concentration; cathode
    c_eq_O2_an      = pp_O2_an * S_O2 # equilibrium concentration; anode
    c_eq_O2_ca      = pp_O2_ca * S_O2 # equilibrium concentration; cathode

    #c_sat_an = 101325 * S_O2
    #c_sat_ca = 101325 * S_H2

    D_H2_KOH = self.av.D_ik[0]

    # permeation flux density // in mol/m² s
    N_perm_H2       = 0 #par.clc_perm_AEL(T,self.p_ca,(c_H2_ca - c_H2_an),
                        #                D_H2_KOH, diff=prm_dff, darc=prm_drc)
    N_perm_O2       = 0 #???            # permeation flux density // in mol/m² s

    ### Volume of liquid species in respective compartment
    V_liq_an = self.av.Vhcell - self.av.V_gas_an
    V_liq_ca = self.av.Vhcell - self.av.V_gas_ca

    A_GL_an, A_GL_ca    = self.av.A_GL
    fG_H2, fG_O2        = self.av.f_G
    N_gn_H2, N_gn_O2    = self.av.n_gn # // in mol/s
    A_sep = self.pcll.separator_area
    ### Calc.concentrations | Haug2017, eq. 4,5
    # H2 in anodic compartment
    c_H2_an = ( (1/(self.av.VL_an + A_GL_an * kL_H2_an))
                * (self.av.VL_an * c_mix_H2 + A_GL_an * kL_H2_an * c_eq_H2_an
                                                    + N_perm_H2 *A_sep ) )
    #dc_H2_an = (self.VL_an/V_liq_an) * c_H2_an # change in conc. ?
    #n_H2_an = (- A_GL_an * kL_H2_an * (c_eq_H2_an - c_H2_an)) # molar flow

    # O2 in anodic compartment
    c_O2_an = ( (1/(self.av.VL_an + A_GL_an * kL_O2_an))
                * (self.av.VL_an * c_mix_O2 + A_GL_an * kL_O2_an * c_eq_O2_an
                            + N_perm_O2 *A_sep + (1 - fG_O2)*N_gn_O2) )
    #dc_O2_an = (self.VL_an/V_liq_an) * c_O2_an
    # ??? c_O2_an = c_sat_an#c_eq_O2_an #c_sat_an
    #n_O2_an = - A_GL_an * kL_O2_an * (c_eq_O2_an - c_O2_an) + fG_O2 * N_gn_O2

    # H2 in cathodic compartment
    c_H2_ca = ( (1/(self.av.VL_ca + A_GL_ca * kL_H2_ca))
                * (self.av.VL_ca * c_mix_H2 + A_GL_ca * kL_H2_ca * c_eq_H2_ca
                            - N_perm_H2 *A_sep + (1 - fG_H2)*N_gn_H2) )
    #dc_H2_ca = (self.VL_ca/V_liq_ca) * c_H2_ca
    # ??? c_H2_ca = c_sat_ca #c_eq_H2_ca #c_sat_ca
    #n_H2_ca = - A_GL_ca * kL_H2_ca * (c_eq_H2_ca - c_H2_ca) + fG_H2 * N_gn_H2


    # O2 in cathodic compartment
    c_O2_ca = ( (1/(self.av.VL_ca + A_GL_ca * kL_O2_ca))
                * (self.av.VL_ca * c_mix_O2 + A_GL_ca * kL_O2_ca * c_eq_O2_ca
                                                    - N_perm_O2 *A_sep) )
    #dc_O2_ca = (self.VL_ca/V_liq_ca) * c_O2_ca
    #n_O2_ca = - A_GL_ca * kL_O2_ca * (c_eq_O2_ca - c_O2_ca)


    c_out = c_H2_an, c_H2_ca, c_O2_an, c_O2_ca #
    #n_out = n_H2_an, n_H2_ca, n_O2_an, n_O2_ca #

    return c_out #, n_out#, pp


def clc_partialpressures(self, n_in):
    '''
    Calc partial pressure of gaseous species in respective compartment
    ---
    returns
        pp_
        ...
    '''
    # TODO: check claculation!!!

    n_H2_an, n_H2_ca, n_O2_an, n_O2_ca = n_in#[4:] #n_in
    #print('n_H2_an: ', n_H2_an)
    #print('---matbal -> pp_an, pp_ca: ', self.pp_an, self.pp_ca)
    #print('---matbal -> pp_H2O: ', self.pp_H2O)
    if n_H2_ca == 0 : #???
    #if N_gn_H2 == 0:
        pp_H2_an = 0
        pp_H2_ca = 0#(pp_ca - p_H2O)
        pp_O2_an = 0#(pp_an - p_H2O)
        pp_O2_ca = 0
    else:
        pp_H2_an = ( n_H2_an / (n_H2_an + n_O2_an)) * (self.av.pp_an - self.av.pp_H2O)
        pp_H2_ca = ( n_H2_ca / (n_H2_ca + n_O2_ca)) * (self.av.pp_ca - self.av.pp_H2O)
        pp_O2_an = ( n_O2_an / (n_O2_an + n_H2_an)) * (self.av.pp_an - self.av.pp_H2O)
        pp_O2_ca = ( n_O2_ca / (n_O2_ca + n_H2_ca)) * (self.av.pp_ca - self.av.pp_H2O)

    return pp_H2_an, pp_H2_ca, pp_O2_an, pp_O2_ca


def clc_molarflows(self, c_in, n_in):#, mb_out):
    '''
    clacl molar flows from concentrations
    '''
    [c_H2_an, c_H2_ca, c_O2_an, c_O2_ca] = c_in
    n_H2_an, n_H2_ca, n_O2_an, n_O2_ca = n_in

    #[A_sep, [VL_an,  VL_ca, VL], [V_hcell, V_gas_an, V_gas_ca],
    # [A_GL_an, A_GL_ca, dp], k_L, [p_H2O, pp_an, pp_ca], [fG_H2, fG_O2], [N_gn_H2, N_gn_O2], [S_H2, S_O2]] = mb_out

    kL_H2_an        = self.av.k_Li[0] #*tfac# mass transfer-coeff.; anode // in m/s
    kL_H2_ca        = self.av.k_Li[1] #*tfac# mass transfer-coeff.; cathode // in m/s
    kL_O2_an        = self.av.k_Li[2] #*tfac# mass transfer-coeff.; anode // in m/s
    kL_O2_ca        = self.av.k_Li[3] #*tfac# mass transfer-coeff.; cathode // in m/s

    ### partial pressures
    pp_H2_an, pp_H2_ca, pp_O2_an, pp_O2_ca= pp =clc_partialpressures(self, n_in)

    S_H2, S_O2 = self.av.S_ik # // in mol/ m3 Pa

    c_eq_H2_an      = pp_H2_an * S_H2 #*testfactor # equilibrium concentration; anode
    c_eq_H2_ca      = pp_H2_ca * S_H2 #*testfactor # equilibrium concentration; cathode
    c_eq_O2_an      = pp_O2_an * S_O2 # equilibrium concentration; anode
    c_eq_O2_ca      = pp_O2_ca * S_O2 # equilibrium concentration; cathode

    A_GL_an, A_GL_ca    = self.av.A_GL
    fG_H2, fG_O2        = self.av.f_G
    N_gn_H2, N_gn_O2    = self.av.n_gn

    n_H2_an = (- A_GL_an * kL_H2_an * (c_eq_H2_an - c_H2_an))
    n_O2_an = - A_GL_an * kL_O2_an * (c_eq_O2_an - c_O2_an) + fG_O2 * N_gn_O2 # ??? -> check, if valid updated
    n_H2_ca = - A_GL_ca * kL_H2_ca * (c_eq_H2_ca - c_H2_ca) + fG_H2 * N_gn_H2
    n_O2_ca = - A_GL_ca * kL_O2_ca * (c_eq_O2_ca - c_O2_ca)

    scl = self.pcll.active_cell_area/self.pcll.active_cell_area_matbal
    n_out = n_H2_an*scl, n_H2_ca*scl, n_O2_an*scl, n_O2_ca*scl
    pp_out = pp_H2_an, pp_H2_ca, pp_O2_an, pp_O2_ca
    return n_out, pp_out



#===============================================================================
'''               ---   calc. Parameters ---    '''
''' adopted from
    /home/dafu_res/0_modeling/py_scripts_v101/mod_3/flws/
    parameter_sensibility_par_applied_fctrs_clc_v3011.py
'''
#===============================================================================

def clc_D_ik(obj, T, w_KOH):

    def Dfit(t, a,b,c):
        F = a*pow(t,2)+b*t+c
        return F

    '''
    curveFit based on values from Haug2017 ExpOP
    '''
    popt_H2 = np.array([ 9.40624977e-13, -5.46323734e-10,  8.12858904e-08])
    popt_O2 = np.array([ 3.73124956e-13, -2.22938720e-10,  3.40223904e-08])
    dH2 = Dfit(T, *popt_H2) #* (1 + (sns * (obj.snsfctr.fD_ik_H2 -1)))
    dO2 = Dfit(T, *popt_O2) #* (1 + (sns * (obj.snsfctr.fD_ik_O2 -1)))
    #return dH2, dO2
    return (dH2, dO2)


def clc_k_Li(obj,T,i, d_b=None, D_ik=None, epsilon=None, fctr=1, testprint=False):
    '''
    alternative calc. of kLi, detailed version (based on Haug)

    returns mass transfer coefficients for H2 and O2 at an/ca
    '''

    g = 9.81 # // in m/s²
    rho_L               = clc_rho_KOH(obj, T, obj.pec.w_KOH)
    eta_L               = viscos_KOH(obj, T) # KOH viscosity Pa s
    #print('eta_L: ', eta_L)
    eps_g_an, eps_g_ca  = obj.av.epsilon #clc_epsilon(obj, i)
    d_b_an, d_b_ca      = obj.av.d_b #clc_d_b(obj, i, i_crits, gamma, beta, rho_H2, rho_KOH)
    #print('d_b: ', d_b_an, d_b_ca)
    D_H2_KOH, D_O2_KOH = obj.av.D_ik

    u_b_an = 0.33*g**(0.76) * (rho_L / eta_L)**0.52 *(d_b_an/2)**1.28
    u_b_ca = 0.33*g**(0.76) * (rho_L / eta_L)**0.52 *(d_b_ca/2)**1.28

    u_sw_an = ((u_b_an * (1 / (1+(eps_g_an / (1 - eps_g_an)**2))) )
                * ( (1-eps_g_an) / (1 +( (1.05) / ((1 + ( 0.0685 / (eps_g_an**2)) )**0.5 -0.5) )) ))
    u_sw_ca = ((u_b_ca * (1 / (1+(eps_g_ca / (1 - eps_g_ca)**2))) )*
            ( (1-eps_g_ca) / (1 +( (1.05) / ((1 + ( 0.0685 / (eps_g_ca**2)) )**0.5 -0.5) )) ))

    Re_an = (rho_L * d_b_an * u_sw_an) / (eta_L )
    Re_ca = (rho_L * d_b_ca * u_sw_ca) / (eta_L )

    Sc_O2 = eta_L / (rho_L * D_O2_KOH)
    Sc_H2 = eta_L / (rho_L * D_H2_KOH)

    #print('Re an,ca: ', Re_an, Re_ca)
    #print('Sc H2,O2: ', Sc_H2, Sc_O2)
    #fctr=0.38 #works good, if db NOT calc. according to HAug
    # fctr = 0.4 without permeation, dbeq=0
    #fctr=0.55
    Sh_H2_an = (2 + ( (0.651 *(Re_an * Sc_H2)**1.72 )
                        / (1 + (Re_an * Sc_H2)**1.22) )*fctr)
    Sh_O2_ca = (2 + ( (0.651 *(Re_ca * Sc_O2)**1.72 )
                        / (1 + (Re_ca * Sc_O2)**1.22) )*fctr)
    Sh_O2_an = (2 + ( (0.651 *(Re_an * Sc_O2)**1.72 )
                        / (1 + (Re_an * Sc_O2)**1.22) )*fctr)
    Sh_H2_ca = (2 + ( (0.651 *(Re_ca * Sc_H2)**1.72 )
                        / (1 + (Re_ca * Sc_H2)**1.22) )*fctr)

    Sh_H2_an_den = ( (1 + (Re_an * Sc_H2)**1.22) )
    Sh_O2_ca_den = ( (1 + (Re_ca * Sc_O2)**1.22) )
    Sh_O2_an_den = ( (1 + (Re_an * Sc_O2)**1.22) )
    Sh_H2_ca_den = ( (1 + (Re_ca * Sc_H2)**1.22) )

    if testprint:
        print('Sh_H2_ca : ', Sh_H2_ca)
        print('Sh_O2_ca : ', Sh_O2_ca)
        print('Sh_H2_an : ', Sh_H2_an)
        print('Sh_O2_an : ', Sh_O2_an)

    kL_H2_ca = Sh_H2_ca * D_H2_KOH / d_b_ca
    kL_O2_ca = Sh_O2_ca * D_O2_KOH / d_b_ca
    kL_H2_an = Sh_H2_an * D_H2_KOH / d_b_an
    kL_O2_an = Sh_O2_an * D_O2_KOH / d_b_an

    db = [d_b_an, d_b_ca, d_b_an, d_b_ca]
    D = [D_H2_KOH, D_H2_KOH, D_O2_KOH, D_O2_KOH]
    Sh = [Sh_H2_an, Sh_H2_ca, Sh_O2_an, Sh_O2_ca]
    #Sh_num = [Sh_H2_an_num, Sh_H2_ca_num, Sh_O2_an_num, Sh_O2_ca_num]
    Sh_den = [Sh_H2_an_den, Sh_H2_ca_den, Sh_O2_an_den, Sh_O2_ca_den]

    kL = [kL_H2_an, kL_H2_ca, kL_O2_an, kL_O2_ca]

    #TODO: check following 12 lines
    kL_n = []
    for n,kli in enumerate(kL):
        if math.isnan(kli):
            kL_n.append(0)
            print(f'kli isnan: [n]: {n}, i: {i}, Sh_den: {Sh_den[n]}')
            print('Re an,ca: ', Re_an, Re_ca)
        else:
            kL_n.append(kli)

    if testprint:
        return tuple(kL_n), tuple(Sh)
    else:
        return tuple(kL_n) #(kL_H2_an, kL_H2_ca, kL_O2_an, kL_O2_ca)


def viscos_KOH(obj, T):
    '''
    dyn. viscosity of KOH
    Haug 2017, appendix ref.:[73]
    '''
    # coefficients | Haug
    cm = np.array([0.9105535967, -0.01062211683, 4.680761561e-5,
                    -9.209312883e-8, 6.814919843e-11])
    eta = 0
    for n in range(len(cm)):
        eta += cm[n]*T**n
    return eta

def clc_Vhcell(obj, T,w_KOH,lx=0.16, ly=0.015, lz=0.145,fctr=1):
    #V_hcell = 00.16 * 0.015 * 0.145
    # geometrical dimensions of half-cell
    # height x width x depth) // in m 3 | Haug2017_mod
    Vhcell = lx * ly * lz*fctr
    return Vhcell

def clc_gamma(obj, T, w_KOH):
    #def surfacetension(conc_in, temp_in):
    '''
    T     // in K
    w_KOH       // in wt%

    zeta - concentration in wt% (valid: 0.02 ... 0.58)
    theta - Temp in °C

    returns
    surface tension // in N/m
    '''
    theta = T - 273.15 # Conversion from K to °C
    if w_KOH >1:
        zeta = w_KOH /100
    else:
        zeta = w_KOH

    a = np.array([[75.4787, -0.138489, -0.336392e-3, 0.475362e-6, -0.264479e-9 ],
                 [-32.8890, 1.34382, -0.910138e-2, 0.396124e-4, -0.573565e-7],
                 [614.527, -12.8736, 0.104855, -0.449076e-3, 0.651193e-6],
                 [-1455.06, 39.851, -0.344234, 0.144383e-2, -0.207599e-5],
                 [1333.62, -38.3316, 0.335129, -0.137313e-2, 0.194911e-5]])


    #def clc_sT(zeta, theta):
    '''
    calculate surface tension
    w.r.t
    temperature (theta) in °C
    fraction of KOH in solution (zeta) in kg/kg

    according to Feldamp 1969 eq. 3
    single result for given zeta and theta
    '''

    outer = 0
    inner = 0
    for j in range(5):                      # range of i values in a_ik
        for k in range(5):                  # range o k values in a_ik

            inner += a[j,k] * theta**(k)    # inner sum
        outer += inner * zeta**(j) * 1e-3   # outer sum
        inner=0
    #return outer
    #return outer # return surface tension in
    return outer


def clc_epsilon(obj,T, i, fctr=1):
    """
    A short description.

    A bit longer description.

    Args:
        variable (type): description

    Returns:
        type: description

    Raises:
        Exception: description

    """

    #def Gas_VOID(i): #i-input in A/cm²
    '''
    calculation of gas voidage; equation and parameters adopted from Haug2017_mod,
    experimental data !
    --- > heavily dependent on cell geometry !!!
    i-input in A/cm² ??
    ---------
    returns
    gas_voidage in anodic and cathodic compartment // in 1
    '''
    ''' valid @ i>= 0.3kA/m²'''

    #fctr=1 #0.48 # good fit:0.49

    ie = i*1e1 # -> ie in kA/m²
    if ie >= getattr(obj, 'i_crit_eps',0.03):
        X1 = np.array([0.59438, 0.76764])
        X2 = np.array([0.59231, 0.73233])
        X3 = np.array([0.75647, 0.73457])
        eps_an,eps_cat = (X1 - (X2 * X3 **(ie))) #* fctr
    else:
        # print('clc_epsilon: i-crit !)')
        eps_an, eps_cat = 0.2, 0.2
    #print('eps: ', eps_an, eps_cat)
    return (eps_an, eps_cat*fctr)


def clc_A_GL(obj, T, i):
    Vhcell = obj.av.Vhcell
    epsilon = obj.av.epsilon
    gamma = obj.av.gamma
    d_b = obj.av.d_b
    p_an = obj.p.anode #p_an,
    p_ca = obj.p.cathode #p_ca
    #def AGL(T,i): # i-input in A/cm²
    ''' calculation of gas-liquid interfacial area'''
    ''' an/ cat???'''
    '''
    p_an = 101325 # anode pressure // in Pa
    p_ca = 101325 # cathode pressure // in Pa
    '''
    #TODO: pan, p_ca ---> global variables
    #V_hcell = 00.16 * 0.015 * 0.145
    # geometrical dimensions of half-cell
    # (height x width x depth) // in m 3 | Haug2017_mod

    #eps_an, eps_ca = Gas_VOID(i) # i-input in A/cm² /// out: eps_an, eps_cat
    eps_an, eps_ca = epsilon
    #gam = surfacetension(w_KOH,T)
    #gam = gamma(T) #


    #d_b_an  = 100 * 1e-6 # 100mikrometer // in m |Haug
    #d_b_ca = gas_bubb0(T,i,gam, haugeq=True)
    d_b_an, d_b_ca = d_b
    #print('d_b (in AGL): ', d_b)

    dp_an       = 4 * gamma / d_b_an
    dp_ca      = 4 * gamma / d_b_ca

    dp = dp_an, dp_ca

    #po_an,po_ca = p_an,p_ca
    #po = po_an,po_ca
    #print('p_an,dp_an: ', p_an,dp_an)
    pp_an = (p_an + dp_an)
    pp_ca = (p_ca + dp_ca)
    #pp =  pp_an, pp_cat

    V_gas_an =  Vhcell * (eps_an *(p_an/pp_an))
    V_gas_ca =  Vhcell * (eps_ca *(p_ca/pp_ca))

    V_b_an = (np.pi/6) * d_b_an**3
    V_b_ca = (np.pi/6) * d_b_ca**3

    S_b_an = np.pi * d_b_an**2
    S_b_ca = np.pi * d_b_ca**2

    A_GL_an = V_gas_an / V_b_an * S_b_an
    A_GL_ca = V_gas_ca / V_b_ca * S_b_ca

    return (A_GL_an, A_GL_ca) #, dp, [V_hcell, V_gas_an, V_gas_ca]



def clc_d_b(obj,T, i, ):
    #i_crits = getattr(obj, 'i_crit_db',0.03)
    gamma = obj.av.gamma
    beta = obj.av.beta
    rho_G = obj.av.rho_H2
    rho_L = obj.av.rho_KOH
    #def gas_bubb0(T,i,*args,haugeq=False, i_crit=0.03): #i-input in A/cm²
    #gam = *gam[0]
    #print(args)
    #if len(args) >0:
    #    gam = args
    #else:
        #print('fill gam')
        #gam = gamma(T)
        #gam = surfacetension(w_KOH,T)

    #print(gam)


    '''vgl. Abb5 in Haug eq.33/34 ! // small current-values???'''
    #print('Line 320: edited current-density limit for gas-bubble diameter !')

    #if i >= obj.i_crit_db: #0.03: # A/cm²
    if i >= getattr(obj, 'i_crit_db',0.03):
        ie = i*1e4 # -> A/m²
        g = 9.81 # // in m/s²
        #rho_G = np.array([0.0749, 0.0727, 0.0706, 0.0686])
        # rho_G: density of hydrogen @ 1bar 50...80°C

        #rho_L = rhoL(T)
        #gam = gamma()[3]
        #beta = 0.125*0.675 #??? ankle...? #edit (2020-08): factor 0.675
        #rho_G = rho_H2
        #rho_L = rho_KOH


        #----------
        #if obj.db_eqh:
        if getattr(obj, 'db_eqh',False):
            #print('-------------------> db_clc_Haug')
            d_b = 593.84 * 1e-6 * (1 + 0.2*i*1e4)**(-0.25) # HAugs calc
        else: #--------------------------------
            #print('-------------------> db_clc_Lit')
            d_b0 = 1.2 * beta * np.sqrt(gamma / (g * (rho_L - rho_G)))
            d_b = d_b0 * (1 + 0.2*ie)**(-0.45)
    else:
        d_b = (200*1e-6)#*1.25

    #print('-->i: ', i, 'd_b: ', d_b)

    d_b_an = 100 *1e-6
    d_b_ca = d_b
    #d_b = 593.84 * 1e-6 * (1 + 0.2*i*1e4)**(-0.25) # HAugs calc
    #print('d_b Haug calc.: ',d_b )
    #d_b = 187*1e-6
    #print('d_b_an, d_b_ca: ', d_b_an, d_b_ca)
    return (d_b_an, d_b_ca)


def clc_beta(*args, fctr=1, sns=False):
    # --> REALLY good: beta = 0.825*np.pi/2
    #1.3#0.675 #1.5 #0.125*0.675 #??? ankle...? #edit (2020-08): factor 0.675
    #print('beta-fctr: ', fctr)
    beta = fctr * np.pi/2 #* (1 + (sns * (obj.snsfctr.fbeta -1)))
    return beta


def clc_f_G(obj, T, i, fctr=3):
    #def f_G(i,n): # i-input: A/cm2
    ''' calculation of hydrogen gas evolution efficiency'''
    ie = i*1e4
    #f_G = np.zeros(n,i)
    '''flow rates of ely also relevant!!!'''

    # f_G1: gas evolution efficiency // in
    # | Vogt (Haug[43]) ||0.5M H2SO4 25°C
    f_G1 = (1 - 1.35 * ie**(-0.095))
    # f_G2: gas evolution efficiency // in
    # | Pierre/Chin Kwie Joe (Haug[44]) ||1M KOH 25°C
    f_G2 = 0.5653 * (1  - np.exp(-0.002061*ie))
    # f_G3: gas evolution efficiency // in
    # | Chin Kwie Joe (Haug[44]) || 6.8M NaOH 80°C)
    f_G3 = 0.16302 *ie**(0.18527)
    # f_G4: gas evolution efficiency // in | Haug
    f_G4 = 0.25744 * ie**(0.14134)

    f_G = f_G1,f_G2,f_G3,f_G4

    if isinstance(fctr,int):
        f_G_ret = f_G[fctr]
    elif isinstance(fctr,tuple):
        f_G_ret = f_G[fctr[0]]*fctr[1]
    else:
        f_G_ret = f_G4 * fctr
    return (f_G_ret, 1)#0.909) #orig:(f_G_ret, 1) # O2->1 #[n]


def clc_S_ik(obj,T,w_KOH):
    #def solub(T):
    ''' solubility in KOH fitted on Fig.2 in Haug2017_exp
    retruns
    solubiities in mol/m3 Pa ???
    '''

    snsO2 = 1#0.8
    snsH2 = 1#0.8

    a_H2 = 5.386916952336441*1e-11
    b_H2 = 2655.181825221195
    c_H2 = 2.2337904127168284*1e-09

    S_H2 = snsH2 * a_H2 * np.exp(b_H2/T) + (c_H2*T)

    a_O2 = 7.667987754500667*1e-10
    b_O2 = 1890.1434701686374
    c_O2 = 1.8224579596926922*1e-09

    S_O2    =snsO2 * a_O2 * np.exp(b_O2/T) + (c_O2*T)
    pas     = 101325 # ref pressure in Pa
    return S_H2/pas, S_O2/pas

def clc_rho_H2(obj,T,w_KOH):
    m_rho_G = -0.00021000000000000017
    b_rho_G = 0.14271150000000946
    rho_G = m_rho_G * T + b_rho_G # density of pure hydrogen
    return rho_G

'''
-> clc_rho_KOH moved to clc_auxvals in ael_aux_vXX
'''
def clc_rho_KOH(obj,T, w_KOH):
   ''' calc density of aqueous KOH-solution'''
   '''valid: 0.01 ... 200 °C /// w = 0 ... 0.5 *100 wt% KOH'''
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
   rho_L_out = rho_int * np.exp(0.86*w_KOH)
   return rho_L_out # in kg/m³

def clc_rho_H2O(T):
    return 999.972 - 7*10**(-3)*(T-273.15-20)**2 # Source=? factors "20" vs "4" ??? (both not valid)

def clc_pp_H2O(obj,pec, T ):
    #w_KOH = 0.31
    #M_KOH = 56.11 *1e-3
    #print('T: ', T)
    #print('w_KOH: ', w_KOH)
    # water vapour part pressure in KOH model J.BALEJ 1985 bis 300°C (eq. 6)  und b_KOH 18
    # check for eq.7 also

    b_KOH     = pec.w_KOH / ((1-pec.w_KOH) * pec.M_KOH)  # MOLALITY KOH at 30w% KOH solution check -> in mol/kg Massenkonzentration (wikipedia)

    p_H2O_KOH = 10** (- 0.01508 * b_KOH - 0.0016788 * b_KOH**2
                + 2.25887*10**(-5)*b_KOH**3
                +    (1 - 0.0012062*b_KOH
                + 5.6024*10**(-4)*b_KOH**2
                - 7.8228*10**(-6)*b_KOH**3) * (35.4462 - 3343.93 / T
                - 10.9 * np.log10(T) + 0.0041645 * T))
    # HAug: 27583 Pa @ 80°C

    return p_H2O_KOH
