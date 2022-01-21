'''
calc. basic values of plant
    -> rated power
    -> number of cells in stack
    -> PID params)
'''
# import plr
import math
import numpy as np
from collections import namedtuple
from importlib import import_module as impm
from dataclasses import dataclass

import epos.auxf.handlingdata as hd
import epos.auxf.faux as fx
from epos.clc import ctrl
#TODO: take into account power of peripherie!!!

@dataclass
class auxvals():
    pass
    #name: str# 'AuxVals'

def clc_pwr_vls(obj, bsc_par, par_dct):
    '''
    Derive all necessary values for plant setup (defining size and power)

    cases:
        0 - (P_N, P0_st_max), (u_N, u_max), A_cell,
            || calc: (i_N, i_max), (P_stack_N, P_stack_max), (n_cells_stack, n_cells_tot),

        1 - (,P0_st_max), (i_N), (u_N, u_max), A_cell
            || calc: (P_N, ), (,i_max), (), , (n_cells_stack, n_cells_tot)

        2 - (P_N, P0_stack_max), (), (i_N), A_cell
            ||

        3 - (P_N, P0_stack_max), (), (i_max), A_cell
            ||
    '''

    '''
    ----------------------------------------------------------------------------
    compact approach:
        loop trough respective dicts and set values
        reference dict needed (for non-mutable values)
    ----------------------------------------------------------------------------
    '''

    obj.clc_m = fx.ini_clc_versions(obj, bsc_par)
    scaling = bsc_par.get('thrml_scaling', False)

    bsc_par = bsc_par['bsc_par'] # Check and remove/correct!!
    obj.pplnt = hd.dct_to_nt(par_dct['plant'], subkey='value') # Plant parameters as namedtuple
    obj.pcll = hd.dct_to_nt(par_dct['cell'], subkey='value')  # Cell parameters as namedtuple
    obj.pop = hd.dct_to_nt(par_dct['operation'], subkey='value')     # Operation parameters as namedtuple
    obj.pec = hd.dct_to_nt(par_dct['electrochemistry'], subkey='value') # Electrochemistry parameters as namedtuple
    obj.bop = hd.dct_to_nt(par_dct['periphery'], subkey='value') # Periphery parameters as namedtuple
    obj.p = hd.dct_to_nt(par_dct['operation']['nominal_electrode_pressure'],
                                subkey='value')
    hd.ini_auxvals(obj, par_dct)

    iv = auxvals() # Ini dataclass in
    pv = auxvals() # Ini dataclass out
    #pv.name = 'AuxVals' !!! Not working properly !

    # print('====> Basic Par: ', bsc_par)
    #print('Dataclass: ', pv)
    #print('par_dct[cell]: ', par_dct['cell'])
    #TODO: belows code redundant (see namedtuples further below!)
    pv.T_N              = par_dct['cell']['temperature_nominal']['value']
    pv.T_max            = par_dct['cell']['temperature_max']['value']
    iv.i_N           = par_dct['cell']['current_density_nominal']['value']  # // in A/m²
    iv.i_ol          = par_dct['cell']['current_density_overload']['value']  # // in A/m²
    iv.i_max         = par_dct['cell']['current_density_max']['value']
    #if not iN:
    #    iN = imx
    iv.u_N           = par_dct['cell']['voltage_nominal']['value']
    iv.u_ol          = par_dct['cell']['voltage_overload']['value']
    iv.u_max         = par_dct['cell']['voltage_max']['value']
    iv.A_cell        = par_dct['cell']['active_cell_area']['value']
    iv.p_N = False
    iv.p_ol = False

    iv.P_plnt_max     = par_dct['plant']['power_of_plant_max']['value']
    iv.P_plnt_ol     = par_dct['plant']['power_of_plant_overload']['value']
    iv.P_plnt_N      = par_dct['plant']['power_of_plant_nominal']['value']
    if not iv.P_plnt_N:
        pv.P_plnt_N_target = bsc_par['rpow_el']
    else:
        pv.P_plnt_N_target = iv.P_plnt_N
    print('iv.P_plnt_N: ',iv.P_plnt_N)

    iv.P_stack_N     = par_dct['plant']['power_of_stack_nominal']['value']
    iv.P_stack_ol    = par_dct['plant']['power_of_stack_overload']['value']
    iv.P_stack_max  = par_dct['plant']['power_of_stack_max']['value']
    iv.pwr_frc_max       = par_dct['operation']['power_fraction_max']['value']
    iv.pwr_frc_min       = par_dct['operation']['power_fraction_min']['value']

    iv.n_clls_plnt   = par_dct['plant']['number_of_cells_in_plant_act']['value']
    iv.n_clls_st     = par_dct['plant']['number_of_cells_in_stack_act']['value']
    iv.n_clls_plnt_max = par_dct['plant']['number_of_cells_in_plant_max']['value']
    iv.n_clls_st_max = par_dct['plant']['number_of_cells_in_stack_max']['value']
    iv.n_st          = par_dct['plant']['number_of_stacks_act']['value']
    iv.n_st_max = par_dct['plant']['number_of_stacks_max']['value']

    iv.Ct_st_ref = par_dct['plant']['heat_capacity_st_ref']['value']
    iv.Ct_st = par_dct['plant']['heat_capacity_st']['value']
    iv.n_clls_st_ref = par_dct['plant']['number_of_cells_in_stack_ref']['value']
    iv.A_c_ref = par_dct['plant']['active_cell_area_ref']['value']
    iv.Rt_st_ref = par_dct['plant']['thermal_resistance_st_ref']['value']
    iv.Rt_st = par_dct['plant']['thermal_resistance_st']['value']
    iv.UA_hx0_st_ref = par_dct['plant']['UA_hx0_st_ref']['value']
    iv.UA_hx0_st = par_dct['plant']['UA_hx0_st']['value']
    iv.dm_clnt_ref = par_dct['plant']['dm_clnt_st_ref']['value']

    pv.prod_rate_N      = par_dct['plant']['flowrate_H2_nominal']['value']
    pv.unit_prod_rate_N      = par_dct['plant']['flowrate_H2_nominal']['unit']
    #pv.iv_P_rect = obj.bop['power_rectifier_nominal']['value']
    pv.P_pmp_ely = par_dct['periphery']['power_pump_ely_nominal']['value']
    pv.P_pmp_clnt = par_dct['periphery']['power_pump_coolant_nominal']['value']
    pv.cp_coolant = par_dct['periphery']['cp_coolant']['value']
    pv.cp_ely = par_dct['periphery']['cp_ely']['value']
    pv.dT_min_ely = par_dct['periphery']['dT_min_ely']['value']
    pv.dT_min_coolant = par_dct['periphery']['dT_min_coolant']['value']
    iv.dp_ely_cycle = par_dct['periphery']['dp_ely_cycle']['value']
    iv.dp_coolant_cycle = par_dct['periphery']['dp_coolant_cycle']['value']

    iv.dm_clnt_max = par_dct['periphery']['massflow_coolant_max']['value']

    iv.dV_ely_nom = par_dct['periphery']['volumetricflow_ely_nominal']['value']
    #pv.iv_dm_ely_max = par_dct['periphery']['massflow_ely_max']['value']
    iv.corrfctr_Ct_st = par_dct['plant']['corrfctr_Ct_st'].get('value',1)
    ############################################################################


    ####
    ### get values from dict
    #cell_dct = par_dct['cell']

    ### calc i_N (u_N, args=(T_N, p, pp, i_ini, ))
        # obj: A_cell
    # p = obj.p
    pp = pp = obj.clc_m.flws.partial_pressure(obj, obj.pec, pv.T_N, obj.p)
    obj.clc_m.aux.clc_auxvals(obj, pv.T_N)
    # pv.rho_H2O = obj.clc_m.flws.xflws.clc_rho_H2O(pv.T_N)
    pv.rho_H2O = obj.av.rho_H2O
    #print('tec_el (?): ',obj.sup_par['tec_el'])
    #print('object: ', obj.__dict__)
    if bsc_par['tec_el'].lower() == 'ael':
        # pv.rho_ely = obj.clc_m.flws.xflws.clc_rho_KOH(obj,pv.T_N, obj.pec.w_KOH)
        pv.rho_ely = obj.av.rho_ely
    elif bsc_par['tec_el'].lower() == 'pem':
        #pv.rho_ely = obj.clc_m.flws.xflws.clc_rho_H2O(pv.T_N)
        pv.rho_ely = pv.rho_H2O
    else:
        print('No valid tech')
    ### clc u_rev and u_tn @ nominal temperature
    pv.dE_rev, pv.U_tn = obj.clc_m.plr.cv_rev(obj, obj.pec, pv.T_N, pp)[1:]
    # =========================================================================

    #print('CAUTION: chekc line 130 in bsc')
    #pv.dE_rev = abs(pv.dE_rev)
    ############################################################################
    ### target
    fctrs_flw = {'kg/s': 1/obj.pec.M_H2,
                'kg/h': 1/obj.pec.M_H2/3600}

    fctrs_flw_e = {'kg/s': 1/3600*39.4, # to kW
                'kg/h': 39.4}
    fctr = fctrs_flw.get(pv.unit_prod_rate_N, False)
    fctr_e = fctrs_flw_e.get(pv.unit_prod_rate_N, False)

    pv.i_target_flw = 2*obj.pec.F  * fctr* pv.prod_rate_N/ (iv.A_cell) # // in A/m2
    # pv.P_cell_target =
    pv.P_plnt_N_target_flw = fctr_e* pv.prod_rate_N

    ############################################################################

    ### limits

    ## cell

    pv.i_N_lim = minnz([iv.i_N, iv.i_ol, iv.i_max])
    pv.i_ol_lim = minnz([iv.i_ol, iv.i_max])

    pv.u_N_lim = minnz([iv.u_N, iv.u_ol, iv.u_max])
    pv.u_ol_lim = minnz([iv.u_ol, iv.u_max])


    ### calc. maximum power of cell
    if pv.u_N_lim and pv.u_N_lim < pv.u_ol_lim:                 # calc i_N based on given u_N | lim: given i_N
        pv.i_N, pv.p_N, pv.u_N = clc_i(obj, pv.T_N, obj.p, pp,
                                        u_val=pv.u_N_lim, i_lim=pv.i_N_lim)    # // A/m², W/m²
        pv.i_ol, pv.p_ol, pv.u_ol = clc_i(obj, pv.T_N, obj.p, pp,
                                        u_val=pv.u_ol_lim, i_lim=pv.i_ol_lim)    # // A/m², W/m²
        #print(f'Dev. in u_N: 1-> {pv.u_N}, 2-> {u_out}' )
        overload_possible = True
    elif pv.u_ol_lim:
        overload_possible= False
        pv.i_N, pv.p_N, pv.u_N = clc_i(obj, pv.T_N, obj.p, pp,
                                        u_val=iv.u_ol_lim, i_lim=pv.i_ol_lim)    # // A/m², W/m²

    elif pv.i_N_lim and pv.i_N_lim < pv.i_ol_lim: # otherwise use limit of i
        overload_possible = True
        pv.u_N = obj.clc_m.voltage_cell(obj, pec, pv.T_N, pv.i_N_lim, obj.p, pp=None)
        pv.u_ol = obj.clc_m.voltage_cell(obj, pec, pv.T_N, pv.i_ol_lim, obj.p, pp=None)
    else:
        overload_possible = False
        pv.u_N = obj.clc_m.voltage_cell(obj, pec, pv.T_N, pv.i_ol_lim, obj.p, pp=None)
        #print('---')
        #print('i_N: ', pv.i_N)
    prnt_attr(pv, 'i_N')
    prnt_attr(pv, 'u_N')
    ############################################################################

    pv.p_N = pv.i_N * pv.u_N        # Specific power of single cell, nominal
    pv.p_ol = pv.i_ol * pv.u_ol     # Specific power of single cell, overload

    if (not iv.A_cell) & (pv.p_N != False):
        pv.A_cell = pv.P_cell_N/pv.p_N
    else:
        pv.A_cell = iv.A_cell



    pv.P_st_N_lim = minnz([iv.P_stack_N, iv.P_stack_ol, iv.P_stack_max])
    pv.P_plnt_N_lim = minnz([pv.P_plnt_N_target, pv.P_plnt_N_target_flw,
                            iv.P_plnt_N, iv.P_plnt_ol, iv.P_plnt_max])
    pv.P_st_ol_lim = minnz([iv.P_stack_ol, iv.P_stack_max])
    pv.P_plnt_ol_lim = minnz([iv.P_plnt_ol, iv.P_plnt_max])

    print('test n_cell = ', pv.P_st_N_lim/(pv.p_N*pv.A_cell)*1e3)

    pv.n_clls_st_lim = minnz([round((pv.P_st_N_lim/(pv.p_N*pv.A_cell))*1e3),
                                iv.n_clls_st,iv.n_clls_st_max])
    pv.n_clls_plnt_lim = minnz([math.ceil(pv.i_target_flw/pv.i_N),
                                    iv.n_clls_plnt,iv.n_clls_plnt_max])
    pv.n_st_lim = minnz([iv.n_st,iv.n_st_max,
                        math.floor(pv.n_clls_plnt_lim/pv.n_clls_st_lim)])


    prnt_attr(pv, 'P_st_N_lim')
    prnt_attr(pv, 'P_plnt_N_lim')
    prnt_attr(pv, 'P_st_ol_lim')
    prnt_attr(pv, 'P_plnt_ol_lim')
    prnt_attr(pv, 'n_clls_st_lim')
    prnt_attr(pv, 'n_clls_plnt_lim')
    prnt_attr(pv, 'n_st_lim')


    pv.P_cell_N = clc_Pcell(P_spec=pv.p_N, A_cell=pv.A_cell,
                            P_stack=pv.P_st_N_lim, P_plnt=pv.P_plnt_N_lim,
                            n_cells_st=pv.n_clls_st_lim,
                            n_cells_plnt=pv.n_clls_plnt_lim,
                            P_lim = pv.p_N*1e-3*pv.A_cell)
    # ??? spec. P valid as P_lim input ????



    if overload_possible:
        pv.P_cell_ol = clc_Pcell(P_spec=pv.p_ol, A_cell=pv.A_cell,
                                P_stack=pv.P_st_ol_lim, P_plnt=pv.P_plnt_ol_lim,
                                n_cells_st=pv.n_clls_st_lim,
                                n_cells_plnt=pv.n_clls_plnt_lim,
                                P_lim = pv.P_cell_N*iv.pwr_frc_max)
        pv.fctr_ol = pv.P_cell_ol/pv.P_cell_N
    pv.i_N_act = pv.P_cell_N/(pv.A_cell * pv.u_N) *1e3 # kW/(m²*V) -> J/(m²s*V)

    pv.n_clls_plnt = minnz([math.floor(pv.P_plnt_N_lim/pv.P_cell_N),
                                        pv.n_clls_plnt_lim])
    pv.n_clls_st = max([math.ceil(pv.n_clls_plnt/pv.n_st_lim), pv.n_clls_st_lim])
    pv.n_st = math.ceil(pv.n_clls_plnt/pv.n_clls_st) #pv.n_st_lim
    pv.n_clls_st = math.floor(pv.n_clls_plnt/pv.n_st)
    pv.n_clls_plnt = pv.n_st * pv.n_clls_st

    #### adjust cell area
    pv.A_cell = (pv.P_plnt_N_lim/pv.n_clls_plnt)/(pv.p_N*1e-3)
    pv.P_cell_N = pv.p_N*1e-3 * pv.A_cell
    if overload_possible:
        pv.P_cell_ol = pv.P_cell_N * pv.fctr_ol


    print('P_plnt_act = ', pv.P_cell_N*pv.n_clls_plnt)
    prnt_attr(pv, 'n_clls_plnt_lim')
    prnt_attr(pv, 'n_clls_plnt')
    prnt_attr(pv, 'n_clls_st')
    prnt_attr(pv, 'n_st')
    pv.n_clls_plnt =  pv.n_clls_st * pv.n_st
    pv.P_st_N_act = pv.n_clls_st * pv.P_cell_N
    pv.P_plnt_N_act = pv.n_clls_st * pv.P_cell_N * pv.n_st
    if overload_possible:
          pv.P_st_ol_act = pv.n_clls_st * pv.P_cell_ol
          pv.P_plnt_ol_act = pv.P_st_ol_act * pv.n_st
    else:
        pv.P_st_ol_act = False
        pv.P_plnt_ol_act = False
    ############################################################################

    print('pv.ov_u_ol, pv.dE_rev, U_tn: ', pv.u_ol, pv.dE_rev, pv.U_tn)
    eff_u_LHV = abs(pv.dE_rev)/pv.u_ol
    eff_u_HHV = abs(pv.U_tn)/pv.u_ol
    if getattr(pv, 'iv_eff_u_HHV',None):
        if pv.iv_eff_u_HHV < eff_u_HHV:
            eff_Stack = pv.iv_eff_u_HHV
    else:
        eff_Stack = eff_u_HHV


    # ==========================================================================

    # check validity of power fraction max
    #TODO: what about min frc ???

    # pwr_frc_theo = pv.P_st_ol_act / pv.P_st_N_act
    pv.pwr_frc_max = minnz([iv.pwr_frc_max, (pv.P_st_ol_act / pv.P_st_N_act), (pv.P_cell_ol / pv.P_cell_N)])
    pv.pwr_frc_min = max([iv.pwr_frc_min, ])
    # if pwr_frc_theo != pv.iv_pwr_frc_max:
    #     print(f'Deviation in permissible power fraction w.r.t. maximum Stack power (theo: {pwr_frc_theo}) | (params: {pv.iv_pwr_frc_max})')
    #     pv.pwr_frc_max = pwr_frc_theo
    # else:
    #     pv.pwr_frc_max = pv.iv_pwr_frc_max

    #if not pv.iv_P_plnt_ol:
    #     pv.P_plnt_ol = pv.P_plnt_N * pv.pwr_frc_max
    #if not pv.iv_P_stack_ol:
    #     pv.P_stack_ol = pv.P_stack_N * pv.pwr_frc_max

    ############################################################################
    pv.P_st_occur_max = max(pv.P_st_ol_act, pv.P_st_N_act)
    pv.P_loss_max = pv.P_st_occur_max*(1-eff_Stack)



    '''
    ============================================================================
    Thermal Management
    --> calc power of pumps !
    P_pmp = P_clc /(eta_opt_pmp * eta_opt_mot) * fctr_scl
    P_clc = dV*dp
    dV_cool = kA*Q_dot / (cp dT rho)
    dV_ely =
    '''
    def clc_massflow(Pv, cp, dT):
        '''
        from clc_dimensions_hex_v01

        returns
        -------
        m_dot: flowrate
            massflow in kg/s
        '''
        m_dot = Pv/(cp*dT)
        return m_dot

    def clc_pwr_pump(V_dot=None, m_dot=None, rho=None, dp=None):
        '''
        from clc_dimensions_hex_v01
        returns
        -------
        P: float
            power of pump in kW
        '''
        if V_dot is None:
            V_dot = m_dot/ rho
        P = V_dot * dp # m3/s * Pa = J/s = W
        return P*1e-3 # // in kW

    # pv.ov_cp_coolant = pv.cp_coolant
    if not pv.cp_coolant:
        pv.cp_coolant = obj.clc_m.flws.xflws.clc_cp_H2O(obj, pv.T_N)
        #raise NotImplementedError

    # pv.ov_cp_ely = pv.iv_cp_ely
    if not pv.cp_ely:
        raise NotImplementedError

    # pv.ov_dT_min_coolant = pv.iv_dT_min_coolant
    if not pv.dT_min_coolant:
        raise NotImplementedError

    # pv.ov_dT_min_ely = pv.iv_dT_min_ely
    if not pv.dT_min_ely:
        raise NotImplementedError

    print('cp_ely: ', pv.cp_ely)
    print('cp_clnt: ', pv.cp_coolant)
    print('rho_H2O: ', pv.rho_H2O)
    print('rho_ely: ', pv.rho_ely)
    # pv.ov_P_pmp_ely = pv.iv_P_pmp_ely
    if not pv.P_pmp_ely:
        # pv.ov_dm_ely = clc_massflow(pv.P_loss_max, pv.ov_cp_ely, pv.ov_dT_min_ely)
        # if pv.iv_dm_ely:
        if iv.dV_ely_nom:
            pv.V0_ely = iv.dV_ely_nom
            # pv.dm_ely = pv.V0_ely * pv.rho_ely
            #pv.ov_dm_ely = pv.iv_dm_ely * pv.n_clls_st
        else:
            pv.V0_ely = 0
        pv.P_pmp_ely = clc_pwr_pump(V_dot=pv.V0_ely* pv.n_clls_st,
                                        dp=iv.dp_ely_cycle)

        #raise NotImplementedError

    # pv.ov_P_pmp_clnt = pv.iv_P_pmp_clnt
    if not pv.P_pmp_clnt:
        if iv.dm_clnt_max > 0:
            pv.dm_clnt = iv.dm_clnt_max
        else:
            pv.dm_clnt = clc_massflow(pv.P_loss_max*1e3, pv.cp_coolant, pv.dT_min_coolant)


        pv.P_pmp_clnt = clc_pwr_pump(m_dot=pv.dm_clnt, rho=pv.rho_H2O, dp=iv.dp_coolant_cycle)
        pv.V0_clnt = pv.dm_clnt / pv.rho_H2O
        #raise NotImplementedError
    pv.C_cw    = pv.dm_clnt * obj.bop.cp_coolant * par_dct['periphery']['corrfctr_coolant']['value']
    pv.ratio_UAHx_Ccw = iv.UA_hx0_st_ref/(iv.dm_clnt_ref * 4184)
    # pv.UA_hx   = pv.C_cw*par_dct['periphery']['corrfctr_UA_hx']['value']
    pv.UA_hx   = pv.ratio_UAHx_Ccw * pv.C_cw

    print('ratio UAHx: ', pv.ratio_UAHx_Ccw)
    print('UAHx: ', pv.UA_hx)
    print('dm_cw= ', pv.dm_clnt)
    print('C_cw= ', pv.C_cw)
    # print('CHECK line 362 !!! Hardcoded COOLANT MASSFLOW !')
    # pv.ov_dm_clnt = pv.ov_dm_clnt

    ### Scale thermal parameters
    prnt_attr(pv,'A_cell')
    prnt_attr(pv,'n_clls_st')

    fctr_d = 0.8
    # pv.fctr_scl = (pv.A_cell * ) / (iv.A_c_ref * )
    pv.fctr_len = pv.n_clls_st / iv.n_clls_st_ref
    pv.fctr_circ = ( (np.sqrt(pv.A_cell/np.pi)*1/fctr_d) /
                        (np.sqrt(iv.A_c_ref/np.pi)*1/fctr_d) )
    print('fctr_len: ', pv.fctr_len)
    print('fctr_circ: ', pv.fctr_circ)
    if iv.Ct_st == False:
        if scaling:
            pv.Ct_st = iv.Ct_st_ref * (pv.fctr_len * pv.fctr_circ) * iv.corrfctr_Ct_st
        else:
            pv.Ct_st = iv.Ct_st_ref
    print(f'Scaling = {scaling} |  scl-fctr = ',(pv.fctr_len * pv.fctr_circ))
    if iv.UA_hx0_st == False: # has no effect,due to adjustment based on C_cw
        if scaling:
            # pv.UA_hx0_st = iv.UA_hx0_st_ref * (pv.fctr_len* pv.fctr_circ)
            pv.UA_hx0_st = pv.UA_hx
        else:
            pv.UA_hx0_st = iv.UA_hx0_st_ref
    if iv.Rt_st == False:
        if scaling:
            pv.Rt_st = iv.Rt_st_ref /(pv.fctr_len * pv.fctr_circ)
        else:
            pv.Rt_st = iv.Rt_st_ref
    ### setup/tune PID
    # vals_pid_tuning = tune_pid(obj) # Returns: (Kp, KI, Kd)
    '''
    ============================================================================
    '''

    # pv.P_st_act
    # pv.P_plnt_DC
    # pv.P_plnt_AC
    '''
    ============================================================================
    '''


    # =========================================================================
    # actual values

    #par_dct['cell']['temperature_nominal']['value'] = pv.
    #par_dct['cell']['temperature_max']['value']

    par_dct['cell']['current_density_nominal']['value']     = pv.i_N_act
    par_dct['cell']['current_density_overload']['value']    = pv.i_ol # // in A/m²
    #par_dct['cell']['current_density_max']['value']         = pv.ov_
    #if not iN:
    #    iN = imx
    par_dct['cell']['voltage_nominal']['value']         = pv.u_N
    par_dct['cell']['voltage_overload']['value']        = pv.u_ol
    #par_dct['cell']['voltage_max']['value']
    par_dct['cell']['active_cell_area']['value']        = pv.A_cell

    par_dct['plant']['power_of_plant_overload']['value']    = pv.P_plnt_ol_act
    par_dct['plant']['power_of_plant_act']['value']         = pv.P_plnt_N_act
    par_dct['plant']['power_of_stack_act']['value']         = pv.P_st_N_act
    par_dct['plant']['power_of_stack_overload']['value']    = pv.P_st_ol_act
    #par_dct['plant']['power_of_stack_max']['value']
    par_dct['operation']['power_fraction_max']['value']     = pv.pwr_frc_max
    par_dct['operation']['power_fraction_min']['value']     = pv.pwr_frc_min

    par_dct['plant']['number_of_cells_in_plant_act']['value'] = pv.n_clls_plnt
    par_dct['plant']['number_of_cells_in_stack_act']['value'] = pv.n_clls_st
    #par_dct['plant']['number_of_cells_in_plant_max']['value']
    #par_dct['plant']['number_of_cells_in_stack_max']['value']
    par_dct['plant']['number_of_stacks_act']['value'] = pv.n_st

    par_dct['plant']['flowrate_H2_nominal']['value'] = 0
    par_dct['plant']['heat_capacity_st']['value'] = pv.Ct_st

    par_dct['plant']['UA_hx0_st']['value'] = pv.UA_hx0_st
    # par_dct['plant']['UA_hx0_st']['value'] = pv.UA_hx # adjusted based on C_cw

    par_dct['plant']['thermal_resistance_st']['value'] = pv.Rt_st
    par_dct['plant']['scaling_factor_len']['value'] = pv.fctr_len
    par_dct['plant']['scaling_factor_circ']['value'] = pv.fctr_circ

    par_dct['periphery']['power_rectifier_nominal']['value'] = pv.P_st_N_act
    par_dct['periphery']['power_pump_ely_nominal']['value'] = pv.P_pmp_ely
    par_dct['periphery']['power_pump_coolant_nominal']['value'] = pv.P_pmp_clnt
    par_dct['periphery']['volumetricflow_coolant_nominal']['value'] = pv.V0_clnt
    par_dct['periphery']['volumetricflow_ely_nominal']['value'] = pv.V0_ely

    par_dct['periphery']['massflow_coolant_max']['value'] = pv.dm_clnt * par_dct['periphery']['corrfctr_coolant']['value']


    # par_dct['periphery']['pid_parameters']['value'] = vals_pid_tuning # (Kp, KI, Kd)
    par_dct['electrochemistry']['u_tn']['value'] = pv.U_tn #


    # =========================================================================
    # =========== nominal values:
    print('Nominal values for simulation:')
    nml = ['u_N', 'u_ol', 'u_max', 'i_N', 'i_ol', 'i_max',
            'P_plnt_N', 'P_plnt_ol', 'P_plnt_max',
            'P_stack_N', 'P_stack_ol', 'P_stack_max',
            'n_clls_st', 'n_clls_plnt', 'n_st',
            'pwr_frc', 'prod_rate_N', 'prod_rate_max',
            'A_cell',
            'P_pmp_clnt', 'P_pmp_ely', 'dm_clnt', 'dm_ely'
            ]

    #l_sym = 80
    #pos_h0 = #hatch0
    #pos_txt0 = # text position 0
    #pos_h1 =
    for n,nm in enumerate(nml):
        attr = getattr(pv, nm, None)
        attr0 = getattr(pv, 'iv_'+nm, None)
        attr1 = getattr(pv, 'ov_'+nm, None)

        if not attr:
            attr = '-'
        if not attr1:
            attr1 = '-'
        if isinstance(attr, float):
            attr = round(attr,2)
        if isinstance(attr1, float):
            attr1 = round(attr1,2)
        if n == 0:
            print(f'|      val          |     par   |  clc  |  act  |')
        print(f'--> {nm} = {attr0} | {attr} | {attr1}')
    # =========================================================================

    ### check validity of plant power (vs. prod-rate)



    #n_cells_st = P0_Stack_max / (p_max * A_cell)
    #P_stack_N = p_N * A_cell * n_cells_st
    #P_stack_max = p_max * A_cell * n_cells_st

    ### calc i_N (P_N, args=(T_N,p,pp, i_ini))
        # obj: A_cell, n_cells_plnt


    ##ä# calc i_N (u_N, P_N, args=(T_N,p,pp, i_ini))
        # P0_stack_max
        # obj: A_cell

    ### calc i_max (u_max, args=(T_N, p, pp, i_ini, ))
        # obj: A_cell

    return par_dct #{}#pv


def minnz(lst_in):
    '''
    returns minimum value of list, which is not False or 0
    '''
    lst_f = [i for i in lst_in if ((i>0) and (i is not False))]
    if lst_f:
        min_o = min(lst_f)
    else:
        min_o = False
    return min_o

def clc_i(obj, T, p, pp, u_val=None, P_=None, n_cells=None, i_lim=None):
    '''
    i_N: float
        nominal current density // in A/m²
    P_N: float
        nominal power (either stack or plant)
    n_cells: int
        number of cells in Stack or Plant w.r.t. P_N

    Returns
    -------
    nominal current density in A/m²
    '''
    i_ini = 1000

    if not i_lim:
        ilim = [(0,1e5)] # // in A/m²
    else:
        ilim = [(0,i_lim)]

    if u_val:
        ppout = obj.clc_m.pwr.bsc_opt(obj, obj.pec, T, i_ini, ilim,p, pp,
                                        u_mx=u_val)
        i_ = ppout[0][0]
        P_ = u_val * i_

    elif P_ and n_cells:
        P_cell = P_ / n_cells # Nominal cell power
        ppout = obj.clc_m.pwr.bsc_opt(obj, obj.pec, T, i_ini, ilim, p, pp,
                                        P_in=P_cell, u_mx=None,)
        i_ = ppout[0][0]
        u_ = ppout[-1][-1] # ppout[-1] -> return from voltage_cell()

    return i_, P_, u_val


def clc_Pcell(P_spec=None, A_cell=None,
                P_stack=None, P_plnt=None,
                n_cells_st=None, n_cells_plnt=None, P_lim=False):

    d = {'0':'P0',
            '1': 'p*A_cell',
            '2': 'div: P_St/n_clls_St',
            '3': 'div: P_plnt / n_clls_plnt'}
    p0 = P_lim #plr_clc.pwr_cell(uN_cell, iN_cell) // in kW
    p1 = P_spec * A_cell *1e-3 # // -> in kW
    p2 = division(P_stack , n_cells_st)
    p3 = division(P_plnt , n_cells_plnt)

    #arr_p_N = np.array([p0, p1, p2, p3])
    #p_res = arr_p_N[arr_p_N >0]
    pi = [p0, p1, p2, p3]
    print('pi = ', pi)
    p_res = [i for i in pi if i > 0]
    if len(p_res):
        p_N = min(p_res)
        idx = pi.index(p_N)
        print(f'Derived cell power via: {d[str(idx)]} ')
        #for p_ in p_res:
        #    print(f'P_N: {pi.index(p_)} -> {p_}')
        #p_N = p_N[-1]
    else:
        print('...cell power could not be determined...')
        p_N = None
    print('p_N: ', p_N)

    return p_N

def clc_Acell(pv):
    A1 = pv.P0_Stack_max/(pv.n_clls_st * pv.P_cell_max)
    A2 = pv.P_cell_max / (pv.i_max * pv.u_max)
    A2 = pv.P_cell_N / (pv.i_N * pv.u_N)

    Ai = [A1, A2, A3]
    A_res = [i for i in Ai if i > 0]

    if len(A_res):
        A_o = min(A_res)
        #idx = pi.index(p_N)
        #print(f'Derived cell power via: {d[str(idx)]} ')
        #for p_ in p_res:
        #    print(f'P_N: {pi.index(p_)} -> {p_}')
        #p_N = p_N[-1]
    else:
        print('...cell power could not be determined...')
        A_o = None
    print('A_o: ', A_o)

    return A_o

def prnt_attr(obj, nm):
    attr = getattr(obj, nm, None)
    print(nm+': ', attr)
    return



def scale_aux_components():
    '''
    Define sizes of components
    '''
    ### Ely pump

    ### Coolant pump

    ### Heat exchanger ?

    ### water Tank

    ### Heat capacity of Cell/Stack

    ### Thermal resistance of Stack

    return

def clc_min_lop():
    '''
    determine minimum power
    '''

    return

def tune_pid(obj,):
    pid = ctrl.PID_controller()
    pid.reset()
    #pid.set_SP(setpoint)

    # pem_aux.clc_auxvals(sim, T_in[j])
    # plr_out = pem_plr.voltage_cell(sim, sim.pec, T_in[j], i_vals[j], sim.p)
    # flws_out = pem_flws.materialbalance(sim, T_in[j], i_vals[j], m_H2O_in_an,
    #                                    sim.p, c_in, n_in)
    # -> gnrl_function()
    # = T(t)

    # bump test -> 0...P_max


    obj.logger.warning('Still manual pid-tuning active')
    #pid.tune_aut()
    pid.tune_man(1,0.05,0.5)
    return pid.get_tuning_par()

###############################################################################

def clc_pwr_vals_old(obj, bsc_par, par_dct):
    '''
    old / original version of power vals calc
    '''
    '''
    ver = bsc_par['clc_ver']
    tec = bsc_par['tec_el'].lower()
    plr_clc = impm('epos.clc.' +tec+ '_plr_' + ver['plr'])
    pwr_clc = impm('epos.clc.' +tec+ '_pwr_' + ver['pwr'])
    gnrl_pwr_clc = impm('epos.clc.gnrl_pwr_' + ver['pwr'])
    flws_clc = impm('epos.clc.' +tec+ '_flws_' + ver['flws'])
    '''
    obj.clc_m = fx.ini_clc_versions(obj, bsc_par)

    #print('par_dct[cell]: ', par_dct['cell'])
    #TODO: belows code redundant (see namedtuples further below!)
    T_N         = par_dct['cell']['temperature_nominal']['value']
    T_max       = par_dct['cell']['temperature_max']['value']
    iN          = par_dct['cell']['current_density_nominal']['value']  # // in A/m²
    imx         = par_dct['cell']['current_density_max']['value']
    #if not iN:
    #    iN = imx
    uN_cell         = par_dct['cell']['voltage_nominal']['value']
    umx_cell         = par_dct['cell']['voltage_max']['value']
    A_cell = par_dct['cell']['active_cell_area']['value']

    pwr_plnt_max = par_dct['plant']['power_of_plant_max']['value']
    pwr_plnt_N  = par_dct['plant']['power_of_plant_nominal']['value']
    pwr_st_N    = par_dct['plant']['power_of_stack_nominal']['value']
    pwr_st_max  = par_dct['plant']['power_of_stack_max']['value']
    p_frc_max   = par_dct['operation']['power_fraction_max']['value']

    n_cells_tot = par_dct['plant']['number_of_cells_in_plant_act']['value']
    n_cells_max = par_dct['plant']['number_of_cells_in_plant_max']['value']
    n_cells_st  = par_dct['plant']['number_of_cells_in_stack_act']['value']
    n_cells_max  = par_dct['plant']['number_of_cells_in_stack_max']['value']
    n_st        = par_dct['plant']['number_of_stacks_act']['value']




    ### Nominal current density and cell voltage
    if not iN:
        if imx:
            i_ini = 0.5*imx
        else:
            i_ini = 100

        #p_ca = par_dct['operation']['nominal_electrode_pressure']['cathode']['value']
        #p_an = par_dct['operation']['nominal_electrode_pressure']['anode']['value']
        #p_in = [p_ca, p_an]
        dummy = None

        #dct_elchem = par_dct['electrochemistry']
        #pec = hd.dct_to_nt(par_dct['electrochemistry'], subkey='value')

        #dct_press = par_dct['operation']['nominal_electrode_pressure']
        p = hd.dct_to_nt(par_dct['operation']['nominal_electrode_pressure'],
                            subkey='value')
        obj.plant = hd.dct_to_nt(par_dct['plant'], subkey='value') # Plant parameters as namedtuple
        obj.pcll = hd.dct_to_nt(par_dct['cell'], subkey='value')  # Cell parameters as namedtuple
        obj.plnt_op = hd.dct_to_nt(par_dct['operation'], subkey='value')     # Operation parameters as namedtuple
        pec = obj.pec = hd.dct_to_nt(par_dct['electrochemistry'], subkey='value') # Electrochemistry parameters as namedtuple


        #print('pressures: ', p.anode, p.cathode)
        '''
        # TODO: redundant code as in simulation.ElSim.setup_sim
        elchem_dct = {}
        for key, sdct in par_dct['electrochemistry'].items():
            elchem_dct[key] = sdct['value']
        NT = namedtuple('NT',elchem_dct)
        pec = NT(**elchem_dct)
        #print('pec: ', pec)
        '''
        # create auxvals object
        obj.av = auxvals()
        #obj = None # dummy
        pp_in = obj.clc_m.flws.partial_pressure(obj, pec, T_N, p)
        #obj.clc_m.plr.clc_bubble_cvrg(obj, pec, T_N, i_ini, p, pp_in)
        #print('pp_in: ', pp_in)
        imax = [(0,imx)]
        pout = obj.clc_m.pwr.op_opt(obj, pec, T_N, i_ini, imax,
                                p, pp_in,
                                u_mx=uN_cell, ifun=obj.clc_m.plr.voltage_cell, ini=True)
        iN = pout[0]
        uN_cell = pout[-1][-1] # internal attr of obj function
        print('---> pout: ', pout)

    if (not n_cells_st) & n_st:
        n_cells_st = math.ceil(n_cells_tot / n_st)
    elif (not n_st) & n_cells_st:
        n_st = math.ceil(n_cells_tot / n_cells_st)
    else:
        n_cells_tot = n_cells_st * n_st

    ### possibilities to calc nominal cell power
    p0_N = 0 #plr_clc.pwr_cell(uN_cell, iN_cell)
    p1_N = iN * uN_cell * A_cell *1e-3 # // -> in kW
    p2_N = division(pwr_st_N , n_cells_st)
    p3_N = division(pwr_plnt_N , n_cells_tot)

    p0_max = 0 #pmx_cell
    p1_max = imx * umx_cell * A_cell *1e-3 # // in kW
    p2_max = division(pwr_st_max ,n_cells_st)
    p2_max = division(pwr_st_N * p_frc_max, n_cells_st)

    arr_p_N = np.array([p0_N, p1_N, p2_N,])
    arr_p_max = np.array([p0_max, p1_max, p2_max,])

    p_N = arr_p_N[arr_p_N >0]
    if len(p_N):
        p_N = p_N[-1]
    else:
        print('Nominal cell power could not be determined...')
    print('p_N: ', p_N)

    p_max = arr_p_max[arr_p_max > 0]
    if len(p_max) >1:
        p_max = p_max[-1]
    else:
        print('Maximum cell power could not be determined...')

    if not p_frc_max:
        p_frc_max = p_max / p_N
        par_dct['operation']['power_fraction_max'] = p_frc_max

    if not A_cell:
        A_cell = p_N / (iN_cell * uN_cell)

    if not pwr_plnt_N:
        pwr_plnt_N = p_N * n_cells_tot
        pwr_plnt_N_1 = p_N * n_cells_st * n_st
        if not (pwr_plnt_N_1 == pwr_plnt_N): # check, if calc. is consistent
            print('Deviation in calculation of nominal power')
        #else:
        #    par_dct['plant']['power_of_plant']['values']['nominal'] = pwr_plnt_N

    if not pwr_plnt_max:
        pwr_plnt_max = p_max * n_cells_tot
        pwr_plnt_max_1 = p_max * n_cells_st * n_st
        if not (pwr_plnt_max_1 == pwr_plnt_max): # check, if calc. is consistent
            print('Deviation in calculation of nominal power')
        else:
            par_dct['plant']['power_of_plant_max']['value'] = float(pwr_plnt_max)

    if not pwr_st_N:
        pwr_st_N = pwr_plnt_N_1 / n_st
        pwr_st_N_1 = p_N * n_cells_st
        if not (pwr_st_N_1 == pwr_st_N): # check, if calc. is consistent
            print('Deviation in calculation of nominal power')
        #else:
        #    par_dct['plant']['power_of_stack']['values']['nominal'] = pwr_st_N

    if not pwr_st_max:
        pwr_st_max = pwr_plnt_max_1 / n_st
        pwr_st_max_1 = p_max * n_cells_st
        if not (pwr_plnt_max_1 == pwr_plnt_max): # check, if calc. is consistent
            print('Deviation in calculation of nominal power')
        else:
            par_dct['plant']['power_of_stack_max']['value'] = pwr_st_max

    if (not n_st) & (not n_cells_st):

        n_cells_st = math.ceil(pwr_st_N / (p_N))
        pwr_st_N = p_N * n_cells_st
        ### check power again
        while pwr_st_N > pwr_st_max:
            #print('----p_N: ', p_N)
            n_cells_st -=1
            pwr_st_N = p_N * n_cells_st

        ### TODO: what about max values???
        n_st = math.ceil(pwr_plnt_N / pwr_st_N)
        pwr_plnt_N = pwr_st_N * n_st
        n_cells_tot = math.ceil(n_st * n_cells_st)

    print('---n_st: ', n_st)
    print('---n_cells_st: ', n_cells_st)
    par_dct['plant']['number_of_stacks_act']['value'] = n_st
    par_dct['plant']['number_of_cells_in_stack_act']['value'] = n_cells_st
    par_dct['plant']['number_of_cells_in_plant_act']['value'] = n_cells_tot
    par_dct['cell']['active_cell_area']['value'] = A_cell
    par_dct['plant']['power_of_stack_nominal']['value'] = pwr_st_N
    par_dct['plant']['power_of_plant_nominal']['value'] = pwr_plnt_N

    return par_dct


def check_polar(T=[303,333,353], i_st=0, i_end=30000, numi=100, p_in=[[101325],[101325]]):
    '''
    calc and plot polarisation characteristics
    '''
    i = np.linspace(i_st, i_end, numi)
    res = np.zeros((len(T),len(i)))
    labs=[]
    for j in range(T):
        for k in range(len(i)):
            res[j,k] = plr_fun(T[j],i[k],p)

    return

### auxilliary calculations
################################################################################

def division(n, d):
    return n / d if d else 0

'''
def clc_rated_power(par_dct,):

    ini plant capacity in ScenarioSetup
    calc
    - nominal power of stack
    - nominal power of plant


    # parameters needed !

    T_nom = par_dct['cell']['temperature']['values']['nominal']
    T_max = par_dct['cell']['temperature']['values']['max']
    i_nom = par_dct['cell']['current_density']['values']['nominal']  # // in A/m²
    i_max = par_dct['cell']['current_density']['values']['max']    # // in A/m²

    partial_pressures = None
    concentrations = None
    # clc polar for i_rated / u_max
    # clc polar for i_max / u_max
    ### ini polar:
    #-> if not 'dE0_an' in par_dct['electrochemistry']: calc E0


    different ways to specify plant pwr and connected values
    - > u_max ? -> i_max (plr)
        - > if not setting u_max (None or val > 5) then no limit is implemented



    # Stack power
    nom_pwr_st = par_dct['plant']['power_of_stack']['values']['nominal']
    if nom_pwr_st == 0:
        nom_pwr_st = None
        par_dct['plant']['power_of_stack']['values']['nominal'] = clc_

    max_pwr_st = par_dct['plant']['power_of_stack']['values']['max']
    if max_pwr_st ==0:
        max_pwr_st = None
        par_dct['plant']['power_of_stack']['values']['max'] =


    # Plant power
    nom_pwr_plnt = par_dct['plant']['power_of_plant']['values']['nominal']
    if nom_pwr_plnt ==0:
        nom_pwr_plnt = None
        par_dct['plant']['power_of_plant']['values']['nominal'] =

    max_pwr_plnt = par_dct['plant']['power_of_plant']['values']['max']
    if max_pwr_plnt ==0:
        max_pwr_plnt = None
        par_dct['plant']['power_of_plant']['values']['max'] =


    # Number of cells
    No_cells = par_dct['plant']['number_of_cells_in_stack']['act']
    if No_cells ==0:
        No_cells = None
        par_dct['plant']['number_of_cells_in_stack']['act'] =

    # Number of stacks
    No_stacks = par_dct['plant']['number_of_stacks']['act']
    if No_stacks ==0:
        No_stacks = None
        par_dct['plant']['number_of_stacks']['act'] =


    pfrc = par_dct['operation']['maximum_power_fraction']
    if pfrc == 0:
    p_max = i_max * u_max
    p_nom = i_nom * u_nom




    ### thermal resistance of stacks

    ### heat capacity of stacks

    ### PID parameters - temp cntrl

    ### PID parameters -


    return par_dct
'''
