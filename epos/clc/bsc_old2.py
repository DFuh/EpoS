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
    obj.pplnt = hd.dct_to_nt(par_dct['plant'], subkey='value') # Plant parameters as namedtuple
    obj.pcll = hd.dct_to_nt(par_dct['cell'], subkey='value')  # Cell parameters as namedtuple
    obj.pop = hd.dct_to_nt(par_dct['operation'], subkey='value')     # Operation parameters as namedtuple
    obj.pec = hd.dct_to_nt(par_dct['electrochemistry'], subkey='value') # Electrochemistry parameters as namedtuple
    obj.bop = hd.dct_to_nt(par_dct['periphery'], subkey='value') # Periphery parameters as namedtuple
    obj.p = hd.dct_to_nt(par_dct['operation']['nominal_electrode_pressure'],
                                subkey='value')
    hd.ini_auxvals(obj, par_dct)

    pv = auxvals() # Ini dataclass
    #pv.name = 'AuxVals' !!! Not working properly !

    #print('Dataclass: ', pv)
    #print('par_dct[cell]: ', par_dct['cell'])
    #TODO: belows code redundant (see namedtuples further below!)
    pv.T_N              = par_dct['cell']['temperature_nominal']['value']
    pv.T_max            = par_dct['cell']['temperature_max']['value']
    pv.iv_i_N           = par_dct['cell']['current_density_nominal']['value']  # // in A/m²
    pv.iv_i_ol          = par_dct['cell']['current_density_overload']['value']  # // in A/m²
    pv.iv_i_max         = par_dct['cell']['current_density_max']['value']
    #if not iN:
    #    iN = imx
    pv.iv_u_N           = par_dct['cell']['voltage_nominal']['value']
    pv.iv_u_ol          = par_dct['cell']['voltage_overload']['value']
    pv.iv_u_max         = par_dct['cell']['voltage_max']['value']
    pv.iv_A_cell        = par_dct['cell']['active_cell_area']['value']
    pv.iv_p_N = False
    pv.iv_p_ol = False

    pv.iv_P_plnt_max     = par_dct['plant']['power_of_plant_max']['value']
    pv.iv_P_plnt_ol     = par_dct['plant']['power_of_plant_overload']['value']
    pv.iv_P_plnt_N      = par_dct['plant']['power_of_plant_nominal']['value']
    pv.iv_P_stack_N     = par_dct['plant']['power_of_stack_nominal']['value']
    pv.iv_P_stack_ol    = par_dct['plant']['power_of_stack_overload']['value']
    pv.iv_P_stack_max  = par_dct['plant']['power_of_stack_max']['value']
    pv.iv_pwr_frc_max       = par_dct['operation']['power_fraction_max']['value']
    pv.iv_pwr_frc_min       = par_dct['operation']['power_fraction_min']['value']

    pv.iv_n_clls_plnt   = par_dct['plant']['number_of_cells_in_plant_act']['value']
    pv.iv_n_clls_st     = par_dct['plant']['number_of_cells_in_stack_act']['value']
    pv.iv_n_clls_plnt_max = par_dct['plant']['number_of_cells_in_plant_max']['value']
    pv.iv_n_clls_st_max = par_dct['plant']['number_of_cells_in_stack_max']['value']
    pv.iv_n_st          = par_dct['plant']['number_of_stacks_act']['value']
    pv.iv_n_st_max = par_dct['plant']['number_of_stacks_max']['value']

    pv.iv_Ct_st_ref = par_dct['plant']['heat_capacity_st_ref']['value']
    pv.iv_Ct_st = par_dct['plant']['heat_capacity_st']['value']
    pv.iv_n_clls_st_ref = par_dct['plant']['number_of_cells_in_stack_ref']['value']
    pv.A_c_ref = par_dct['plant']['active_cell_area_ref']['value']
    pv.iv_Rt_st_ref = par_dct['plant']['thermal_resistance_st_ref']['value']
    pv.iv_Rt_st = par_dct['plant']['thermal_resistance_st']['value']
    pv.iv_UA_hx0_st_ref = par_dct['plant']['UA_hx0_st_ref']['value']
    pv.iv_UA_hx0_st = par_dct['plant']['UA_hx0_st']['value']

    pv.prod_rate_N      = par_dct['plant']['flowrate_H2_nominal']['value']

    #pv.iv_P_rect = obj.bop['power_rectifier_nominal']['value']
    pv.iv_P_pmp_ely = par_dct['periphery']['power_pump_ely_nominal']['value']
    pv.iv_P_pmp_clnt = par_dct['periphery']['power_pump_coolant_nominal']['value']
    pv.iv_cp_coolant = par_dct['periphery']['cp_coolant']['value']
    pv.iv_cp_ely = par_dct['periphery']['cp_ely']['value']
    pv.iv_dT_min_ely = par_dct['periphery']['dT_min_ely']['value']
    pv.iv_dT_min_coolant = par_dct['periphery']['dT_min_coolant']['value']
    pv.iv_dp_ely_cycle = par_dct['periphery']['dp_ely_cycle']['value']
    pv.iv_dp_coolant_cycle = par_dct['periphery']['dp_coolant_cycle']['value']

    pv.iv_dm_clnt_max = par_dct['periphery']['massflow_coolant_max']['value']

    pv.iv_dV_ely_nom = par_dct['periphery']['volumetricflow_ely_nominal']['value']
    #pv.iv_dm_ely_max = par_dct['periphery']['massflow_ely_max']['value']

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

    ### calc. maximum power of cell

    if pv.iv_u_N:                 # calc i_N based on given u_N | lim: given i_N
        i_lim = minnz([pv.iv_i_N, pv.iv_i_ol, pv.iv_i_max])
        pv.i_N, pv.p_N, pv.u_N = clc_i(obj, pv.T_N, obj.p, pp,
                                        u_val=pv.iv_u_N, i_lim=i_lim)    # // A/m², W/m²
        #print(f'Dev. in u_N: 1-> {pv.u_N}, 2-> {u_out}' )
    elif pv.iv_i_N:
        pv.i_N = pv.iv_i_N
        pv.u_N = obj.clc_m.voltage_cell(obj, pec, pv.T_N, pv.i_N, obj.p, pp=None)
    else:
        pv.i_N = pv.iv_i_max
        pv.u_N = obj.clc_m.voltage_cell(obj, pec, pv.T_N, pv.i_N, obj.p, pp=None)
        #print('---')
        #print('i_N: ', pv.i_N)

    if pv.iv_u_ol or pv.iv_u_max: # calc i_ol based on given u_N | lim: given i_N
        i_lim = minnz([pv.iv_i_ol, pv.iv_i_max])
        u_lim = minnz([pv.iv_u_ol, pv.iv_u_max])
        pv.i_ol, pv.p_ol, pv.u_ol = clc_i(obj, pv.T_N, obj.p, pp, u_val=u_lim,
                                    P_=None, n_cells=None, i_lim=i_lim)  # A/m² , W/m²
        #print(f'Dev. in u_max: 1-> {pv.u_max}, 2-> {u_out}' )
        #pv.i_lim = pv.i_ol
    else:
        pv.i_ol, pv.u_ol = pv.i_N, pv.u_N

    pv.p_N = pv.i_N * pv.u_N        # Specific power of single cell, nominal
    pv.p_ol = pv.i_ol * pv.u_ol     # Specific power of single cell, overload
    #if not pv.i_N:
    #    pv.i_N = min([i for i in [pv.i_ol, pv.i_max] if i])
    #    pv.i_max, pv.p_max = pv.i_N, pv.p_N
        #u_max = 0

    pv.P_cell_N = clc_Pcell(P_spec=pv.p_N, A_cell=pv.iv_A_cell,
                            P_stack=pv.iv_P_stack_N, P_plnt=pv.iv_P_plnt_N,
                            n_cells_st=pv.iv_n_clls_st, n_cells_plnt=pv.iv_n_clls_plnt, P_lim = pv.iv_p_N)
    if (not pv.iv_A_cell) & (pv.p_N != False):
        pv.A_cell = pv.P_cell_N/pv.p_N
    else:
        pv.A_cell = pv.iv_A_cell

    pv.P_cell_ol = clc_Pcell(P_spec=pv.p_ol, A_cell=pv.iv_A_cell,
                            P_stack=pv.iv_P_stack_ol, P_plnt=pv.iv_P_plnt_ol,
                            n_cells_st=pv.iv_n_clls_st, n_cells_plnt=pv.iv_n_clls_plnt, P_lim = pv.iv_p_ol)

    ############################################################################


    pv.n_clls_st_lim = minnz([pv.iv_n_clls_st, pv.iv_n_clls_st_max, np.inf])
    pv.n_st_lim = minnz([pv.iv_n_st, pv.iv_n_st_max, np.inf])
    pv.n_clls_plnt_lim = minnz([pv.iv_n_clls_plnt, pv.iv_n_clls_plnt_max,
                                pv.n_st_lim * pv.n_clls_st_lim])

    # pv.P_plnt_lim = minnz([pv.iv_P_plnt_N])
    ##### lim values
    pv.P_plnt_N_lim = minnz([pv.iv_P_plnt_N, pv.iv_P_plnt_ol, pv.iv_P_plnt_max, pv.n_clls_plnt_lim*pv.P_cell_N])
    pv.P_plnt_ol_lim = minnz([pv.iv_P_plnt_ol, pv.iv_P_plnt_max, pv.n_clls_plnt_lim*pv.P_cell_ol])

    # pv.P_target_plnt_N = minnz([pv.P_target_plnt_N, pv.P_target_plnt_ol])

    pv.P_st_N_lim = minnz([pv.iv_P_stack_N, pv.iv_P_stack_ol, pv.iv_P_stack_max, pv.n_clls_st_lim*pv.P_cell_N])
    pv.P_st_ol_lim = minnz([pv.iv_P_stack_ol, pv.iv_P_stack_max, pv.n_clls_st_lim*pv.P_cell_ol])

    pv.n_clls_plnt = math.floor(pv.P_plnt_N_lim / pv.P_cell_N)
    print('pv.n_clls_plnt: ',pv.n_clls_plnt)


    print('n_clls_st_lim: ', pv.n_clls_st_lim)
    pv.n_st_clc = minnz([math.ceil(pv.n_clls_plnt/ pv.n_clls_st_lim),math.floor(pv.P_plnt_N_lim/pv.P_st_N_lim)])
    if pv.n_st_clc <= 0:
        pv.n_st_clc = 1
    print('pv.n_st_lim: ',pv.n_st_lim)

    if pv.n_st_clc > pv.n_st_lim:
        pv.n_st = pv.n_st_lim
    else:
        pv.n_st = pv.n_st_clc

    pv.n_clls_st_clc =math.floor(pv.n_clls_plnt/pv.n_st)
    print('pv.n_clls_st_clc: ', pv.n_clls_st_clc)

    if pv.n_clls_st_clc > pv.n_clls_st_lim:
        pv.n_clls_st = pv.n_clls_st_lim
    else:
        pv.n_clls_st = pv.n_clls_st_clc

    pv.P_st_N = pv.P_cell_N * pv.n_clls_st
    print('pv.P_st_N: ',pv.P_st_N)
    print('pv.n_st: ',pv.n_st)

    pv.P_plnt_N = pv.P_st_N * pv.n_st

    pv.n_clls_plnt = pv.n_st*pv.n_clls_st

    ### final clc N
    pv.P_st_N_fin = pv.n_clls_st * pv.P_cell_N
    pv.P_plnt_N_fin = pv.P_st_N_fin * pv.n_st


    ### clc ol


    pv.P_cell_ol = clc_Pcell(P_spec=pv.p_ol, A_cell=pv.A_cell,
                            P_stack=pv.P_st_ol_lim, P_plnt=pv.P_plnt_ol_lim,
                            n_cells_st=pv.n_clls_st, n_cells_plnt=pv.n_clls_plnt, P_lim = pv.P_cell_ol)
    pv.P_st_ol = pv.n_clls_st * pv.P_cell_ol
    if pv.P_st_ol > pv.P_st_ol_lim:
        pv.P_st_ol = pv.P_st_ol_lim
        print('sth. might went wrong in ol calc...')

    pv.P_plnt_ol = pv.n_st * pv.P_st_ol
    if pv.P_plnt_ol > pv.P_plnt_ol_lim:
        pv.P_plnt_ol = pv.P_st_ol_lim
        print('sth. might went wrong in ol calc (plnt)...')

    pv.pwr_frc_min = minnz([pv.iv_pwr_frc_min,])
    pv.pwr_frc_max = minnz([pv.iv_pwr_frc_max, pv.P_cell_ol/pv.P_cell_N ]) #pv.P_cell_ol/pv.P_cell_N])
    print('pwr_frc_max (iv): ', pv.pwr_frc_max)
    print('pwr_frc_max (clc): ', pv.P_cell_ol/pv.P_cell_N )
    print('pwr_frc_max (out): ', pv.pwr_frc_max)

    pv.P_cell_ol = pv.P_cell_N * pv.pwr_frc_max
    pv.P_st_ol_fin = pv.n_clls_st * pv.P_cell_ol
    pv.P_plnt_ol_fin = pv.P_st_ol_fin * pv.n_st

    # pv.P_target_stack_N = minnz([pv.P_target_stack_N, pv.P_target_stack_ol])

    # pv.P_target_stack_N = minnz([pv.P_target_stack_N, math.ceil(pv.P_target_plnt_N/pv.n_st_lim)])
    # pv.P_target_stack_ol = minnz([pv.P_target_stack_ol, math.ceil(pv.P_target_plnt_ol/pv.n_st_lim)])


    ######



    if False:
        pv.n_clls_st_lim = minnz([pv.n_clls_st_lim, math.ceil(pv.P_target_stack_N/pv.P_cell_N)])
        pv.n_st_lim = minnz([pv.n_st_lim, math.ceil(pv.P_target_plnt_N/pv.P_cell_N)])
        pv.n_clls_plnt_lim = minnz([pv.iv_n_clls_plnt, pv.iv_n_clls_plnt_max,
                                    pv.n_st_lim * pv.n_clls_plnt_lim])

        pv.P_stack_N = pv.n_clls_st_lim * pv.P_cell_N
        pv.P_plnt_N = pv.P_stack_N * pv.n_st_lim
        pv.P_stack_ol = pv.P_target_stack_N * pv.pwr_frc_max
        pv.P_plnt_ol = pv.P_stack_ol * pv.n_st_lim



    ### clc i_act
    pv.ov_i_N, pv.ov_p_N, pv.ov_u_N = clc_i(obj, pv.T_N, obj.p, pp,
                                            u_val=None, P_=pv.P_cell_N/pv.A_cell,
                                            n_cells=1, i_lim=pv.i_N)
    pv.ov_i_ol, pv.ov_p_ol, pv.ov_u_ol = clc_i(obj, pv.T_N, obj.p, pp,
                                            u_val=None, P_=pv.P_cell_ol/pv.A_cell,
                                            n_cells=1, i_lim=pv.i_ol)

    #print('======================================')
    #print('P_Stack: ', pv.P_stack_N)
    #print('P_cell: ', pv.P_cell_N)
    #print('i_N: ', pv.i_N)
    #print('======================================')

    ### Stack power and number of cells
    #if pv.iv_n_clls_st:
    #    pv.n_clls_st_lim = pv.iv_n_clls_st
    #elif pv.iv_n_clls_st_max:
    #    pv.n_clls_st_lim = pv.iv_n_clls_st_max

        #if 0 < (pv.iv_n_clls_st_max*1) < pv.n_clls_st:

    # pv.P_stack_N = minnz([pv.iv_P_stack_N, pv.iv_P_stack_ol, pv.iv_P0_stack_max])

    #print('pv.P_stack_N: ', pv.P_stack_N)
    if False:
        if not pv.P_stack_N:
            pv.P_stack_N = pv.n_clls_st * pv.P_cell_N
        else:
            pv.n_clls_st = minnz([math.ceil(pv.P_stack_N / pv.P_cell_N), pv.n_clls_st_lim])

    if False:
        ### overload
        pv.P_stack_ol = minnz([pv.iv_P_stack_ol, pv.iv_P0_stack_max])
        if not pv.P_stack_ol:
            pv.P_stack_ol = pv.n_clls_st * pv.P_cell_ol
        else:
            pv.n_clls_st = minnz([math.ceil(pv.P_stack_ol / pv.P_cell_ol), pv.n_clls_st_lim])

        if (pv.P_stack_ol / pv.P_stack_N) >(pv.iv_pwr_frc*1):
        #pv.pwr_frc_lim = min(pv.P_stack_ol / pv.P_stack_N, pv.pwr_frc_max)
            pv.P_stack_ol = pv.P_stack_N * pv.iv_pwr_frc
            pv.n_clls_st = math.floor(pv.P_stack_ol / pv.P_cell_ol)
        #prnt_attr(pv, 'n_clls_st')
    ### Plant power
    if False:
        pv.n_st_lim = minnz([pv.iv_n_st, pv.iv_n_st_max])
        pv.P_plnt_N = minnz([pv.iv_P_plnt_N, pv.iv_P_plnt_ol])

    if False:
        if not pv.P_plnt_N:
            pv.P_plnt_N = pv.n_st * pv.P_stack_N
        else:
            pv.n_st = minnz([math.ceil(pv.P_plnt_N / pv.P_stack_N), pv.n_st_lim])

    # pv.n_clls_plnt_lim = minnz([pv.iv_n_clls_plnt, pv.iv_n_clls_plnt_max])
    if False:
        pv.n_clls_plnt = pv.n_st_lim * pv.n_clls_st_lim
        pv.n_clls_plnt = minnz([pv.n_clls_plnt_lim, pv.n_clls_plnt])
    prnt_attr(pv, 'n_st_lim')

    ### actual values
    '''
    replace lines below by more efficient code
    '''
    pv.ov_P_cell_N = pv.P_cell_N
    pv.ov_P_cell_ol = pv.P_cell_ol
    pv.ov_P_stack_N = pv.ov_P_cell_N * pv.n_clls_st
    pv.ov_P_stack_ol = pv.ov_P_cell_ol * pv.n_clls_st
    pv.ov_P_plnt_ol = pv.ov_P_stack_ol * pv.n_st
    pv.ov_P_plnt_N = pv.ov_P_stack_N * pv.n_st
    pv.ov_u_N = pv.u_N
    pv.ov_u_ol = pv.u_ol
    pv.ov_i_N = pv.i_N
    pv.ov_i_ol = pv.i_ol

    pv.n_clls_plnt = pv.n_st * pv.n_clls_st
    # pv.n_clls_st = pv.n_clls_st
    prnt_attr(pv, 'n_clls_plnt')
    if False:
        P_plnt_N0 = pv.n_clls_plnt * pv.P_cell_N
        P_plnt_N1 = pv.n_st_lim * pv.P_stack_N
        if P_plnt_N0 != P_plnt_N1:
            print(f'Deviation in Plant power: N0= {P_plnt_N0} | N1= {P_plnt_N1}')
            pv.ov_P_plnt_N = minnz([P_plnt_N0, P_plnt_N1])#P_plnt_N1
        else:
            pv.ov_P_plnt_N = P_plnt_N1
        pv.ov_P_plnt_ol = pv.ov_P_stack_ol * pv.n_st_lim
        pv.pwr_frc = pv.ov_P_stack_ol / pv.P_stack_N

    #pv.ov_P_plnt_N  = pv.P_stack_act * pv.n_st
    #pv.ov_P_plnt_ol = pv.P_stack_ol * pv.n_st

    #prnt_attr(pv, 'n_clls_st')
    #print('test: ', type(getattr(pv, 'n_clls_st')))

    ### voltage efficiency (worst case)

    print('pv.ov_u_ol, pv.dE_rev, U_tn: ', pv.ov_u_ol, pv.dE_rev, pv.U_tn)
    eff_u_LHV = abs(pv.dE_rev)/pv.ov_u_ol
    eff_u_HHV = abs(pv.U_tn)/pv.ov_u_ol
    if getattr(pv, 'iv_eff_u_HHV',None):
        if pv.iv_eff_u_HHV < eff_u_HHV:
            eff_Stack = pv.iv_eff_u_HHV
    else:
        eff_Stack = eff_u_HHV
    pv.P_loss_max = pv.ov_P_stack_ol*(1-eff_Stack)

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

    pv.ov_cp_coolant = pv.iv_cp_coolant
    if not pv.ov_cp_coolant:
        pv.ov_cp_coolant = obj.clc_m.flws.xflws.clc_cp_H2O(obj, pv.T_N)
        #raise NotImplementedError

    pv.ov_cp_ely = pv.iv_cp_ely
    if not pv.ov_cp_ely:
        raise NotImplementedError

    pv.ov_dT_min_coolant = pv.iv_dT_min_coolant
    if not pv.ov_dT_min_coolant:
        raise NotImplementedError

    pv.ov_dT_min_ely = pv.iv_dT_min_ely
    if not pv.ov_dT_min_ely:
        raise NotImplementedError

    print('cp_ely: ', pv.ov_cp_ely)
    print('cp_clnt: ', pv.ov_cp_coolant)
    print('rho_H2O: ', pv.rho_H2O)
    print('rho_ely: ', pv.rho_ely)
    pv.ov_P_pmp_ely = pv.iv_P_pmp_ely
    if not pv.ov_P_pmp_ely:
        # pv.ov_dm_ely = clc_massflow(pv.P_loss_max, pv.ov_cp_ely, pv.ov_dT_min_ely)
        # if pv.iv_dm_ely:
        if pv.iv_dV_ely_nom:
            pv.ov_V0_ely = pv.iv_dV_ely_nom
            pv.ov_dm_ely = pv.ov_V0_ely * pv.rho_ely
            #pv.ov_dm_ely = pv.iv_dm_ely * pv.n_clls_st
        else:
            pv.ov_V0_ely = 0
        pv.ov_P_pmp_ely = clc_pwr_pump(V_dot=pv.ov_V0_ely* pv.n_clls_st,
                                        dp=pv.iv_dp_ely_cycle)

        #raise NotImplementedError

    pv.ov_P_pmp_clnt = pv.iv_P_pmp_clnt
    if not pv.ov_P_pmp_clnt:
        if pv.iv_dm_clnt_max > 0:
            pv.ov_dm_clnt = pv.iv_dm_clnt_max
        else:
            pv.ov_dm_clnt = clc_massflow(pv.P_loss_max*1e3, pv.ov_cp_coolant, pv.ov_dT_min_coolant)


        pv.ov_P_pmp_clnt = clc_pwr_pump(m_dot=pv.ov_dm_clnt, rho=pv.rho_H2O, dp=pv.iv_dp_coolant_cycle)
        pv.ov_V0_clnt = pv.ov_dm_clnt / pv.rho_H2O
        #raise NotImplementedError
    pv.C_cw    = pv.ov_dm_clnt * obj.bop.cp_coolant * par_dct['periphery']['corrfctr_coolant']['value']
    pv.UA_hx   = pv.C_cw*par_dct['periphery']['corrfctr_UA_hx']['value']
    # print('CHECK line 362 !!! Hardcoded COOLANT MASSFLOW !')
    # pv.ov_dm_clnt = pv.ov_dm_clnt

    ### Scale thermal parameters
    pv.fctr_scl = (pv.A_cell * pv.n_clls_st) / (pv.A_c_ref * pv.iv_n_clls_st_ref)
    if pv.iv_Ct_st == False:
        pv.Ct_st = pv.iv_Ct_st_ref * pv.fctr_scl
    if pv.iv_UA_hx0_st == False:
        pv.UA_hx0_st = pv.iv_UA_hx0_st_ref * pv.fctr_scl
    if pv.iv_Rt_st == False:
            pv.Rt_st = pv.iv_Rt_st_ref /pv.fctr_scl
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

    # ==========================================================================

    # check validity of power fraction max
    #TODO: what about min frc ???
    '''
    pwr_frc_theo = pv.P_stack_ol / pv.P_stack_act
    if pwr_frc_theo != pv.iv_pwr_frc_max:
        print(f'Deviation in permissible power fraction w.r.t. maximum Stack power (theo: {pwr_frc_theo}) | (params: {pv.iv_pwr_frc_max})')
        pv.pwr_frc_max = pwr_frc_theo
    else:
        pv.pwr_frc_max = pv.iv_pwr_frc_max
    '''
    #if not pv.iv_P_plnt_ol:
    #     pv.P_plnt_ol = pv.P_plnt_N * pv.pwr_frc_max
    #if not pv.iv_P_stack_ol:
    #     pv.P_stack_ol = pv.P_stack_N * pv.pwr_frc_max
    # =========================================================================
    # actual values

    #par_dct['cell']['temperature_nominal']['value'] = pv.
    #par_dct['cell']['temperature_max']['value']

    par_dct['cell']['current_density_nominal']['value']     = pv.ov_i_N
    par_dct['cell']['current_density_overload']['value']    = pv.ov_i_ol # // in A/m²
    #par_dct['cell']['current_density_max']['value']         = pv.ov_
    #if not iN:
    #    iN = imx
    par_dct['cell']['voltage_nominal']['value']         = pv.ov_u_N
    par_dct['cell']['voltage_overload']['value']        = pv.ov_u_ol
    #par_dct['cell']['voltage_max']['value']
    par_dct['cell']['active_cell_area']['value']        = pv.A_cell

    par_dct['plant']['power_of_plant_overload']['value']    = pv.ov_P_plnt_ol
    par_dct['plant']['power_of_plant_act']['value']         = pv.ov_P_plnt_N
    par_dct['plant']['power_of_stack_act']['value']         = pv.ov_P_stack_N
    par_dct['plant']['power_of_stack_overload']['value']    = pv.ov_P_stack_ol
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
    par_dct['plant']['thermal_resistance_st']['value'] = pv.Rt_st
    par_dct['plant']['scaling_factor']['value'] = pv.fctr_scl

    par_dct['periphery']['power_rectifier_nominal']['value'] = pv.ov_P_stack_ol
    par_dct['periphery']['power_pump_ely_nominal']['value'] = pv.ov_P_pmp_ely
    par_dct['periphery']['power_pump_coolant_nominal']['value'] = pv.ov_P_pmp_clnt
    par_dct['periphery']['volumetricflow_coolant_nominal']['value'] = pv.ov_V0_clnt
    par_dct['periphery']['volumetricflow_ely_nominal']['value'] = pv.ov_V0_ely

    par_dct['periphery']['massflow_coolant_max']['value'] = pv.ov_dm_clnt * par_dct['periphery']['corrfctr_coolant']['value']
    par_dct['plant']['UA_hx0_st']['value'] = pv.UA_hx

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
    lst_f = [i for i in lst_in if i]
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
    p0 = P_lim #plr_clc.pwr_cell(uN_cell, iN_cell)
    p1 = P_spec * A_cell *1e-3 # // -> in kW
    p2 = division(P_stack , n_cells_st)
    p3 = division(P_plnt , n_cells_plnt)

    #arr_p_N = np.array([p0, p1, p2, p3])
    #p_res = arr_p_N[arr_p_N >0]
    pi = [p0, p1, p2, p3]
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
    attr = getattr(obj, nm)
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
