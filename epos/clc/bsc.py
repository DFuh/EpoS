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

import epos.aux.handlingdata as hd
import epos.aux.faux as fx
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

    obj.clc_m = fx.ini_clc_versions(obj, bsc_par)
    obj.pplnt = hd.dct_to_nt(par_dct['plant'], subkey='value') # Plant parameters as namedtuple
    obj.pcll = hd.dct_to_nt(par_dct['cell'], subkey='value')  # Cell parameters as namedtuple
    obj.pop = hd.dct_to_nt(par_dct['operation'], subkey='value')     # Operation parameters as namedtuple
    obj.pec = hd.dct_to_nt(par_dct['electrochemistry'], subkey='value') # Electrochemistry parameters as namedtuple
    obj.p = hd.dct_to_nt(par_dct['operation']['nominal_electrode_pressure'],
                                subkey='value')
    hd.ini_auxvals(obj, par_dct)

    pv = auxvals() # Ini dataclass
    #pv.name = 'AuxVals' !!! Not working properly !

    #print('Dataclass: ', pv)
    #print('par_dct[cell]: ', par_dct['cell'])
    #TODO: belows code redundant (see namedtuples further below!)
    pv.T_N         = par_dct['cell']['temperature_nominal']['value']
    pv.T_max       = par_dct['cell']['temperature_max']['value']
    pv.i_N          = par_dct['cell']['current_density_nominal']['value']  # // in A/m²
    pv.i_ol          = par_dct['cell']['current_density_overload']['value']  # // in A/m²
    pv.i_max         = par_dct['cell']['current_density_max']['value']
    #if not iN:
    #    iN = imx
    pv.u_N         = par_dct['cell']['voltage_nominal']['value']
    pv.u_ol         = par_dct['cell']['voltage_overload']['value']
    pv.u_max         = par_dct['cell']['voltage_max']['value']
    pv.A_cell = par_dct['cell']['active_cell_area']['value']

    pv.P_plnt_ol = par_dct['plant']['power_of_plant_max']['value']
    pv.P_plnt_N  = par_dct['plant']['power_of_plant_nominal']['value']
    pv.P_stack_N    = par_dct['plant']['power_of_stack_nominal']['value']
    pv.P_stack_ol  = par_dct['plant']['power_of_stack_overload']['value']
    pv.P0_stack_max  = par_dct['plant']['power_of_stack_max']['value']
    pv.pwr_frc_max   = par_dct['operation']['power_fraction_max']['value']

    pv.n_clls_plnt = par_dct['plant']['number_of_cells_in_plant_act']['value']
    pv.n_clls_st  = par_dct['plant']['number_of_cells_in_stack_act']['value']
    pv.n_clls_plnt_max = par_dct['plant']['number_of_cells_in_plant_max']['value']
    pv.n_clls_st_max  = par_dct['plant']['number_of_cells_in_stack_max']['value']
    pv.n_st        = par_dct['plant']['number_of_stacks_act']['value']

    pv.prod_rate_N = par_dct['plant']['flowrate_H2_nominal']['value']
    ### get values from dict
    #cell_dct = par_dct['cell']

    ### calc i_N (u_N, args=(T_N, p, pp, i_ini, ))
        # obj: A_cell
    p = obj.p
    pp = pp = obj.clc_m.flws.partial_pressure(obj, obj.pec, pv.T_N, p)


    if not pv.P_plnt_ol:
         pv.P_plnt_ol = pv.P_plnt_N * pv.pwr_frc_max
    if not pv.P_stack_ol:
         pv.P_stack_ol = pv.P_stack_N * pv.pwr_frc_max
    ### calc. power of cell
    #print('u_N: ', pv.u_N)
    if pv.u_N: # calc i_N based on given u_N | lim: given i_N
        i_lim = min([i for i in [pv.i_N, pv.i_ol, pv.i_max] if i])
        pv.i_N, pv.p_N, u_out = clc_i(obj, pv.T_N, p, pp,
                                        u_val=pv.u_N, i_lim=i_lim)    # // A/m², W/m²
        print(f'Dev. in u_N: 1-> {pv.u_N}, 2-> {u_out}' )
    else:
        pv.i_N = pv.i_max
        pv.u_N = obj.clc_m.voltage_cell(obj, pec, pv.T_N, pv.i_N, p, pp=None)
        #print('---')
        #print('i_N: ', pv.i_N)

    if pv.u_ol: # calc i_ol based on given u_N | lim: given i_N
        i_lim = min([i for i in [pv.i_ol, pv.i_max] if i])
        pv.i_ol, pv.p_ol, u_out = clc_i(obj, pv.T_N, p, pp, u_val=pv.u_max,
                                    P_=None, n_cells=None, i_lim=i_lim)  # A/m² , W/m²
        print(f'Dev. in u_max: 1-> {pv.u_max}, 2-> {u_out}' )
    else:
        pv.i_ol, pv.u_ol = pv.i_N, pv.u_N

    pv.p_N = pv.i_N * pv.u_N
    pv.p_ol = pv.i_ol * pv.u_ol
    #if not pv.i_N:
    #    pv.i_N = min([i for i in [pv.i_ol, pv.i_max] if i])
    #    pv.i_max, pv.p_max = pv.i_N, pv.p_N
        #u_max = 0

    pv.P_cell_N = clc_Pcell(P_spec=pv.p_N, A_cell=pv.A_cell,
                            P_stack=pv.P_stack_N, P_plnt=pv.P_plnt_N,
                            n_cells_st=pv.n_clls_st, n_cells_plnt=pv.n_clls_plnt)
    if (not pv.A_cell) & (pv.p_N != False):
        pv.A_cell = pv.P_cell_N/pv.p_N

    pv.P_cell_ol = clc_Pcell(P_spec=pv.p_ol, A_cell=pv.A_cell,
                            P_stack=pv.P_stack_ol, P_plnt=pv.P_plnt_ol,
                            n_cells_st=pv.n_clls_st, n_cells_plnt=pv.n_clls_plnt)

    print('======================================')
    print('P_Stack: ', pv.P_stack_N)
    print('P_cell: ', pv.P_cell_N)
    print('i_N: ', pv.i_N)
    print('======================================')
    if not pv.n_clls_plnt:
        pv.n_clls_plnt = math.ceil(pv.P_plnt_N / pv.P_cell_N)
        if pv.n_clls_plnt_max & (pv.n_clls_plnt_max < pv.n_clls_plnt):
            pv.n_clls_plnt = pv.n_clls_plnt_max

    # TODO: use max or nom ???
    if not pv.n_clls_st:
        #n_clls_st0 = pv.P0_stack_max / pv.P_cell_ol
        #n_clls_st1 = pv.P_stack_N / pv.P_cell_N
        #print(f'Dev. in n_clls_st: 1-> {n_clls_st0}, 2-> {n_clls_st1}' )
        pv.n_clls_st = math.ceil(pv.P_stack_N / pv.P_cell_N)
        if pv.n_clls_st_max & (pv.n_clls_st_max < pv.n_clls_st):
            pv.n_clls_st = pv.n_clls_st_max



    pv.P_stack_act = pv.n_clls_st * pv.P_cell_N
    if not pv.P_stack_N:
        pv.P_stack_N = pv.P_cell_N * pv.n_clls_st

    if not pv.P_stack_ol:
        pv.P_stack_max = pv.P_cell_max * pv.n_clls_st

    if not pv.n_st:
        n_st0 = math.ceil(division(pv.P_plnt_N, pv.P_stack_N))
        n_st1 = math.ceil(pv.P_plnt_N/pv.P0_stack_max)
        print(f'Dev. in n_st: 1-> {n_st0}, 2-> {n_st1}' )
        pv.n_st = n_st0

    print('cells in plant old, new: ', pv.n_clls_plnt, (pv.n_st * pv.n_clls_st))
    pv.P_plnt_act = pv.n_st * pv.P_stack_act

    #if not pv.A_cell:
    #    pv.A_cell = clc_Acell(pv)
        #A_cell = P0_Stack_max/(n_clls_st * P_cell_max)
    '''
    if not pv.u_N:
        pv.i_N, pv.p_N, pv.u_N = clc_i(obj, pv.T_N, p, pp,
                            u_val=None, P_=pv.P_cell_N, n_cells=1, i_lim=None)
    if not pv.u_max:
        pv.i_max, pv.p_max, pv.u_max = clc_i(obj, pv.T_N, p, pp, u_val=None, P_=pv.P_cell_ol, n_cells=1, i_lim=None)

    if not pv.i_N:
        pv.i_N, pv.p_N, pv.u_N = clc_i(obj, pv.T_N, p, pp, u_val=pv.u_N, P_=None, n_cells=None, i_lim=None)
    if not pv.i_max:
        pv.i_max, pv.p_max, pv.u_max = clc_i(obj, pv.T_N, p, pp, u_val=pv.u_max, P_=None, n_cells=None, i_lim=None)



    if not pv.P_stack_N:
        pv.P_stack_N = pv.P_cell_N * pv.n_clls_st

    #if not pv.n_clls_plnt:
    #    pv.n_clls_plnt = pv.n_st * pv.n_clls_st
    '''


    # ==========================================================================

    # check validity of power fraction max
    #TODO: what about min frc ???
    pwr_frc_theo = pv.P_stack_N / pv.P_stack_N
    if pwr_frc_theo != pv.pwr_frc_max:
        print(f'Deviation in permissible power fraction w.r.t. maximum Stack power (theo: {pwr_frc_theo}) | (params: {pv.pwr_frc_max})')
        pv.pwr_frc_max = pwr_frc_theo

    # =========================================================================
    # =========== nominal values:
    print('Nominal values for simulation:')
    nml = ['u_N', 'u_ol', 'u_max', 'i_N', 'i_ol', 'i_max',
            'P_plnt_N', 'P_plnt_act', 'P_plnt_max',
            'P_stack_N', 'P_stack_act', 'P_stack_ol', 'P_stack_max',
            'n_clls_st', 'n_clls_plnt', 'n_st',
            'pwr_frc_max', 'prod_rate_N', 'prod_rate_max'
            ]

    for nm in nml:
        attr = getattr(pv, nm, None)
        if not attr:
            attr = 'Not defined'
        print(f'--> {nm} = {attr}')
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

    return {}#pv


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
                n_cells_st=None, n_cells_plnt=None):

    d = {'0':0,
            '1': 'p*A_cell',
            '2': 'div: P_St/n_clls_St',
            '3': 'div: P_plnt / n_clls_plnt'}
    p0 = 0 #plr_clc.pwr_cell(uN_cell, iN_cell)
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
