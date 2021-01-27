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

import epos.aux.handlingdata as hd

#TODO: take into account power of peripherie!!!

def clc_pwr_vals(bsc_par, par_dct):

    ver = bsc_par['clc_ver']
    tec = bsc_par['tec_el'].lower()
    plr_clc = impm('epos.clc.' +tec+ '_plr_' + ver['plr'])
    pwr_clc = impm('epos.clc.' +tec+ '_pwr_' + ver['pwr'])
    gnrl_pwr_clc = impm('epos.clc.gnrl_pwr_' + ver['pwr'])
    flws_clc = impm('epos.clc.' +tec+ '_flws_' + ver['flws'])

    #print('par_dct[cell]: ', par_dct['cell'])
    T_N         = par_dct['cell']['temperature']['values']['nominal']
    T_max       = par_dct['cell']['temperature']['values']['max']
    iN          = par_dct['cell']['current_density']['values']['nominal']  # // in A/m²
    imx         = par_dct['cell']['current_density']['values']['max']
    uN_cell         = par_dct['cell']['cell_voltage']['values']['nominal']
    umx_cell         = par_dct['cell']['cell_voltage']['values']['max']
    A_cell = par_dct['cell']['active_cell_area']['value']

    pwr_plnt_max = par_dct['plant']['power_of_plant']['values']['max']
    pwr_plnt_N  = par_dct['plant']['power_of_plant']['values']['nominal']
    pwr_st_N    = par_dct['plant']['power_of_stack']['values']['nominal']
    pwr_st_max  = par_dct['plant']['power_of_stack']['values']['max']
    p_frc_max   = par_dct['operation']['maximum_power_fraction']

    n_cells_tot = par_dct['plant']['number_of_cells_in_plant']['act']
    n_cells_st  = par_dct['plant']['number_of_cells_in_stack']['act']
    n_st        = par_dct['plant']['number_of_stacks']['act']




    ### Nominal current density and cell voltage
    if not iN:
        if imx:
            i_ini = 0.5*imx
        else:
            i_ini = 1

        p_ca = par_dct['operation']['nominal_electrode_pressure']['cathode']['value']
        p_an = par_dct['operation']['nominal_electrode_pressure']['anode']['value']
        p_in = [p_ca, p_an]
        dummy = None

        dct_elchem = par_dct['electrochemistry']
        pec = hd.dct_to_nt(dct_elchem, subkey='value')
        dct_press = par_dct['operation']['nominal_electrode_pressure']
        p = hd.dct_to_nt(dct_press, subkey='value')
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
        #obj = None # dummy
        pp_in = flws_clc.partial_pressure_smpl(dummy, pec, opv, T_N,)
        #print('pp_in: ', pp_in)
        pout = gnrl_pwr_clc.op_opt(dummy, pec, T_N, i_ini, imx,
                                p, pp_in,
                                u_mx=uN_cell, ifun=plr_clc.voltage_cell, ini=True)
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
    p2_max = pwr_st_N * p_frc_max

    arr_p_N = np.array([p0_N, p1_N, p2_N,])
    arr_p_max = np.array([p0_max, p1_max, p2_max,])

    p_N = arr_p_N[arr_p_N >0]
    if len(p_N):
        p_N = p_N[-1]
    else:
        print('Nominal cell power could not be determined...')

    p_max = arr_p_max[arr_p_max > 0]
    if len(p_max) >1:
        p_max = p_max[-1]
    else:
        print('Maximum cell power could not be determined...')

    if not p_frc_max:
        p_frc_max = p_max / p_N
        par_dct['operation']['maximum_power_fraction'] = p_frc_max

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
            par_dct['plant']['power_of_plant']['values']['max'] = pwr_plnt_max

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
            par_dct['plant']['power_of_stack']['values']['max'] = pwr_st_max

    if (not n_st) & (not n_cells_st):

        n_cells_st = math.ceil(pwr_st_N / (p_N))
        pwr_st_N = p_N[0] * n_cells_st
        ### check power again
        while pwr_st_N > pwr_st_max:
            #print('----p_N: ', p_N)
            n_cells_st -=1
            pwr_st_N = p_N[0] * n_cells_st

        ### TODO: what about max values???
        n_st = math.ceil(pwr_plnt_N / pwr_st_N)
        pwr_plnt_N = pwr_st_N * n_st
        n_cells_tot = math.ceil(n_st * n_cells_st)

    print('---n_st: ', n_st)
    print('---n_cells_st: ', n_cells_st)
    par_dct['plant']['number_of_stacks']['act'] = n_st
    par_dct['plant']['number_of_cells_in_stack']['act'] = n_cells_st
    par_dct['plant']['number_of_cells_in_plant']['act'] = n_cells_tot
    par_dct['cell']['active_cell_area']['value'] = A_cell
    par_dct['plant']['power_of_stack']['values']['nominal'] = pwr_st_N
    par_dct['plant']['power_of_plant']['values']['nominal'] = pwr_plnt_N

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
