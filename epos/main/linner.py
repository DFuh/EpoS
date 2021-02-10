'''
inner loop
-> main calculations
'''

import numpy as np

def subloop(obj, data_in, tnum, time_incr_clc, ):
    (date, t_abs, t_diff,
    T_st, m_c, m_ely, i_cell, # Temp. of Stack, massflow coolant, massflow water/electrolyte
    u_cell, u_an, u_ca, u_dgr,
    P_in, P_act, P_st, P_rct, P_aux,
    pp_H2_sep_ca, pp_H2_sep_an, pp_O2_sep_ca, pp_O2_sep_an,
    n_H2_an, n_H2_ca, n_O2_an, n_O2_ca,
    n_per_H2, n_per_O2,
    n_H2O,
    x_H2inO2, x_O2inH2,
    d_mem,                              ) = data_in
    # --------------------------------------------------------------------------
    #clc modules
    plr_clc, flws_clc, dgr_clc, pwr_clc, thrm_clc = obj.clc_m
    # --------------------------------------------------------------------------
    # pre clc

    sl_no_error = True
    ###########################################################################
    # clc loop
    m = 1
    while( m < tnum ) & (sl_no_error):
        t_diff[m] = date[m] - date[m-1] #time_incr_clc           # Redundant  ?
        t_abs[m] = t_abs[m-1] + t_diff[m]   #

        if m==2:
            plr_clc.testf(m, obj.pec)

        #i_cell = 2
        ### Calc stack temperature and coolant flowrate
        #T_st[m], m_c[m], m_ely[m] = thrm_clc.heatbalance(obj,  Tconst=True)
        #print('T[m]: ', T[m])
        #m_c[m] =

        ### clc pressure at electrodes

        ### Calc densities of flows
        #rho_?

        ### Calc auxilliary power consumption/ demand (BoP)
        # based on feed-water supply, gas-dryer (flows of product gas)
        #P_aux[m] =

        ### Calc power demand of rectifier
        #P_rect[m] =

        ### Maximum power gradient
        #pow_grad
        '''
        testtuple = flws_clc.materialbalance(obj, T_st[m],  i_cell[m], m_c[m], p_an[m], p_ca[m])
        print('testtuple output a: ', testtuple.n_H2_out_ca)
        print('testtuple output b: ', testtuple.x_H2_out_ca)
        print('-test auxvals: ', testtuple.pp_H2_mem_ca)
        '''
        ### Absolute pressure at electrodes
        #p_ca[m]
        #p_an[m]


        ### Partial pressures of species  ?
        #pp_H2_ca =
        #pp_O2_ca =
        #pp_H2_an =
        #pp_O2_an =
        #pp_H2O =

        ### Concentration of species
        #c_H2_ca =
        #c_O2_ca =
        #c_H2_an =
        #c_O2_an =


        ### Optimal power point + u->i
        #i[m] =
        #u_cell[m] =
        #u_an[m] =
        #u_ca[m] =

        ### ctrl: plnt on/off?
        # coolant pump / heater running
        #
        # plnt_running = True/ False


        ### Actual stack power (cons)
        #P_st[m] =

        #P_frc ??? // P_pe ??? (in old version: m_pwr_v29 |line 244+256)

        ### Actual plat power (cons)
        #P_act[m] =

        ### mass balance


        ### cell ageing



        m += 1
    ###########################################################################
    # --------------------------------------------------------------------------
    data_out = np.array([   date, t_abs, t_diff,
                            T_st, m_c, m_ely, i_cell,
                            u_cell, u_an, u_ca, u_dgr,
                            P_in, P_act, P_st, P_rct, P_aux,
                            pp_H2_sep_ca, pp_H2_sep_an, pp_O2_sep_ca, pp_O2_sep_an,
                            n_H2_an, n_H2_ca, n_O2_an, n_O2_ca,
                            n_per_H2, n_per_O2,
                            n_H2O,
                            x_H2inO2, x_O2inH2,
                            d_mem
                            ])

    return data_out
