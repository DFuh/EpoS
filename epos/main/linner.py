'''
inner loop
-> main calculations
'''

import numpy as np

def subloop(obj, data_in, tnum, time_incr_clc, ):
    (date, t_abs, t_diff,
    T_st, m_c, m_ely, i_cell,
    u_cell, u_an, u_ca, u_dgr,
    P_in, P_act, P_st, P_rct, P_aux,
    p_an, p_ca,
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
        t_diff[m] = time_incr_clc           # Redundant  ?
        t_abs[m] = t_abs[m-1] + t_diff[m]   #

        plr_clc.testf(m, obj.pec)

        ### Calc stack temperature and coolant flowrate
        #T_st[m] =
        #m_c[m] =

        ### clc pressure at electrodes

        ### Calc densities of flows
        #rho_?

        ### Calc auxilliary power consumption/ demand (BoP)
        #P_aux[m] =


        ### Calc power demand of rectifier
        #P_rect[m] =

        ### Maximum power gradient
        #pow_grad

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
        #pp_H2_ca =
        #pp_O2_ca =
        #pp_H2_an =
        #pp_O2_an =


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
                            p_an, p_ca,
                            n_H2_an, n_H2_ca, n_O2_an, n_O2_ca,
                            n_per_H2, n_per_O2,
                            n_H2O,
                            x_H2inO2, x_O2inH2,
                            d_mem,
                            ])

    return data_out
