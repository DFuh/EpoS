'''
inner loop
-> main calculations
'''

import numpy as np

def subloop(data_in, tnum):
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
    # pre clc

    ###########################################################################
    # clc loop
    m = 0
    while( m < tnum ):

        data_in[1:,m] = m
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
