'''
inner loop
-> main calculations
'''

import numpy as np
from datetime import datetime

def subloop(obj, data_in, tnum, time_incr_clc, ini=False):
    '''
    (date, t_abs, t_diff,
    T_st, m_ely, m_c, i_cell, # Temp. of Stack, massflow coolant, massflow water/electrolyte
    u_cell, u_an, u_ca, u_dgr,
    P_in, P_act, P_st, P_rct, P_aux,
    pp_H2_sep_ca, pp_H2_sep_an, pp_O2_sep_ca, pp_O2_sep_an,
    n_H2_an, n_H2_ca, n_O2_an, n_O2_ca,
    n_per_H2, n_per_O2,
    n_H2O,
    x_H2inO2, x_O2inH2,
    d_mem,                              ) = data_in
    '''
    ntd = data_in
    # --------------------------------------------------------------------------
    #clc modules
    plr, flws, dgr, pwr, thrm, strg, aux = obj.clc_m
    # --------------------------------------------------------------------------
    # pre clc

    sl_no_error = True

    print('linner ntd.date: ', datetime.fromtimestamp(ntd.date[0]/1e9).strftime("%Y%m%d -%H%M%S"))
    ###########################################################################
    # clc loop
    m = 1
    while( m < tnum+1 ) & (sl_no_error):

        ntd.t_diff[m] = (ntd.date[m] - ntd.date[m-1]) /1e9 #time_incr_clc // in s           # Redundant  ?
        if 0 <= m <= 2:
            print('tdiff: ', ntd.t_diff[m])
        if (not ini)or(m>1):
            ntd.t_abs[m] = ntd.t_abs[m-1] + ntd.t_diff[m]   #

        if (m==2)&False:
            plr.testf(m, obj.pec)

        # ==========
        # CAUTION: Stack vs. cell vs. plant level (Power, ...)
        # ->> Voltage, etc.             | on cell-level)
        # ->> m_ely, m_c, P_heat, n_i   | on plant level
        # ==========
        # ntd.T_st[m] = ntd.T_st[m-1]+1
        # ntd.t_abs[m] = ntd.t_abs[m-1]+10
        if True:
            #i_cell = 2
            ntd.m_ely[m-1] = 0.1*obj.av.stckfctr


            ### Calc stack temperature and coolant flowrate
            (ntd.T_st[m], ntd.m_ely[m],
            ntd.m_c[m], P_heat) = thrm.heatbalance(obj, ntd.T_st[m-1],
                                                        ntd.m_ely[m-1], ntd.m_c[m-1],
                                                        ntd.u_cell[m-1], ntd.i_cell[m-1],
                                                        ntd.n_H2_ca[m-1], ntd.n_O2_an[m-1],
                                                        ntd.n_H2O_cns[m-1],
                                                        (ntd.t_abs[m-1], ntd.t_abs[m]),
                                                        ntd=ntd, Tconst=False)# True)
            '''
            ??? split heatbalance in:
                temp-clc
                auxvals (flws)
                temp-ctrl
            ???
            '''
            aux.clc_auxvals(obj, ntd.T_st[m])
            #print('T_st[m]: ', T_st[m])
            #m_c[m] =

            ### clc pressure at electrodes

            ### Calc densities of liquid flows
            ## -->clc auxpars ?
            #rho_?
            #flws.xflws.clc_flws_auxpars(obj, ntd.T_st[m]) #???
            # ---> moved inside heatbalance



            ### Calc auxilliary power consumption/ demand (BoP)
            # based on feed-water supply, gas-dryer (flows of product gas)
            #print(f'm_ely: {ntd.m_ely[m]}|| m_clnt: {ntd.m_c[m]}')
            ntd.P_aux[m] = pwr.clc_pwr_bop(obj, ntd.m_ely[m], ntd.m_c[m],
                                            ntd.n_H2_ca[m-1], P_heat)
            P_avail = ntd.P_in[m] - (ntd.P_aux[m]) # Plant level
            #print(f'ntd.P_in[m]={ntd.P_in[m]} ||ntd.P_aux[m]={ntd.P_aux[m]}')
            #print('P_avail = ', P_avail)
            ###
            #pwr.cntrl_pow_clc(obj, pec, T, i, p, pp, P_avail)

            '''
            CHECK: assignment w.r.t ->> m vs. m-1  (pp, t_diff, ...)
            '''
            # print('P_avail (linner): ', P_avail)
            # print('P_min (linner): ', obj.av.power_stack_min*obj.pplnt.number_of_stacks_act)
            if P_avail > (obj.av.power_stack_min*obj.pplnt.number_of_stacks_act):
                pp = ntd.pp_H2_ca[m-1], ntd.pp_O2_an[m-1], obj.av.pp_H2O
                #print('pp (linner): ', pp)
                (ntd.P_st[m], ntd.P_rct[m],
                ntd.i_cell[m], ntd.u_cell[m]) = pwr.cntrl_pow_clc(obj, obj.pec,
                                                    ntd.T_st[m], ntd.i_cell[m-1],
                                                    obj.p, pp,
                                                    P_avail, ntd.P_st[m-1],
                                                    ntd.u_cell[m-1], ntd.t_diff[m])
                print( 'cntrl_pow_clc: P_St={0}, P_rct={1},i={2},u={3}'.format(ntd.P_st[m], ntd.P_rct[m],
                                                    ntd.i_cell[m], ntd.u_cell[m]))
            else:
                (ntd.P_st[m], ntd.P_rct[m], ntd.i_cell[m], ntd.u_cell[m]) = (0,0,0,0)
            '''
            CAUTION: P_st -> Power of ONE Stack
            '''
            # print('P_act0 = ', ntd.P_act[m])
            # print('P_rct = ', ntd.P_rct[m])
            ntd.P_act[m] = (ntd.P_st[m] * obj.pplnt.number_of_stacks_act
                            + ntd.P_rct[m])
            # print('P_act1 = ', ntd.P_act[m])
            # print('stckfctr = ',obj.av.stckfctr)
            ### Maximum power gradient
            #pow_grad

            n_in = (ntd.n_H2_an[m-1]/obj.av.stckfctr,
                    ntd.n_H2_ca[m-1]/obj.av.stckfctr,
                    ntd.n_O2_an[m-1]/obj.av.stckfctr,
                    ntd.n_O2_ca[m-1]/obj.av.stckfctr)
            c_in = ntd.c_H2_an[m-1], ntd.c_H2_ca[m-1], ntd.c_O2_an[m-1], ntd.c_O2_ca[m-1]

            flws.materialbalance(obj,ntd.T_st[m],  ntd.i_cell[m],
                                          ntd.m_ely[m], obj.p, c_in, n_in,
                                          stf=obj.av.stckfctr,
                                          ntd=ntd, sns=False, m=m)
            #print('flws output n_H2_ca: ', ntd.n_H2_ca)
            #print('flws output x_H2_an: ', ntd.x_H2_an)
            #print('flws output pp_H2_ca: ', ntd.pp_H2_ca)

            if ((any(np.isnan([ntd.T_st[m], ntd.i_cell[m], ntd.u_cell[m]]))) &
                (any(np.isnan([ntd.T_st[m-1], ntd.i_cell[m-1], ntd.u_cell[m-1]])))):
                sl_no_error = False

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
    '''
    data_out = np.array([   date, t_abs, t_diff,
                            T_st, m_ely, m_c, i_cell,
                            u_cell, u_an, u_ca, u_dgr,
                            P_in, P_act, P_st, P_rct, P_aux,
                            pp_H2_sep_ca, pp_H2_sep_an, pp_O2_sep_ca, pp_O2_sep_an,
                            n_H2_an, n_H2_ca, n_O2_an, n_O2_ca,
                            n_per_H2, n_per_O2,
                            n_H2O,
                            x_H2inO2, x_O2inH2,
                            d_mem
                            ])
    '''
    data_out = ntd

    return data_out
