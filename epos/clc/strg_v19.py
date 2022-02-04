'''
storage state calculation
'''
import numpy as np
import pandas as pd
import epos.auxf.faux as fx
xstrg = fx.dyn_aux_import(__file__, __name__)

def clc_strg_state_iso(obj, df_lst, ):
    '''
    Calculate state of storage isolated;
    no rebound/feedback to EL-simulation

    df_lst - list (of dataframes w.r.t. different years of clc)
    '''
    obj.logger.info('Run Storage Model (iso)')
    ### ini strg object
    storage = xstrg.STRG(obj)

    full_df = pd.concat(df_lst)
    # print('full_df (strg): ', full_df.head(10))
    full_df['date'] = pd.to_datetime(full_df.date)
    # print(full_df.head(10))
    dates = full_df.date
    full_df = full_df.set_index('date')

    obj.logger.warning('Hardcoded Temp. and Press. (input); Storage Model (iso)')
    T_in = 313 # Temp of gas treatment
    p_in = 101325 # Pressure of Electrolyzer output // in Pa

    # print('dates (in strg-df): ', dates)
    # print('tdiff (00): ', dates.iloc[0], dates.iloc[-1])
    # print('time-diff (0): ', full_df.date_num.iloc[-1]-full_df.date_num.iloc[0])
    #full_df['date'] = pd.to_datetime(full_df.date)
    # print('tdiff (1): ', (dates.iloc[-1]-dates.iloc[0]).total_seconds())
    # print('tdiff (2): ', (full_df.index[-1]-full_df.index[0]).total_seconds())
    seconds_tot = (dates.iloc[-1]-dates.iloc[0]).total_seconds()
    t_span=[0,seconds_tot]
    # print('t_span: ', t_span)
    # print('full_df: ', full_df.head(4))


    bal_df = full_df.copy()
    bal_df['flow_H2_prd'] = bal_df['n_H2_ca'] * obj.pec.M_H2 # Convert mol/s to kg/s
    bal_df['flow_H2_cns'] = bal_df['dm_H2_dmnd'] # kg/s
    flow_cns = bal_df.flow_H2_cns.to_numpy()
    flow_prd = bal_df.flow_H2_prd.to_numpy()
    tdiff = full_df.t_abs.diff().to_numpy()
    tdiff[0] = 0
    # storage.cap_m
    (bal_df['m_strg_sc'],
    bal_df['m_dot_H2_grid'],
    bal_df['m_dot_H2_strg']) =  xstrg.fill_strg(storage.m, storage.m_max, (flow_prd-flow_cns), tdiff)
    # TOD: add flow_demand !
    # flow_cns = full_df['n_H2_ca'].to_numpy()/3600

    # bal_df['flow_H2_prd_cums'] = bal_df.flow_H2_prd.cumsum()
    # bal_df['flow_H2_cns_cums'] = bal_df.flow_H2_cns.cumsum()

    # bal_df['flow_H2_resid'] = bal_df.flow_H2_prd_cums - bal_df.flow_H2_cns_cums
    # bal_df['m_strg_theo'] = bal_df.flow_H2_resid.cumsum()

    #full_df['flow_H2_resid'] = np.where(
    #                    (full_df.flow_H2_cns_cums - full_df.flow_H2_prd_cums) >0,
    #                    full_df.flow_H2_cns, 0)
    # bal_df['flow_H2_ext'] = np.where(bal_df.m_strg_theo<0,-bal_df.flow_H2_resid,0)
    # flow_H2_ext = bal_df['flow_H2_ext'].to_numpy()

    # flw_sum_prd = flow_prd + flow_H2_ext
    # flow_bal = flow_prd - flow_cns
    # dm_in = np.where(flow_bal>0, flow_bal, 0)
    # dm_out = np.where(flow_bal<0, abs(flow_bal), 0)
    dm_in = np.where(bal_df.m_dot_H2_strg>0, bal_df.m_dot_H2_strg, 0)
    dm_out = np.where(bal_df.m_dot_H2_strg<0, abs(bal_df.m_dot_H2_strg), 0)

    ret, m_dot_in_act = storage.clc_state_mstp(T_in, p_in, t_span, dm_in, dm_out,
                                        max_step=np.inf, heatloss_wall=True)
    full_df = full_df.assign(rho_strg=ret[0],
                    u_strg=ret[1], # J/kg
                    T_Strg=ret[2], # K
                    p_strg=ret[3], # Pa
                    t_Strg=ret[4], # h
                    m_strg=ret[5], # kg
                    U_strg=ret[6], # kJ
                    Q_strg=ret[7], # kW
                    P_cmp=ret[8], # kW
                    dm_H2_to_strg_act=m_dot_in_act,   # kg/s
                    dm_H2_to_strg=dm_in,        # kg/s
                    dm_H2_from_strg=dm_out,     # kg/s
                    dm_H2_ext=bal_df['m_dot_H2_grid'].to_numpy()
                    )

    return full_df

def clc_strg_state_dyn():
    '''
    Calculate state of storage dynamically;
    rebound/feedback to EL-simulation
    '''

    return
