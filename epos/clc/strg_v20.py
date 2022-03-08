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
    -> new version of calc sliced by years

    df_lst - list (of dataframes w.r.t. different years of clc)

    '''
    obj.logger.info('Run Storage Model (iso)')
    lssim = len(df_lst)
    lst_out = []
    for n,df in enumerate(df_lst):

        if n == 0: # ini-vals provided by strg-par
            T0_in = None
            p0_in = None
        else:
            T0_in = T0pre # lst_out[n-1].T_strg.iloc[-1]
            p0_in = p0pre # lst_out[n-1].p_strg.iloc[-1]
        obj.logger.info('Start Storage Simu no. %s / %s', str(n+1),str(lssim))
        try:
        #if True:

            df = clc_strg_sngl(obj, df,T0_in=T0_in, p0_in=p0_in)

            T0pre = df.T_strg.iloc[-1]
            p0pre = df.p_strg.iloc[-1]
        except Exception as e:
        #if False:
            # obj.logger.info(' --- Storage Simu no. %s failed', str(n/lssim))
            obj.logger.info(' --- Storage Simu no. %s failed; >> Cause: %s ', str(n), e)
            T0pre = None
            p0pre = None
        obj.logger.info('End Storage Simu no. %s / %s', str(n+1),str(lssim))
        lst_out.append(df)
    return lst_out

def clc_strg_sngl(obj, full_df, T0_in, p0_in):
    '''
    Single run of strg-simu
    '''
    ### ini strg object
    storage = xstrg.STRG(obj, T0=T0_in, p0=p0_in)

    # full_df = pd.concat(df_lst)
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
    # print('tdiff (1): ', (dates.iloc[-1]-dates.iloc[0]))
    # print('tdiff (2): ', (full_df.index[-1]-full_df.index[0]).total_seconds())

    seconds_tot = (dates.iloc[-1]-dates.iloc[0]).total_seconds()
    t_span=[0,seconds_tot]
    # print('t_span: ', t_span)
    # print('full_df: ', full_df.head(4))



    # bal_df = full_df[] #.copy()
    # bal_df['flow_H2_prd'] = bal_df['n_H2_ca'] * obj.pec.M_H2 # Convert mol/s to kg/s
    # bal_df['flow_H2_cns'] = bal_df['dm_H2_dmnd'] # kg/s
    #flow_cns = bal_df.flow_H2_cns.to_numpy()
    # flow_prd = bal_df.flow_H2_prd.to_numpy()*obj.pstrg.fctr_scl_input
    flow_prd = full_df.n_H2_ca.to_numpy() * obj.pec.M_H2*obj.pstrg.fctr_scl_input
    flow_cns = full_df.dm_H2_dmnd.to_numpy()
    # flow_prd = bal_df.flow_H2_prd.to_numpy()*obj.pstrg.fctr_scl_input

    tdiff = full_df.t_abs.diff().to_numpy()
    tdiff[0] = 0
    # storage.cap_m

    # clmns = ['n_H2_ca', 'dm_H2_dmnd']
    l = [col for col in ['m_strg_sc', 'm_dot_H2_grid', 'm_dot_H2_strg'] if col in full_df.columns]
    if len(l) < 3:
        obj.logger.info(' --- Run fill_strg --- ')
        (m_strg_sc,
        m_dot_H2_grid,
        m_dot_H2_strg) =  xstrg.fill_strg(storage.m - storage.m_min, storage.m_max-storage.m_min, (flow_prd-flow_cns), tdiff)
    else:
        obj.logger.info(' --- data from simple calc exists --- ')
        m_strg_sc = full_df.m_strg_sc.to_numpy()
        m_dot_H2_grid = full_df.m_dot_H2_grid.to_numpy()
        m_dot_H2_strg = full_df.m_dot_H2_strg.to_numpy()
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
    dm_in = np.where(m_dot_H2_strg>0, m_dot_H2_strg, 0)
    dm_out = np.where(m_dot_H2_strg<0, abs(m_dot_H2_strg), 0)

    obj.logger.info(' --- start clc_state_mstp --- ')

    max_stp = np.inf
    obj.logger.info(' solve_ivp max_step= %s', str(max_stp))
    ret, m_dot_in_act = storage.clc_state_mstp(T_in, p_in, t_span, dm_in, dm_out,
                                        max_step=max_stp, heatloss_wall=True, simu_obj=obj)
    full_df = full_df.assign(rho_strg=ret[0],
                    u_strg=ret[1], # J/kg
                    T_strg=ret[2], # K
                    p_strg=ret[3], # Pa
                    t_strg=ret[4], # h
                    m_strg=ret[5], # kg
                    U_strg=ret[6], # kJ
                    Q_strg=ret[7], # kW
                    P_cmp=ret[8], # kW
                    #dm_H2_to_strg_act=m_dot_in_act,   # kg/s
                    dm_H2_to_strg=dm_in,        # kg/s
                    dm_H2_from_strg=dm_out,     # kg/s
                    dm_H2_ext=m_dot_H2_grid
                    )
    return full_df

def clc_strg_state_dyn():
    '''
    Calculate state of storage dynamically;
    rebound/feedback to EL-simulation
    '''

    return
