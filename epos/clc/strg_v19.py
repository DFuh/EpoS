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
    storage = xstrg.STRG()

    full_df = pd.concat(df_lst)
    full_df['date'] = pd.to_datetime(full_df.date)
    dates = full_df.date
    full_df = full_df.set_index('date')

    T_in = 313 # Temp of gas treatment
    p_in = 101325 # Pressure of Electrolyzer output // in Pa

    print('tdiff (00): ', dates.iloc[0], dates.iloc[-1])
    print('time-diff (0): ', full_df.date_num.iloc[-1]-full_df.date_num.iloc[0])
    #full_df['date'] = pd.to_datetime(full_df.date)
    print('tdiff (1): ', (dates.iloc[-1]-dates.iloc[0]).total_seconds())
    print('tdiff (2): ', (full_df.index[-1]-full_df.index[0]).total_seconds())
    seconds_tot = (dates.iloc[-1]-dates.iloc[0]).total_seconds()
    t_span=[0,seconds_tot]
    print('t_span: ', t_span)
    # print('full_df: ', full_df.head(4))
    flow_prd = full_df['n_H2_ca'].to_numpy() * obj.pec.M_H2 # Convert mol/s to kg/s

    # TOD: add flow_demand !
    flow_cns = full_df['n_H2_ca'].to_numpy()/3600

    flow_bal = flow_prd-flow_cns
    dm_in = np.where(flow_bal>0, flow_bal, 0)
    dm_out = np.where(flow_bal<0, abs(flow_bal), 0)

    ret, m_dot_act = storage.clc_state_mstp(T_in, p_in, t_span, dm_in, dm_out,
                                        max_step=np.inf, heatloss_wall=True)
    full_df = full_df.assign(rho_strg=ret[0],
                    u_strg=ret[1],
                    T_Strg=ret[2],
                    p_strg=ret[3],
                    t_Strg=ret[4],
                    m_strg=ret[5],
                    U_strg=ret[6],
                    Q_strg=ret[7],
                    w_cmp=ret[8])

    return full_df

def clc_strg_state_dyn():
    '''
    Calculate state of storage dynamically;
    rebound/feedback to EL-simulation
    '''

    return
