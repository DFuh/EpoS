#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 09:49:29 2022

@author: dafu_res
"""
'''
adopted from elTeco: materialbalance (and adapted to EpoS)
-> remove, when elTeco is an importable package

TODO:
    - handle different c_electr data and columns ?
    - optional: specify additional cost (or emission) component

'''

import numpy as np
import pandas as pd
# class Gases():
#         def __init__(self,M=None,T=None, p=None ):
#             super().__init__()
#             self.rho = self.clc_density(M, T, p )


def clc_density(R,M,T,p):
    # R = 8.314 # // in J/ mol K
    rho = (M*p) / (R * T)
    return rho

def make_matbal_df(obj, df_lst, meda_lst, yr_lst):
    '''
    Parameters
    ----------
    obj
    df_lst
    meda_lst
    yr_lst

    Returns
    -------
    None.
    '''
    matbal_pth = obj.flpth_out_basic + '_matbal.csv'
    # mb_data_lst = []
    # mb_out = None
    mb_out = pd.concat(df_lst, ignore_index=True)
    mb_out.index = yr_lst
    # for num, df_in in enumerate(df_lst): #self.years): #self.pth_lst):

        # df_out = clc_materialbalance(obj, df_in, meda_lst[num]['year'])
        # if isinstance(mb_out, pd.DataFrame):
        #     mb_out = mb_out.append(df_out, ignore_index=True)
        # else:
        #     mb_out = df_out.copy()
        #mb_df.to_csv(mb_pth)
    mb_out.to_csv(matbal_pth, index=False)
    # mb_out = mb_out.set_index('year')
    # dct_mb = mb_out.T.to_dict() # TODO: unefficient!!! dict->df->dict (see:materialbalance)
    return

def make_hires_matbal(obj, df_lst):#, yr_lst):
    '''
    Make/extract high-resolution matbal-data.
    - sample results (df)
    - apply materialbalance
    - write to csv


    Parameters
    ----------
    obj : object of class [elSim]
        Simulation instance.
    df_lst : list
        list of elSim data output one df per year.

    Returns
    -------
    None

    '''
    obj.logger.info('Start hires matbal')
    lst_mbdfs = []
    for i,df in enumerate(df_lst):
        if not type(df.index)==pd.DatetimeIndex:
            obj.logger.info('Adjust df: DatetimeIndex (hires matbal)')
            if 'date' in df.columns:
                df['date'] = pd.to_datetime(df.date)
                df = df.set_index('date')
        val_resample = '7D'
        dfs = df.resample(val_resample) # Slice dataframe according to selected time-period (default: 1week = 7Days)
        obj.logger.info('Resampled df > %s < (hires matbal)', val_resample)
        for [idx,dfi] in dfs:
            dfo = clc_materialbalance(obj, dfi, idx)
            lst_mbdfs.append(dfo)

    ### Make final df
    obj.logger.info('Make final df (hires matbal)')
    mb_out = pd.concat(lst_mbdfs)#, ignore_index=True)
    # mb_out.index = yr_lst
    ### Write to csv
    matbal_pth = obj.flpth_out_basic + '_matbal_hires.csv'
    obj.logger.info('Write  hire-matbal-df (hires matbal)')
    mb_out.to_csv(matbal_pth, index=False)

    return


def clc_materialbalance(obj, df, yr, stats=True, sig_stats=True):
        '''
        calculate amounts of educts/ products
        and operation stats

        Parameters
        ----------
        obj :
        df :
        yr :
        stats : Bool:
            The default is True.
        sig_stats : Bool
            The default is True.
        '''

        #obj.logger.warning('Hardcoded Temp. and Press. (input); teco_matbal')

        '''

        WARNING !

        -> Temp. and Pressure for V_H2 and V_O2 hardcoded !!!
        '''

        T_in = 313 # Temp of gas treatment
        p_in = 101325 # Pressure of Electrolyzer output // in Pa
        obj.av.rho_Hydrogen = clc_density(R=obj.pec.R, M=obj.pec.M_H2, T=T_in, p=p_in)
        obj.av.rho_Oxygen = clc_density(R=obj.pec.R, M=obj.pec.M_O2, T=T_in, p=p_in)
        #df = ?
        #df = df.reset_index()
        #print('df.head: ', df.head())
        #df['dt_s'] = (df.index - df.index.shift(1)).dt.seconds #(df.date-df.date.shift(1)).dt.seconds #total_seconds()

        ### Include energy prices and cO2-emissionfactors


        # df['dt_s'] = (df.date-df.date.shift(1)).dt.seconds #total_seconds()
        df['dt_s'] = df.t_diff
        df['dt_hr'] = df.dt_s / 3600
        m_H2 = np.nansum(df.n_H2_ca * df.dt_s * obj.pec.M_H2) # amount of produced Hydrogen // in kg
        V_H2 = m_H2 / obj.av.rho_Hydrogen
        m_O2 = np.nansum(df.n_O2_an * df.dt_s * obj.pec.M_O2) # amount of produced Oxygen // in kg
        V_O2 = m_O2 / obj.av.rho_Oxygen
        m_H2O = np.nansum( abs(df.n_H2O_cns) * df.dt_s * obj.pec.M_H2O) # amount of consumed Water // in kg

        if 'dm_H2_dmnd' in df.columns:
            m_H2_dmnd = np.nansum(df.dm_H2_dmnd * df.dt_s) # amount of demanded Hydrogen // in kg
        else:
            m_H2_dmnd = 0
        if 'dm_H2_ext' in df.columns:
            m_H2_ext = sum(df.dm_H2_ext * df.dt_s)
        else:
            m_H2_ext = 0
        if 'm_dot_H2_grid' in df.columns:
            m_H2_to_grid_sc = sum(np.where(df.m_dot_H2_grid>=0,df.m_dot_H2_grid,0) * df.dt_s)
            m_H2_frm_grid_sc = sum(np.where(df.m_dot_H2_grid<0,df.m_dot_H2_grid,0) * df.dt_s)
        else:
            m_H2_to_grid_sc = 0
            m_H2_frm_grid_sc = 0
        if 'm_dot_H2_strg' in df.columns:
            m_H2_to_strg_sc = sum(np.where(df.m_dot_H2_strg>=0,df.m_dot_H2_strg,0) * df.dt_s)
            m_H2_frm_strg_sc = sum(np.where(df.m_dot_H2_strg<0,df.m_dot_H2_strg,0) * df.dt_s)
        else:
            m_H2_to_grid_sc = 0
            m_H2_frm_grid_sc = 0
        #if 'm_dot_H2_grid' in df.columns:
        #    m_H2_ext = sum(df.m_dot_H2_grid * df.dt_s)
        #else:
        #    m_H2_ext = 0

        # print(df.P_in.head(5))
        # print(df.P_act.head(5))
        arr_E_util_act = (df.P_act + df.P_aux) * df.dt_hr
        arr_E_util_st = df.P_st * df.dt_hr
        arr_E_util_DC = df.P_st * df.dt_hr * obj.pplnt.number_of_stacks_act
        E_util_act = np.nansum(arr_E_util_act) # amount of utilized energy // in kWh
        E_util_st = np.nansum(arr_E_util_st) # amount of utilized energy in one Stack (DC) // in kWh
        E_util_DC = np.nansum(arr_E_util_DC) # amount of utilized energy in all Stacks (DC) // in kWh
        arr_E_in = df.P_in* df.dt_hr
        E_in = np.nansum(arr_E_in) # amount of available energy from EE // in kWh

        # TODO: Include compressor-efficiency ?
        # eff_cmp = f(P_cmp)
        if 'P_cmp' in df.columns:
            arr_E_cmp = df.P_cmp * df.dt_hr
            E_cmp = np.nansum(arr_E_cmp)
        else:
            E_cmp = 0
            arr_E_cmp = False


        # bsc_par = self.simu_obj.par['basic']


        # cE_util = bsc_par.get('column_name_electricity_costs_E_util', False)
        cE_util = 'c_electr' # obj.prms.get('nm_col_c_electr', False)

        # f_emiss_util = bsc_par.get('column_name_emission_factor_E_util', False)
        f_emiss_util =  'f_emiss_spc' #obj.prms.get('nm_col_f_emiss', False)

        # cE_cmp = bsc_par.get('column_name_electricity_costs_E_cmp', False)
        cE_cmp = cE_util

        # f_emiss_cmp = bsc_par.get('column_name_emission_factor_E_cmp', False)
        f_emiss_cmp = f_emiss_util
        #print('df.columns: ', df.columns)
        #print('col f_emiss_util: ', f_emiss_util)
        #print('col f_emiss_cmp: ', f_emiss_cmp)
        #print('col c_electr (util): ', cE_util)
        #print('col c_electr (cmp): ', cE_cmp)
        CE_util = np.nansum(arr_E_util_act * df[cE_util]) if cE_util in df.columns else 0
        CE_util = CE_util/1e3 # E_util in kWh // c_agora in €/MWh
        CE_cmp = np.nansum(arr_E_cmp * df[cE_cmp]) if cE_cmp in df.columns else 0
        CE_cmp = CE_cmp/1e3 # E_cmp in kWh // c_agora in €/MWh
        emiss_E_util = np.nansum(arr_E_util_act * df[f_emiss_util]) if f_emiss_util in df.columns else 0
        emiss_E_util = emiss_E_util/1e3 # f_em,iss in g/kWh -> emiss_E_util in kg
        emiss_E_cmp = np.nansum(arr_E_cmp * df[f_emiss_cmp]) if f_emiss_cmp in df.columns else 0
        emiss_E_cmp = emiss_E_cmp/1e3 # f_emiss in g/kWh -> emiss_E_cmp in kg


        E_LHV_H2 = m_H2 * obj.pec.LHV_H2_m
        E_HHV_H2 = m_H2 * obj.pec.HHV_H2_m

        P_N_el = obj.prms['parameters_tec_el']['plant']['power_of_plant_act'].get('value', False)
        P_max_el = df.P_act.max()
        P_max_gen = df.P_in.max()

        if stats:
            #TODO: distinguish between P_st and P_act !!!
            t_op_el = np.nansum(np.where(df.P_act >0, 1,0) * df.dt_hr)# operation time of electrolyser
            t_fl_el = E_util_act/P_N_el # None # full load hours of electrolyser
            # print(f't_simu: t_max={df.dt_hr.max()}, Min= {df.dt_hr.min()}')
            t_simu = df.t_abs.max()/3600 - df.t_abs.min()/3600
        else:
            t_op_el = None# operation time of electrolyser
            t_fl_el = None # full load hours of electrolyser

        if sig_stats:
            t_op_gen = np.nansum(np.where(df.P_in >0, 1,0) * df.dt_hr)# operation time of generation
            t_fl_gen = E_in/P_max_gen # full load hours of ee plant(s)
        else:
            t_op_gen = None # operation time of ee plant(s)
            t_fl_gen = None # full load hours of ee plant(s)

        keys = ('year E_HHV_H2 E_LHV_H2 m_H2 V_H2 m_H2_ext m_O2 V_O2 m_H2O m_H2_dmnd'
                +' E_util_act E_util_DC E_util_st E_in E_cmp'
                +' CE_util CE_cmp emiss_E_util emiss_E_cmp'
                +' t_op_el t_op_gen t_fl_el t_fl_gen'
                +' t_simu'
                +' P_N_el P_max_el P_max_gen'
                +' m_H2_to_grid_sc m_H2_frm_grid_sc m_H2_to_strg_sc m_H2_frm_strg_sc')
        vals = [yr, E_HHV_H2, E_LHV_H2, m_H2, V_H2, m_H2_ext, m_O2, V_O2, m_H2O,
                    m_H2_dmnd, E_util_act, E_util_DC, E_util_st, E_in, E_cmp,
                    CE_util, CE_cmp, emiss_E_util, emiss_E_cmp,
                    t_op_el, t_op_gen, t_fl_el, t_fl_gen,
                    t_simu,
                    P_N_el, P_max_el, P_max_gen,
                    m_H2_to_grid_sc, m_H2_frm_grid_sc, m_H2_to_strg_sc, m_H2_frm_strg_sc]
        mb_dct = {}
        for key, val in zip(keys.split(' '), vals):
            mb_dct[key] = [val]
        # print(df.head())
        # print('mb_dict: ', mb_dct)
        mb_df = pd.DataFrame.from_dict(mb_dct)
        return mb_df
