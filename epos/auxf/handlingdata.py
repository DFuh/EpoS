'''
handling data
'''
import os
import numpy as np
import pandas as pd
from collections import namedtuple
from dataclasses import dataclass


import epos.auxf.readingfiles as rf
import epos.auxf.handlingfiles as hf
from epos.clc import teco_matbal_v20 as tecomb

def check_sig_data(obj,):
    '''
    read metadata of sig and check, if data are consistent with them
    '''
    sig_metadata = {}
    for sc_key,sc_dct in obj.scen_dict.items(): #e.g.{"simu_90:{...}, simu_81:{...}"}

        # read sig file -> metadata (dict), data (df)
        par_sig_in = obj.par_sig[sc_dct['input']]
        pth_sig = hf.mk_abspath(obj, par_sig_in['relpth_sig_data'])
        # sc_dct['sig'] yields name of sig (specified on sig_params)
        (sig_mtd,
        sig_df) = rf.read_in_dataset(obj, pth_sig,
                                search_key=par_sig_in['searchkey_sig_metadata'])#obj.sig_par['searchkey_sig_metadata'])
        print('sIG METADATA: ', sig_mtd)
        # get data properties
        prop_dict = get_properties_df(sig_df)
        if not prop_dict['dt_continous']:
            obj.logger.warning('Time-range in sig-dataset is not continous!')
        del prop_dict['dt_continous']
        # check properties
        if sig_mtd is not None:
            for key,val in prop_dict.items():
                #print('Key, val: ', key, val)
                if key not in sig_mtd.keys():
                    #print('Keys in metadata: ', sig_mtd.keys())
                    #print('Key, Val of prop_dict: ', key, ',', val)
                    obj.logger.info('Updating sig-metadata: %s : %s', str(key), str(val))
                    sig_mtd[key] = val
                elif val != sig_mtd[key]:
                    #print('Key, Val of prop_dict: ', key, ',', val)
                    #print('Key, Val of ,metadata_dict: ', key, ',', sig_mtd[key])
                    obj.logger.info('Updating sig-metadata: %s : %s', str(key), str(val))
                    sig_mtd[key] = val
            sig_metadata[sc_key]=sig_mtd
        else:
            sig_metadata[sc_key]=None
    return sig_metadata

def get_properties_df(df_in):
    '''
    extract properties of given df
    -> startdate (first date of datarows)
    -> enddate (last date of datarows)
    -> timeincrement
    -> years contained in df
    '''

    df_c = df_in.copy()
    sd      = df_c.Date.min().strftime("%Y-%m-%d %H:%M:%S") # startdate
    ed      = df_c.Date.max().strftime("%Y-%m-%d %H:%M:%S") # enddate
    tdiff = list(df_c.Date.diff()[1:] / np.timedelta64(1, 's')) #df_c.Date.diff()
    result = all(element == tdiff[0] for element in tdiff)

    years = df_c.Date.dt.year.drop_duplicates().to_list() # list of years in sig-df
    #years = df_c.Date.dt.year.unique() # list of years in sig-df

    # get enddate
    # get list of years
    return {'start_date': sd, 'end_date': ed, 'time_incr': int(tdiff[0]), 'dt_continous':result, 'years': years}


# ----------------------------------------------------
def mk_pth_out():

    return

def ini_data_output(obj):
    '''
    mk path
    make dir
    make files (wrt years)
    '''


    # make path
    '''
    edit for storing data in external directory

    '''
    # relpth_data_out = hf.mirror_output_path(basename=obj.cwd,
    #                                        ref_pth=obj.prms['refpth_input_data'], #reference
    #                                        filepath=obj.prms['relpth_sig_data'], # real sig path
    #                                        bsc_pth_out=obj.prms['bsc_pth_data_out'],
    #                                        tday=obj.today_ymd,name=obj.prms['scen_name'])

    #if obj.prms['external_dir_data_output']:
    #    obj.pth_data_out = os.path.join(obj.prms['external_dir_data_output'], '')#None)#relpth_data_out)

    #else:
    #    obj.pth_data_out = os.path.join(obj.cwd, 'data', relpth_data_out)
    #obj.pth_data_out = os.path.join('~/calc_out', relpth_data_out)


    # make directory
    hf.mk_dir(full_path=obj.pth_data_out,
                add_suffix=None,
                no_duplicate=False)

    # ini files for each year
    startdate_0 = pd.Timestamp(obj.prms['metadata_sig']['start_date']).strftime('%Y-%m-%d %H:%M:%S')
    enddate_0 = pd.Timestamp(obj.prms['metadata_sig']['end_date']).strftime('%Y-%m-%d %H:%M:%S')
    l = len(obj.prms['metadata_sig']['years'])
    lst_pths_out = []
    for n,yr in enumerate(obj.prms['metadata_sig']['years']):
        sd = pd.Timestamp(str(yr)).strftime('%Y-%m-%d %H:%M:%S')
        dt = pd.Timedelta(str(obj.prms['metadata_sig']['time_incr'])+'s')
        ed = (pd.Timestamp(str(yr+1))-dt).strftime('%Y-%m-%d %H:%M:%S')
        if sd < startdate_0:
            sd = startdate_0
        if ed > enddate_0:
            ed = enddate_0
        dates = sd,dt,ed
        #for n,yr in enumerate(simu_inst.sig.years):
        if n <10:
            nm = f'__results_{yr}_0{n}.csv'
        else:
            nm = f'__results_{yr}_{n}.csv'

        #obj.flpth_data_out = obj.path_data_out+'/'+str(obj.tag) + nm
        flpth_out = os.path.join(obj.pth_data_out, str(obj.tag)+nm)
        obj.flpth_out_basic = os.path.join(obj.pth_data_out, str(obj.tag))
        df0_out, df0_cd, df_key_lst, unit_lst = mk_df_data_output(obj, dates)

        ###add row with units edit: DF, 20201203
        '''
        # -->>> TODO: get valid units !
        '''
        # units=[['[1]']*len(df0_out.columns)] # Spaceholder
        '''
        <<<---
        '''
        if False:
            unit_lst[0] = '>units<'
            df_units = pd.DataFrame(data=[unit_lst],columns=df0_out.columns, index=['unit'])



        if n == 0: #
            df0 = df0_out.copy()
            lst0_df_keys = df_key_lst
            # df_wr = df0_cd.copy()
        else:
            pass
        df_wr = df0_out.copy()
        df_out = df_wr # pd.concat([df_units, df_wr])        # df_wr])
        ### following line creates output-file
        hf.mk_output_file(obj, yr, n, l, flpth_out, df_out, dates)
        lst_pths_out.append(flpth_out)
    #print('df0: ', df0.head())
    return df0_cd, lst0_df_keys, lst_pths_out

def slice_df_by_years(df):
    df['years'] = df.index.year
    df_lst = []
    for key, dfi in df.groupby(by='years'):
        df_lst.append(dfi)
    return df_lst


def mk_df_data_output(obj, dates):
    '''
    initialize output dataframe
    '''
    key_lst = []
    full_lst = []
    unit_lst = []
    #print(inst.s_parameters.data_output['varkeys'])
    default_data = [0]*len(obj.prms['output_parameters']) #np.zeros(len(key_lst))
    for j,(key, val) in enumerate(obj.prms['output_parameters'].items()):
        #print(key, '--', val)
        full_lst.append(key)
        if val['store'] & obj.prms['select_stored_values']:
            key_lst.append(key)
        default_data[j] = val.get("initial_value", 0)
        unit_lst.append(val.get("unit", 'n/a'))
    #print('key-list for df0: ', key_lst)
    # default_data = [0]*len(full_lst) #np.zeros(len(key_lst))
    #default_data[0] = pd.to_datetime(str(obj.sig.df['Date'].min().year))#'2020-01-01 00:00:00'
    default_data[0] = dates[0] #'2020-01-01 00:00:00'
    # print('default_data: ', default_data)
    #df0     = pd.DataFrame(data=[default_data], columns=full_lst) # full df
    #NT0 = namedtuple('NT0', full_lst)
    #nt0 = NT0(default_data)
    df0     = pd.DataFrame( columns=full_lst) # full df
    df0_cd     = pd.DataFrame( data=[default_data], columns=full_lst) # full df
    #df_out  = pd.DataFrame(data=[default_data], columns=key_lst) # output (store) df

    if not obj.prms['select_stored_values']:
        key_lst = full_lst
    #df0['date'] = '2020-01-01 00:00:00'
    return df0, df0_cd, key_lst, unit_lst#, full_lst


def ini_auxvals(obj, par):
    '''
    initialize dataclass to store/handle auxilliary values
    -> get initial values from par-dict
    '''
    @dataclass
    class AuxVals():
        print('...Ini AuxVals()')
        d_mem = par['electrochemistry'].get('d0_mem', {}).get('value',None)
        c0_memrxn = par['electrochemistry'].get('c0_memrxn',{}).get('value',np.zeros(9)) # initial concentrations for membrane degradation reaction
        #if obj.pplnt.power_gradient_stack_max_pos:
        #    dPdt_p = obj.pplnt.power_gradient_stack_max_pos
        dPdt_p = obj.pplnt.power_of_stack_nominal * obj.pop.power_gradient_max_positive
        dPdt_n = obj.pplnt.power_of_stack_nominal * obj.pop.power_gradient_max_negative
        power_stack_min = obj.pplnt.power_of_stack_nominal * obj.pop.power_fraction_min
        #if obj.pplnt.power_gradient_stack_max_neg:
        #    dPdt_n = obj.pplnt.power_gradient_stack_max_neg
        if obj.pcll.voltage_gradient_max_pos:
            dudt_p = obj.pcll.voltage_gradient_max_pos
        if obj.pcll.voltage_gradient_max_neg:
            dudt_n = obj.pcll.voltage_gradient_max_neg

        ### for plr clc
        cnc_O2_mem_an = 0 #c_O2_mem_an
        cnc_H2_mem_ca = 0 #c_H2_mem_ca

        t_diff = 0
        t_op = 0
        t_ol = 0 # Overload time (counter) // in s
        t_uni = 0 # Uninterrupted time of operation // in s
        t_nom = 0

        dRct = 1 # // in 1 !!!

        dU_dgr_abs = 0

    @dataclass
    class PressureVals():
        print('...Ini PressureVals()')
        anode_nom = par['operation']['nominal_electrode_pressure'].get('anode', {}).get('value',None)
        cathode_nom = par['operation']['nominal_electrode_pressure'].get('cathode', {}).get('value',None)
        anode = anode_nom
        cathode = cathode_nom
        pp_O2_mem_ca = par['operation']['initial_pressure_values'].get('pp_O2_ca', {}).get('value',0)
        pp_O2_mem_an = par['operation']['initial_pressure_values'].get('pp_O2_an', {}).get('value',0)
        pp_H2_mem_ca = par['operation']['initial_pressure_values'].get('pp_H2_ca', {}).get('value',0)
        pp_H2_mem_an = par['operation']['initial_pressure_values'].get('pp_H2_an', {}).get('value',0)
    obj.av = AuxVals()
    obj.p = PressureVals()

    return

def ini_snsvals(obj,):
    '''
    initialize dataclass to store/handle sensitivity factors/values
    -> get initial values from par-dict
    '''
    @dataclass
    class AuxVals():
        pass

    obj.snsfctr = AuxVals()

    return


def dct_to_nt(dct_in, subkey=None):
    '''
    namedtuple from dict
    '''
    # TODO: implement check of units!

    ndct = {}
    for key, sdct in dct_in.items():
        if isinstance(sdct, dict):
            if subkey and (subkey in sdct):
                ndct[key] = sdct[subkey]
    NT = namedtuple('NT',ndct)
    #NT = namedtuple('NT',par['electrochemistry'])
    #self.pec = NT(**par['electrochemistry'])
    nt = NT(**ndct)
    return nt



def setup_refvals_nt(ref_dict, testmode):
    '''
    Setup refvals for testing basic calculations

    Returns
    -------
    Namedtuple containing reference values
    '''
    dct = ref_dict.get(testmode, None)
    #if not isinstance(), dict):
        #raise Exception('No refvals found for mode: ', mode)
    #    nt = None
    if dct:
        new_dict = {}
        for key, val in dct.items():#ref_dict[testmode].items():
            for skey, sval in val.items():
                if 'val' in skey:
                    new_dict[skey] = sval
                elif (not 'unit' in skey) or ('comm' in skey) or ('val' in skey):
                    new_key = key+'_'+skey
                    new_dict[new_key] = sval

        NT = namedtuple('NT',new_dict)
        nt = NT(**new_dict)
        return nt
    return dct


def read_data(obj,):
    # pths = [val for key,val in obj.prms.items()
    #          if (('relpth_' in key) and (obj.prms['refpth_input_data'] in val))]
    df_lst = []
    spec_dct = {}
    pths = []
    colnms = []
    skeys=[]
    nms = []
    # print(obj.prms['refpth_input_data'])
    for key,val in obj.prms.items():
        if ((val is not None) and ('relpth_' in key)
            #and (str(obj.prms['refpth_input_data']) in str(val))
            and ('bumptest_data' not in key)
            and ('sig_par' not in key)):
            pths.append(hf.mk_abspath(obj,pth=val))
            nm = key.replace('relpth_', '').replace('_data','')
            colnms.append(obj.prms['nm_col_'+nm])
            skeys.append(obj.prms.get('searchkey_'+nm+'_metadata',None))
            nms.append(nm)

    for pth, colnm, skey, nm in zip(pths, colnms, skeys, nms):
        specs, df = rf.read_in_dataset(obj, abspth_to_fl=pth,
                                    search_key=skey)
        print(nm, colnm)
        if not colnm in df.columns:
            print('cols in df: ', df.columns)
            print('given col-name: ', colnm)
            colnm = input('Please insert valid column-name to be used:')
            obj.prms['nm_col'+nm] = colnm
        spec_dct['metadata_'+nm] = specs
        df_lst.append(df)

    # for i,dfi in enumerate(df_lst):
        # print(nms[i])
        # print(dfi.head(3))
        # print('---')
    idx = nms.index('sig')
    df0 = df_lst[idx].copy()
    del df_lst[idx]

    for dfi,nm in zip(df_lst, nms):
        df0 = merge_and_fill(dfi, df0, indi=False)
    return specs, df0

def extract_data(obj, df_lst, meda_lst,
                pths_orig_data=[],
                extr_keys=[], extr_nms=[], basic_pth=None):
    '''
    extract data from full results-df
        loop trough df list (year based)
        and apply:
            extraction
            matbal

    '''
    matbal_df_lst = []
    yr_lst = []

    if any(isinstance(i, list) for i in extr_keys):
        nested_keys = True
        N = len(extr_nms)
        extr_df_lst = [ [] for _ in range(N) ]
        extr_meda_lst = [ [] for _ in range(N) ]
        extr_pth_lst = [ [] for _ in range(N) ]
        # print('extr_df_lst = ', extr_df_lst)

    else:
        nested_keys = False
        extr_df_lst = []
        extr_meda_lst = []
        extr_pth_lst = []

    for j, (df, meda) in enumerate(zip(df_lst, meda_lst)):
        ### matbal
        df_mb_out = tecomb.clc_materialbalance(obj, df, meda['year'])
        matbal_df_lst.append(df_mb_out)
        yr_lst.append(meda['year'])

        ### extraction
        if nested_keys:
            for k, key_lst in enumerate(extr_keys):
                key_lst = [key for key in key_lst if key in list(df.columns)]
                # print(f'k = {k} // key_lst = ', key_lst)
                extr_df_lst[k].append(df[key_lst].copy())
                mdi = meda.copy()
                mdi['extraction_for'] = '|'.join(key_lst)
                nm = '_extr_'+ extr_nms[k]+'_' #+''.join(keys_lst)
                extr_meda_lst[k].append(mdi)
                pth = pths_orig_data[j].replace('results', nm)
                extr_pth_lst[k].append(pth)
        elif not nested_keys:
            meda['extraction_for'] = '|'.join(extr_keys)
            extr_meda_lst.append(meda)
            nm = '_extr_'+ ''.join(extr_nms)+'_'
            # nm = extr_nms[0]+'_'+''.join(extr_keys)
            extr_pth_lst.append(pths_orig_data[j].replace('results', nm))

            ### make sure, keys are in df
            key_lst = [key for key in extr_keys if key in list(df.columns)]

            extr_df_lst.append(df[key_lst].copy())

        # print('Final extr_pth_lst: ', extr_pth_lst)

    return matbal_df_lst, yr_lst, extr_df_lst, extr_meda_lst, extr_pth_lst

def add_teco_data_to_df(obj, df_lst_in, metad_lst):
    '''
    add data (electricity costs, emission factors etc.) to simu-df

    '''

    ### read add_data
    add_df0 = rf.read_in_dataset(obj, rel_flpth=self.prms['relpth_c_electr_data'],
                                search_key=self.prms['searchkey_c_electr_metadata'])
    add_df0 = [self.prms['relpth_c_electr_data']]
    add_df1 = rf.read_in_dataset(obj, rel_flpth=self.prms['relpth_f_emiss_data'],
                                search_key=self.prms['searchkey_f_emiss_metadata'])

    df_lst_out = []

    for meda, dfi in zip(metad_lst, df_lst):
        meda['year']

    return

def merge_and_fill(data_df, ref_df=None, freq='10s', indi=True):
    '''
    merge dataframes by date and fill nan-values
    '''


    # df['Date'] = pd.to_datetime(df.Date, format='%Y-%m-%d %H:%M:%S')
    if not 'Date' in data_df.columns:
        print(f'...data-df does not contain column: Date')
        print(f'Columns ins current df: {df.columns}')
        cnm = input('Insert column name to be used: ')
        if cnm:
            data_df.rename(columns = {cnm:'Date'}, inplace = True)
        else:
            print('Could not process df finally... abort...')
            return None
    data_df['Date'] = pd.to_datetime(data_df.Date, format='%Y-%m-%d %H:%M:%S')
    # data_df = data_df.set_index('Date')

    l_cols= 1
    if ref_df is not None:
        l_periods=len(ref_df)
        if 'Date' in ref_df.columns:
            ref_df['Date'] = pd.to_datetime(ref_df.Date, format='%Y-%m-%d %H:%M:%S')
            ref_df.set_index('Date', inplace=True)
        freq = pd.infer_freq(ref_df.index)
    else:
        ref_date = pd.Series(pd.date_range(start_date, periods=l_periods, freq=freq))
        # date2 = pd.Series(pd.date_range('2012-1-1 12:00:00', periods=14, freq='5s'))
        ref_data = np.zeros((l_cols, l_periods))
        # data2 = np.array([np.linspace(30,36,14),np.linspace(40,46,14)])

        ref_df = pd.DataFrame(data=ref_data.T, index=ref_date)
        ref_df = ref_df.reset_index(drop=True)
        ref_cols = ['refnm_']*l_cols
        ref_cols=['refnm_'+str(i) for i in range(l_cols)]
        ref_df.columns=['Date'] + ref_cols

    #changed_df = exdf.chg_resolution(df, tdiff=5, dt_column='Date', end_date='2012-01-01 12:01:00',
    #                                 start_date='2012-01-01 12:00:00', posoffset=None,
    #                                 negoffset=None, filling='ffill')
    # dt_column='Date', end_date='2012-01-01 12:01:00', start_date='2012-01-01 12:00:00')

    ndf2 = pd.merge(ref_df, data_df, on='Date', how='outer', indicator=indi)# True)
    ndf2 = ndf2.sort_values(by='Date')
    ndf2 = ndf2.ffill(axis = 0)
    ndf2.reset_index(drop=True, inplace=True)
    return ndf2
