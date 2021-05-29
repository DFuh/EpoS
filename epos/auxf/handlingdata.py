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

def check_sig_data(obj,):
    '''
    read metadata of sig and check, if data are consistent with them
    '''
    sig_metadata = {}
    for sc_key,sc_dct in obj.scen_dict.items(): #e.g.{"simu_90:{...}, simu_81:{...}"}

        # read sig file -> metadata (dict), data (df)

        # sc_dct['sig'] yields name of sig (specified on sig_params)
        sig_mtd, sig_df = rf.read_in_signal_dataset(obj, rel_flpth=obj.sig_par[sc_dct['sig']]['path'],
                                                    search_key=obj.sig_par[sc_dct['sig']]['searchkey_sig_metadata'])#obj.sig_par['searchkey_sig_metadata'])
        print('sIG METADATA: ', sig_mtd)
        # get data properties
        prop_dict = get_properties_df(sig_df)
        if not prop_dict['dt_continous']:
            obj.logger.warning('Time-range in sig-dataset is not continous!')
        del prop_dict['dt_continous']
        # check properties
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

    # get enddate
    # get list of years
    return {'start_date': sd, 'end_date': ed, 'time_incr': tdiff[0], 'dt_continous':result, 'years': years}

# ----------------------------------------------------
def ini_data_output(obj):
    '''
    mk path
    make dir
    make files (wrt years)
    '''
    # make path
    obj.pth_data_out = hf.mirror_output_path(basename=obj.cwd,
                                            ref_pth=obj.prms['refpth_sig_data'], #reference
                                            filepath=obj.prms['relpth_sig_data'], # real sig path
                                            bsc_pth_out=obj.prms['bsc_pth_data_out'],
                                            tday=obj.today_ymd,name=obj.prms['scen_name'])
    pth_out = os.path.join(obj.cwd, obj.pth_data_out)
    #print('--- Pth_out in ini_data_output: ', pth_out)
    # make directory
    hf.mk_dir(full_path=pth_out,
                add_suffix=None,
                no_duplicate=False)

    # ini files for each year
    startdate_0 = obj.prms['metadata_sig']['start_date']
    enddate_0 = obj.prms['metadata_sig']['end_date']
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
        flpth_out = os.path.join(pth_out, str(obj.tag)+nm)
        df0_out, df_key_lst = mk_df_data_output(obj, dates)

        ###add row with units edit: DF, 20201203
        '''
        # -->>> TODO: get valid units !
        '''
        units=[['[1]']*len(df0_out.columns)] # Spaceholder
        '''
        <<<---
        '''
        df_units = pd.DataFrame(data=units,columns=df0_out.columns, index=['unit'])
        df_out = pd.concat([df_units, df0_out])

        ### following line creates output-file
        hf.mk_output_file(obj, yr, n, l, flpth_out, df_out, dates)
        if n == 0: #???
            df0 = df0_out
            lst0_df_keys = df_key_lst
        lst_pths_out.append(flpth_out)
    return df0, lst0_df_keys, lst_pths_out


def mk_df_data_output(obj, dates):
    '''
    initialize output dataframe
    '''
    key_lst = []
    full_lst = []
    #print(inst.s_parameters.data_output['varkeys'])
    for key, val in obj.prms['output_parameters'].items():
        #print(key, '--', val)
        full_lst.append(key)
        if val['store'] & obj.prms['select_stored_values']:
            key_lst.append(key)
    #print('key-list for df0: ', key_lst)
    default_data = [0]*len(full_lst) #np.zeros(len(key_lst))
    #default_data[0] = pd.to_datetime(str(obj.sig.df['Date'].min().year))#'2020-01-01 00:00:00'
    default_data[0] = dates[0] #'2020-01-01 00:00:00'
    #df0     = pd.DataFrame(data=[default_data], columns=full_lst) # full df
    #NT0 = namedtuple('NT0', full_lst)
    #nt0 = NT0(default_data)
    df0     = pd.DataFrame( columns=full_lst) # full df
    #df_out  = pd.DataFrame(data=[default_data], columns=key_lst) # output (store) df

    if not obj.prms['select_stored_values']:
        key_lst = full_lst
    #df0['date'] = '2020-01-01 00:00:00'
    return df0, key_lst#, full_lst


def ini_auxvals(obj, par):
    '''
    initialize dataclass to store/handle auxilliary values
    -> get initial values from par-dict
    '''
    @dataclass
    class AuxVals():
        d_mem = par['electrochemistry'].get('d0_mem', None)

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

        t_ol = 0 # Overload time (counter) // in s
        t_nom = 0
        dRct = 0 # ?

    obj.av = AuxVals()

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
