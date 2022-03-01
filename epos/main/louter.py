'''
outer loop
'''
import numpy as np
import pandas as pd
import traceback
import time
from collections import namedtuple

import epos.auxf.faux as fx
import epos.main.aux_louter as xlo

from epos.main import linner


#from importlib import import_module as impm

# TODO: implement output selecttion (reducing data) through slice of columns-list

print(__name__, ' imported.')

def mainloop(obj, ):
    '''
    obj         -> simu instance
    clc_m       -> modules for calc
    par         -> dictenary containing all parameters
    par_elch    -> namedtuple with parameters of electrochemistry
    '''

    obj.logger.info('Ini mainloop')

    #clc_m = obj.calc_modules # | tuple
    ### call in simu.__init__() does not work with Pool.map(): pickling error
    plr_clc, flws_clc, dgr_clc, pwr_clc, thrm_clc, strg_clc, aux_clc = obj.clc_m

    # power input data
    #print('+++ Input df:', simu_inst.sig.df)
    #print('+++ Input df types:', simu_inst.sig.df.dtypes)
    input_df = obj.data_input.copy()
    # dmnd_df = obj.data_H2dmnd.copy()
    input_df['Date'] = pd.to_datetime(input_df.Date)
    # input_df['Date_new'] = pd.date_range(start=input_df.Date.min(),end=input_df.Date.max(), periods=len(input_df))
    # input_df = input_df.set_index('Date_new')
    # dmnd_df = dmnd_df.set_index('Date')

    ###

    #input_df = input_df.loc[pd.to_datetime(simu_inst.s_parameters.starttime), pd.to_datetime(simu_inst.s_parameters.stoptime)].copy()
    #print(pd.to_datetime(obj.prms['sig_metadata']['start_date']))
    input_df = input_df.set_index('Date')

    # print('(louter) input-df (0): ', input_df)
    # print('obj.prms[date_start]: ', obj.prms['date_start'])
    # print('obj.prms[date_end]: ', obj.prms['date_end'])
    # print('obj.scn_setup: ', obj.scn_setup)
    if (not obj.scn_setup):
        #input_df = input_df.loc[pd.to_datetime(simu_inst.s_parameters.starttime): pd.to_datetime(simu_inst.s_parameters.stoptime)]
        input_df = input_df.loc[obj.sd:obj.ed]
    date_in   = input_df.index
    #print('Slicing dates: ', obj.sd, obj.ed)
    #print('Data Input (head, louter): ', input_df.head(5))
    #print('Data Input (tail, louter): ', input_df.tail(5))
    ### ini dict for applying round() to output-df
    rnd_dct = {}
    for key, val in obj.prms['output_parameters'].items():
        dec_val = val.get('round_dec', None)
        if dec_val is not None:
            rnd_dct[key] = dec_val

    '''
    else:
        sd_sig = None
        ed_sig = None
    sd_par = pd.to_datetime(obj.prms['starttime'])

    elif (not obj.scn_setup)
    '''
    # print('(louter) input-df: ', input_df)
    # dmnd_df = dmnd_df.loc[sd:ed]
    # print('Input-DF: ', input_df.head(10))
    #input_df['Date'] = pd.to_datetime(input_df.Date)
    #sidx = input_df[input_df.Date==sd].index[0]
    #eidx = input_df[input_df.Date==ed].index[0]
    #print('sidx: ', sidx)
    #print('eidx: ', eidx)
    #input_df = input_df.iloc[sidx:eidx]
    #print('+++ Sliced df:', input_df)

    # date_in   = input_df.index#.to_numpy()
    #date_in = input_df.Date
    #print('+++ ', __name__,'Date_in: ', date_in)
    # print(dmnd_df.head(4))
    # print(dmnd_df.columns)
    # print(obj.prms['nm_col_H2dmnd'])
    # print(dmnd_df['dm_H2_dmnd'])


    power_in = input_df[obj.prms['nm_col_sig']].to_numpy() # sig-input df -> to np.array
    pow_idx = obj.df0.columns.get_loc('P_in')-1 # TODO: Hardcoded !!

    fctr_scl_sig = obj.prms.get('fctr_scl_sig',False)
    if fctr_scl_sig:
        power_in = power_in*fctr_scl_sig
    if obj.prms.get('nm_col_H2dmnd', False):
        if obj.prms['nm_col_H2dmnd'] in input_df.columns:
            H2dmnd_in = input_df[obj.prms['nm_col_H2dmnd']].to_numpy() # sig-input df -> to np.array
            dmnd_idx = obj.df0.columns.get_loc(obj.prms['nm_col_H2dmnd'])-1 # TODO: Hardcoded !!
            fctr_scl_H2dmnd = obj.prms.get('fctr_scl_H2dmnd',False)
            if fctr_scl_H2dmnd:
                H2dmnd_in = H2dmnd_in * fctr_scl_H2dmnd
        else:
            dmnd_idx=None
    if obj.prms.get('nm_col_c_electr', False):
        if obj.prms['nm_col_c_electr'] in input_df.columns:
            c_electr_in = input_df[obj.prms['nm_col_c_electr']].to_numpy() # sig-input df -> to np.array
            c_electr_idx = obj.df0.columns.get_loc('c_electr')-1 # TODO: Hardcoded !!
        else:
            c_electr_idx=None
    if obj.prms.get('nm_col_f_emiss',False):
        if obj.prms['nm_col_f_emiss'] in input_df.columns:
            f_emiss_in = input_df[obj.prms['nm_col_f_emiss']].to_numpy() # sig-input df -> to np.array
            f_emiss_idx = obj.df0.columns.get_loc('f_emiss_spc')-1 # TODO: Hardcoded !!
        else:
            f_emiss_idx=None
    #print('input df: ', input_df.head())
    #print('input df (dmnd): ', input_df[obj.prms['nm_col_H2dmnd']].head())
    #print('input df (c_electr): ', input_df[obj.prms['nm_col_c_electr']].head())
    #print('input df: (f_emiss)', input_df[obj.prms['nm_col_f_emiss']].head())

    # print('pow_idx: ', pow_idx)
    # length of input-power-df
    len_df_pin          = len(power_in)



    # auxilliary variables
    time_incr_clc   = int(obj.prms['time_incr_clc'])

    if obj.scn_setup:
        time_incr_act = 600
        '''
        CAUTION: HARDCODED !
        '''
    else:

        ### time_incr sig
        if obj.prms['metadata_sig'] is not None:
            time_incr_sig_par = int(obj.prms['metadata_sig']['time_incr'])

        else:
            time_incr_sig_par = None

        # print('input_df.index =', type(input_df.index) )
        # input_df['dt'] = (input_df.index - input_df.index.shift(-1)).seconds

        # input_df['dt'] = input_df['Date'].diff().dt.total_seconds()

        # print('input_df: ', input_df.head(20))
        # print('len(input_df) = ', len(input_df))
        # time_incr_sig_act = (input_df.index - input_df.index.shift(-1)).seconds[100]
        '''
        time_incr_sig_act = input_df['dt'].iloc[20]
        if time_incr_sig_act != time_incr_sig_par:
            print('Deviation in time incr! -> par = {time_incr_sig_par} // act = {time_incr_sig_act}')
            time_incr_act = time_incr_sig_act


            time_incr_sig_act_max = input_df['dt'].max()
            time_incr_sig_act = input_df['dt'].iloc[20]
            time_incr_sig_act_min = input_df['dt'].min()
            print('Deviation in time incr! -> par = {time_incr_sig_par} // act = {time_incr_sig_act}')
            print('Time_incr_act max = ', time_incr_sig_act_max)
            print('Time_incr_act min = ', time_incr_sig_act_min)
            print('[0] -> time_incr_par = ', time_incr_sig_par)
            print('[1] -> time_incr_act = ', time_incr_sig_act)
            slct_ti = input('Which time incr is correct ?')
            time_incr_act = [time_incr_sig_par, time_incr_sig_act][int(slct_ti)]

            obj.prms['metadata_sig']['time_incr'] = int(time_incr_act)
        '''
        time_incr_act = time_incr_sig_par
    tnum            = int(int(time_incr_act) /  time_incr_clc) # number of inner loops



    # print('tnum= ', tnum)
    no_error        = True # initial value
    #lst_full_cols   = obj.df0.columns
    #print('l-cols: ', len(obj.df0.columns))
    clmns =  obj.df0.columns.to_list()

    del clmns[1]
    NT0 = namedtuple('NT0', clmns)
    arr_zeros       = np.zeros((len(clmns), tnum+1)) # default array, shape: length of df0.columns, tnum+1 (initial vals)
    #col_idxs        = xlo.get_col_indexes(obj.df0, lst_full_cols) # make list of indexesfor selected values from clc_array

    #TODO: check, if array-length is consistent with df-length


    # initial values
    #arr_data_in = simu_inst.df[0].to_numpy()
    # set new input data
    data_clc_in         = arr_zeros.copy() # Surplus?
    #simu_inst.df0.iloc[-1][0] = pd.to_datetime(simu_inst.df0.iloc[-1][0])
    data_clc_in[1:,0]    = obj.df0.iloc[-1][1:].to_numpy()[1:] # np array
    if obj.scn_setup:
        df_fin = obj.df0.copy()
        # print('df_ini: ', df_fin)
    #data_clc_in[0,0] = 000

    ### set initial Values
    #for key, val in obj.prms
    #idx_T = obj.df0.columns.get_loc()
    #data_clc_in[,0] =

    ### ini storage model for simple calc
    if obj.prms.get('storage_clc_smpl',False) and not obj.scn_setup:
        sc_Strg = strg_clc.xstrg.STRG(obj) #, T0=T0_in, p0=p0_in)
        obj.av.sc_strg_m = sc_Strg.m - sc_Strg.m_min
        obj.av.sc_strg_m_max = sc_Strg.m_max-sc_Strg.m_min

    obj.av.fctr_n_c_A_abs = obj.pplnt.number_of_cells_in_plant_act * obj.pcll.active_cell_area
    obj.av.fctr_n_c_A_st = obj.pplnt.number_of_cells_in_stack_act * obj.pcll.active_cell_area

    ### counter for P_low
    obj.av.cnt_P_lo = 0

    # call inner loop
    obj.logger.info('Starting mainloop ... ')
    obj.logger.info('Degradation enabled: %s', str(obj.av.enable_dgr))
    t0 = time.time()
    idx = 0 # Years-index
    k = 0
    frc_diff = 0
    while( no_error & (k < len_df_pin)):
        frc = round(k/len_df_pin,3)
        #frc_diff += frc
        #if frc_diff >0.05:
        print(f"Progress (l_outer): {frc} %", end="\r")
        #frc_diff = 0
        # datetime column
        #data_clc_in[0,:] = pd.date_range(start=date_in[k], end=date_in[k+1], freq=str(time_incr_clc)+'s') # daterange
        dat_r = pd.date_range(start=(pd.Timestamp(date_in[k])
                                    -pd.to_timedelta(str(time_incr_clc)+'s')),
                                periods=tnum+1,
                                freq=str(time_incr_clc)+'s')
        data_clc_in[0,:] =  dat_r# daterange
        # print('Date: ', date_in[k], type(date_in[k]))
        # print('Date-year: ', date_in[k].year)#astype('datetime64[Y]').astype(int) + 1970)#dt.year)
        #TODO: check, if initial dates are consistent

        try:
            # input power value
            P_in    = power_in[k] # in kW ?
            data_clc_in[pow_idx,:] = P_in
            ### further input values
            if dmnd_idx is not None:
                data_clc_in[dmnd_idx,:] = H2dmnd_in[k] # // in kg/h ???
            if c_electr_idx is not None:
                data_clc_in[c_electr_idx,:] = c_electr_in[k] # // in kg/h ???
            if f_emiss_idx is not None:
                data_clc_in[f_emiss_idx,:] = f_emiss_in[k] # // in kg/h ???
            nt_clc_in = NT0(*data_clc_in)
            # print('nt_clc_in: ', nt_clc_in)
            #data_clc_out = linner.subloop(obj, data_clc_in, tnum, time_incr_clc, )
            nt_clc_out = linner.subloop(obj, nt_clc_in, tnum, time_incr_clc, ini=not(bool(k)))

            # extract data to be stored
            #df_out = simu_inst.df[simu_inst.lst_svals] # store only columns in key_lst
            #data_tb_stored = data_clc_out[col_idxs]

        except:
            no_error = False
            #obj.logger.warning(' -!!!-  \n -> an error occurred in loop number: ', k)
            obj.logger.warning(f'Error occurred in mainloop k={k}')
            traceback.print_exc() # print error message
            #data_clc_out = data_clc_in.copy()
            #data_tb_stored = data_clc_out[col_idxs]
            nt_clc_out = nt_clc_in

        'SEE BELOW! -> check error disabled !'
        # no_error        = xlo.check_err(disabled=False) # returns False, if error

        data_clc_in         = arr_zeros.copy()
        #data_in[:,0]    = data_clc_out[:,-1] # new 'initial' values
        '''
        CHECK: Performance of line below !
        '''
        d_out = nt_clc_out._asdict()
        data_clc_in[:,0]    = np.array([v[-1] for k,v, in d_out.items()])
        #if k <3:
            #print('new data_in: ', data_in)
        # store data
        # print('find index:')
        #dt = date_in[k]
        #print('dt: ', dt)
        #print('lst: ')
        if obj.wrt_output_data:
            try:
                nidx = obj.prms['metadata_sig']['years'].index(date_in[k].year)
            except:
                nidx = obj.prms['metadata_sig']['years'].index(str(date_in[k].year))

            idx_chng = 1 if nidx != idx else 0
            idx = nidx
            # print('idx: ', idx)
            # print(f'idx = {idx}, nidx={nidx}, idx_chng={idx_chng}')
            flpth_data_out = obj.lst_pths_out[idx]
            #pd.DataFrame(data=data_tb_stored.T, index=dat_r).to_csv(flpth_data_out, mode='a', header=False)
            slc = 1 #0 if ((k==0) or idx_chng) else 1

            # print('df: ', pd.DataFrame(d_out, index=dat_r)[slc:].head(2))

            df_out = pd.DataFrame(d_out, index=dat_r)
            df_out = df_out.round(rnd_dct)
            df_out[slc:].to_csv(flpth_data_out, mode='a', header=False)
            #df tbs = pd.DataFrame(data=data_tb_stored.T, )
            #df_tbs['DR'] = dat_r
            #df_tbs.to_csv(flpth_data_out, mode='a', header=False)
            #with open(simu_inst.path_data_out)
            if obj.par_thrm_out:
                fx.parameter_log_thrm(obj, flpth_data_out, obj.par_thrm_a, obj.par_thrm_b, idx=nt_clc_out.date[1:])

        elif obj.scn_setup:
            #print('df_fin: ', df_fin.tail(3))
            #n_df = pd.DataFrame(d_out, index=dat_r)[1:]
            #print('n_df: ',  n_df.tail(3))
            df_fin = pd.concat([df_fin,pd.DataFrame(d_out, index=dat_r)[1:]])

        k +=1
    obj.logger.info('End mainloop ... ')
    '''
    if no_error
        enms = 'without'
    else:
        enms = 'with'
    obj.logger.info(f'Main calculation ended {enms} errors')
    '''
    if obj.scn_setup:
        return df_fin
    else:
        return
