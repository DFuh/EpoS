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

    input_df = input_df.set_index('Date')
    # dmnd_df = dmnd_df.set_index('Date')

    ###

    #input_df = input_df.loc[pd.to_datetime(simu_inst.s_parameters.starttime), pd.to_datetime(simu_inst.s_parameters.stoptime)].copy()
    #print(pd.to_datetime(obj.prms['sig_metadata']['start_date']))

    # TODO: update lines below !!! -> use super-parameters
    sd = pd.to_datetime(obj.prms['metadata_sig']['start_date'])
    ed = pd.to_datetime(obj.prms['metadata_sig']['end_date'])
    #input_df = input_df.loc[pd.to_datetime(simu_inst.s_parameters.starttime): pd.to_datetime(simu_inst.s_parameters.stoptime)]

    input_df = input_df.loc[sd:ed]
    # dmnd_df = dmnd_df.loc[sd:ed]

    #input_df['Date'] = pd.to_datetime(input_df.Date)
    #sidx = input_df[input_df.Date==sd].index[0]
    #eidx = input_df[input_df.Date==ed].index[0]
    #print('sidx: ', sidx)
    #print('eidx: ', eidx)
    #input_df = input_df.iloc[sidx:eidx]
    #print('+++ Sliced df:', input_df)

    date_in   = input_df.index#.to_numpy()
    #date_in = input_df.Date
    #print('+++ ', __name__,'Date_in: ', date_in)
    # print(dmnd_df.head(4))
    # print(dmnd_df.columns)
    # print(obj.prms['nm_col_H2dmnd'])
    # print(dmnd_df['dm_H2_dmnd'])
    power_in = input_df[obj.prms['nm_col_sig']].to_numpy() # sig-input df -> to np.array
    pow_idx = obj.df0.columns.get_loc('P_in')-1 # TOD: Hardcoded !!
    if obj.prms['nm_col_H2dmnd'] in input_df.columns:
        H2dmnd_in = input_df[obj.prms['nm_col_H2dmnd']] # sig-input df -> to np.array
        dmnd_idx = obj.df0.columns.get_loc('dm_H2_dmnd')-1 # TODO: Hardcoded !!
    else:
        dmnd_idx=None
    if obj.prms['nm_col_c_electr']+'_x' in input_df.columns:
        c_electr_in = input_df[obj.prms['nm_col_c_electr']+'_x'] # sig-input df -> to np.array
        c_electr_idx = obj.df0.columns.get_loc('c_electr')-1 # TODO: Hardcoded !!
    else:
        c_electr_idx=None
    if obj.prms['nm_col_f_emiss']+'_x' in input_df.columns:
        f_emiss_in = input_df[obj.prms['nm_col_f_emiss']+'_x'] # sig-input df -> to np.array
        f_emiss_idx = obj.df0.columns.get_loc('f_emiss_spc')-1 # TODO: Hardcoded !!
    else:
        f_emiss_idx=None


    # print('pow_idx: ', pow_idx)
    # length of input-power-df
    len_df_pin          = len(power_in)


    # auxilliary variables
    time_incr_clc   = int(obj.prms['time_incr_clc'])
    tnum            = int(int(obj.prms['metadata_sig']['time_incr']) /  time_incr_clc) # number of inner loops
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
    #data_clc_in[0,0] = 000

    ### set initial Values
    #for key, val in obj.prms
    #idx_T = obj.df0.columns.get_loc()
    #data_clc_in[,0] =

    obj.av.stckfctr = obj.pplnt.number_of_stacks_act*obj.pplnt.number_of_cells_in_stack_act
    # call inner loop
    obj.logger.info('Starting mainloop ... ')
    t0 = time.time()
    idx = 0 # Years-index
    k = 0
    while( no_error & (k < len_df_pin)):
        # datetime column
        #data_clc_in[0,:] = pd.date_range(start=date_in[k], end=date_in[k+1], freq=str(time_incr_clc)+'s') # daterange
        dat_r = pd.date_range(start=(pd.Timestamp(date_in[k])
                                    -pd.to_timedelta(str(time_incr_clc)+'s')),
                                periods=tnum+1,
                                freq=str(time_incr_clc)+'s')
        data_clc_in[0,:] =  dat_r# daterange
        print('Date: ', date_in[k], type(date_in[k]))
        print('Date-year: ', date_in[k].year)#astype('datetime64[Y]').astype(int) + 1970)#dt.year)
        #TODO: check, if initial dates are consistent

        try:
            # input power value
            P_in    = power_in[k] # in kW ?
            data_clc_in[pow_idx,:] = P_in
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
        no_error        = xlo.check_err(disabled=True) # returns False, if error
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
        # dt = date_in[k]
        # print('dt: ', dt)
        #print('lst: ')
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
        pd.DataFrame(d_out, index=dat_r)[slc:].to_csv(flpth_data_out, mode='a', header=False)
        #df tbs = pd.DataFrame(data=data_tb_stored.T, )
        #df_tbs['DR'] = dat_r
        #df_tbs.to_csv(flpth_data_out, mode='a', header=False)
        #with open(simu_inst.path_data_out)

        k +=1

    '''
    if no_error
        enms = 'without'
    else:
        enms = 'with'
    obj.logger.info(f'Main calculation ended {enms} errors')
    '''
    return
