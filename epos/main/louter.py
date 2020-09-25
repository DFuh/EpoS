'''
outer loop
'''
import numpy as np
import pandas as pd
import traceback
from main import aux_louter as xlo
from main import linner

print(__name__, ' imported.')

def mainloop(simu_inst, ):




    # power input data
    print('+++ Input df:', simu_inst.sig.df)
    print('+++ Input df types:', simu_inst.sig.df.dtypes)
    input_df = simu_inst.sig.df.copy()
    input_df = input_df.set_index('Date')
    #input_df = input_df.loc[pd.to_datetime(simu_inst.s_parameters.starttime), pd.to_datetime(simu_inst.s_parameters.stoptime)].copy()
    print(pd.to_datetime(simu_inst.s_parameters.starttime))
    #input_df = input_df.loc[pd.to_datetime(simu_inst.s_parameters.starttime): pd.to_datetime(simu_inst.s_parameters.stoptime)]
    input_df = input_df.loc[simu_inst.s_parameters.starttime: simu_inst.s_parameters.stoptime]
    print('+++ Sliced df:', input_df)
    date_in   = input_df.index#.to_numpy()
    #TODO: check with specified timerange in par
    print('+++ ', __name__,'Date_in: ', date_in)
    power_in = input_df[simu_inst.sig.p_col] # sig-input df -> to np.array

    # length of input-power-df
    len_df_pin          = len(power_in)


    # auxilliary variables
    time_incr_clc   = int(simu_inst.s_parameters.time_increment_clc)
    tnum            = int(simu_inst.sig.timediff /  time_incr_clc) # number of inner loops
    no_error        = True # initial value
    lst_full_cols   = simu_inst.df0.columns
    arr_zeros       = np.zeros((len(lst_full_cols), tnum)) # default array, shape: length of df0.columns, tnum+1 (initial vals)
    col_idxs        = xlo.get_col_indexes(simu_inst.df0, lst_full_cols) # make list of indexesfor selected values from clc_array
    #TODO: check, if array-length is consistent with df-length

    # initial values
    #arr_data_in = simu_inst.df[0].to_numpy()
    # set new input data
    data_clc_in         = arr_zeros.copy()
    #simu_inst.df0.iloc[-1][0] = pd.to_datetime(simu_inst.df0.iloc[-1][0])
    data_clc_in[1:,0]    = simu_inst.df0.iloc[-1][1:].to_numpy() # np array
    data_clc_in[0,0] = 000


    # call inner loop
    print('__________right before loop')
    k = 0
    while( no_error & (k < len_df_pin)):
        # datetime column
        #data_clc_in[0,:] = pd.date_range(start=date_in[k], end=date_in[k+1], freq=str(time_incr_clc)+'s') # daterange
        dat_r = pd.date_range(start=date_in[k], periods=tnum, freq=str(time_incr_clc)+'s')
        data_clc_in[0,:] =  dat_r# daterange
        print('Date: ', date_in[k])
        print('Date-year: ', date_in[k].year)#astype('datetime64[Y]').astype(int) + 1970)#dt.year)
        #TODO: check, if initial dates are consistent

        try:
            # input power value
            P_in    = power_in[k] # in kW ?

            data_clc_out = linner.subloop(data_clc_in, tnum )

            # extract data to be stored
            #df_out = simu_inst.df[simu_inst.lst_svals] # store only columns in key_lst
            data_tb_stored = data_clc_out[col_idxs]

        except:
            no_error = False
            print(' -!!!-  \n -> an error occurred in loop number: ', k)
            traceback.print_exc() # print error message
            data_clc_out = data_clc_in.copy()
            data_tb_stored = data_clc_out[col_idxs]

        no_error        = xlo.check_err(disabled=True) # returns False, if error
        data_in         = arr_zeros.copy()
        data_in[:,0]    = data_clc_out[:,-1] # new 'initial' values
        # store data
        pd.DataFrame(data=data_tb_stored.T, index=dat_r).to_csv(simu_inst.flpth_data_out, mode='a', header=False)
        #with open(simu_inst.path_data_out)

        k +=1

    return
