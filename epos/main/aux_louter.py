'''
auxilliary functions for outer loop
'''

def check_err(disabled=False):#val_in, k, err, av):
    ''' check main loop for basic calculation errors
        from mod_101,
        edited: 2020-09-17
    '''
    if not disabled:
        u = val_in[0][-1] # TODO:check 'position'
        T = val_in[0][-1]
        #print('u,i,k:',u,T,k)
        if ((u == 0.0) & (T == 0.0)) or np.isnan(T) :
            #av.err = True
            print('...abort calculation... Iteration No.:',k)
            return False#av.err
    else:
        print('WARNING: ', check_err.__name__, 'disabled')
    return True

def get_col_indexes(df, cols_lst):
    '''
    return list of indexes for respective column names
    '''
    idx_lst = []
    for item in cols_lst:
        idx_lst.append(df.columns.get_loc(item))
    return idx_lst
