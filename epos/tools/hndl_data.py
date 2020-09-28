'''
data handling
- conversion
- checking
- cleaning ?
'''

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
    df_c['tdiff'] = df_c.Date.diff()

    years = df_c.Date.dt.year.drop_duplicates().to_list() # list of years in sig-df

    # get enddate
    # get list of years
    return {'p_startdate': sd, 'p_enddate': ed, 'p_timediff': df_c.tdiff, 'p_years': years}


def select_df_columns(df_in, lst_colnm):
    '''
    rename columns of given dataframe based on given list
    drop columns, if necessary
    '''
    lst_orig_columns = df_in.columns
    while(len(lst_orig_columns) != len(lst_columns):)
        drp = input('Columns in (orig) dataframe: ', lst_orig_columns, '\n Which columns do you want to drop?')
        try:
            df_out = df_in.drop([drp], axis=1)
        except:
            df_out = df_in.drop(lst_orig_columns[drp])
    return


def mark_df_columns():
    '''
    Mark columns with specific data (power, wind, ...)
    '''
    return
