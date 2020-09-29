'''
reading and writing to files
'''
import json
import math
import pandas as pd

def write_to_json(filename, data):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    return


def symline(symbol='-', length=10):
    return symbol*length

def txt_symline(text='', symbol='-', length=20):
    lnm = len(text)
    otxt = '-' * (length - math.ceil(lnm/2)) + ' '+ text +' '+ '-' * (length - math.ceil(lnm/2) -1)
    return otxt

def write_to_csv(filepath, datasets=[], header=None, headline=[], footline=[], data_headline=[], data_footline=[], sep=',', spcsym=' ', l0=20):
    '''
    write data to csv
    filepath, -> full path to file
    datasets=[], -> list of data(containers) to be stored
    header=None, -> data-header
    headline=[], -> headline(s) to be plotted on top of doc
    footline=[], -> footline(s) to be written on bottom of doc
    data_headline=[], -> line(s) to be added before each dataset
    data_footline=[], -> line(s) to be added after each dataset
    sep=',',
    spcsym=' ',
    l0=20
    '''
    '''
    for i,data in enumerate(datasets):
        if isinstance(data, dict): # dict input
    '''
    with open(filepath, 'w') as f:

            lbr = '\n' # linebreak
            for hl in headline:
                f.write(hl+lbr)

            #for dhl, dfl in zip(data_headline, data_footline):
            #    f.write(dhl+lbr)
            for i,data in enumerate(datasets):
                #f.write(data_headline[i]+lbr)
                #print('data: ', dataset)

                f.write(data_headline[i]+lbr)
                if isinstance(data, dict):
                    for key, value in data.items():
                        print(key, value)
                        if not hasattr(value, 'items'):
                            spaces = spcsym * (l0-len(key))
                            f.write(f'\t {key}'+sep + spaces + f'{value}'+lbr)
                        else:
                            f.write(f'{key}'+lbr)
                            for subkey, subvalue in value:
                                spaces = spcsym * l0-len(subkey)
                                f.write(f'\t {subkey}'+sep + spaces + f'{subvalue}'+lbr)

                elif isinstance(data, pd.DataFrame): # df input
                    if header: # use header, if specified, don't use header from df
                        with open(filepath, 'wr') as f:
                            f.write(header+lbr)
                            header=False
                    data.to_csv(f, index=False, sep=sep) #specify header?
                    #for dfl in data_footline:
                if data_footline[i]:
                    f.write(data_footline[i]+lbr)

                if footline:
                    for fl in footline:
                        f.write(fl+lbr)
        #elif isinstance(data, pd.DataFrame): # df input
        #    data.to_csv(filepath)
    return
