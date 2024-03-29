'''
functions/methods for reading files
'''
import os
import json
import csv
import pandas as pd
from pathlib import Path

def read_json_file(abspth_to_fl=None,pth=None, flnm=None):
    '''
    read json file from full_path
    '''
    if not abspth_to_fl:
        abspth_to_fl = os.path.join(pth,flnm)
    if not os.path.exists(abspth_to_fl):
        print('Invalid filename (readingfile.read_json_file): ', abspth_to_fl)
        file = None
    #if not basename:
        #basename = os.getcwd()
    #    basename = Path(__file__).parents[1]
    else:
        try:
            #if rel_pth:
            #    pth_to_fl = os.path.join(basename, rel_pth)
            #else:
            #    if parent_dir:
            #        pth_to_fl = os.path.join(basename,parent_dir,filename)
            #    else:
            #        pth_to_fl = os.path.join(basename,filename)

            with open(abspth_to_fl, 'r') as f:
                file = json.load(f)
        except Exception as e:
            print('Error in json-file: ', e)
            #print('pth_to_file: ', pth_to_fl)
            #print(read_json_file.__name__,
            #f': Something went wrong! -> \n Basename: {basename}, relative path: {rel_pth}, filename: {filename}')
            raise Exception(f'Could not read jsonfile: {abspth_to_fl}')

            file = None
    return file


def NOread_in_signal_dataset(obj, filename=None):
    '''
    read signal dataset
    '''
    rf.rea

    return



#----------------------------------------------------
def read_in_dataset(obj, abspth_to_fl=None, pth=None, flnm=None,
                        skprws=None, headerrow=[0],
                        search_key='end sig', search_key2='metadata', corrskprws=0, nrws=None):
    if not abspth_to_fl:
        abspth_to_fl = os.path.join(pth,flnm)
    if not os.path.exists(abspth_to_fl):
        print('Invalid filename (readingfile.read_in_dataset): ', abspth_to_fl)
        file = None

    # decide wether to use date from pars or from df

    ### read specs
    if search_key is not None:
        line_specs_end = find_line(abspth_to_fl, search_key, s_key2=search_key2)
        print('line_specs_end: ', line_specs_end)
    else:
        line_specs_end = None
    if skprws is not None:
        specs=None
    elif line_specs_end:
        specs   = read_metadata(abspth_to_fl, line_specs_end) # returns dict

        if corrskprws == 0:
            for nc in range(5):
                out = find_line(abspth_to_fl, '----', num_start=line_specs_end+nc)
                if out is not None:
                    corrskprws = out
            skprws  = corrskprws +1
        else:
            skprws  = line_specs_end + corrskprws # 3 + corrskprws
        print('corr rws: ', corrskprws)
        #print('-->Specs: ', specs)
    else:
        specs   = None
        skprws  = None

    ### read data
    #data = None
    print('flpth, skprws, headr: ', abspth_to_fl, skprws, headerrow)
    df = pd.read_csv(abspth_to_fl, skiprows=skprws, header=headerrow, nrows=nrws)
    if df.empty:
        raise Exception('could not read data')
    else:
        try:
            for colnm in ['Date', 'date']:
                if colnm in df.columns:
                    df['Date'] = pd.to_datetime(df[colnm])
        except:
            print(df.head(3))
            raise Exception("No valid Date-column in df")
    #df = df.set_index('Date') # in louter
    return specs, df#specs, data


def find_line(pth_to_file, search_key, s_key2=None, num_end=50, num_start=0):
    '''
        get line in csv, with search key
        or any specified search text
    '''
    with open(pth_to_file, 'r') as f:
        for num, line in enumerate(f):#,1):
            if num >num_start-1:
                # print(f'[{num}]->| line= {line}')
                if (search_key.lower() in line.lower()):# & (s_key2 in line.lower()):
                    return num
                if num > num_end:
                    return None

def read_metadata(pth, ln):
    '''
    read metadata from csv file
    pth -> full pth to file
    ln -> last line of metadata in csv file
    '''
    d = {} # Ini dict
    with open(pth, 'r') as f:
        rf =csv.reader(f)
        for num, line in enumerate(rf): # Loop through rows/lines of file
            if num < ln:
                if line[0].strip()[0].isalpha(): # True, if the first symbol in line is a character
                    key, value = line[0], line[1] # First element as kex in dict, Second element as value (sep: ',')
                    d[key.strip()] = value.strip() # remove whitepspaces
            else:
                break # Not nice, but avoids full loop through big datasets
    return d

'''
def read_in_parameter_file(obj, pth_par_file):

    #read full parameter file for Simulation

    if not os.path.isfile(pth_par_file):
        flnm = obj.basepath + '/par/' + pth_par_file

    return parameters
'''
