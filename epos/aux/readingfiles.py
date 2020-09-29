'''
functions/methods for reading files
'''
import os
import json
import csv
import pandas as pd

def read_json_file(rel_pth=None, basename=None, filename=None, parent_dir=None):
    '''
    read json file from full_path
    '''
    if not basename:
        basename = os.getcwd()
    try:
        if rel_pth:
            pth_to_fl = os.path.join(basename, rel_pth)
        else:
            if parent_dir:
                pth_to_fl = os.path.join(basename,parent_dir,filename)
            else:
                pth_to_fl = os.path.join(basename,filename)

        with open(pth_to_fl, 'r') as f:
            file = json.load(f)
    except:
        print(read_json_file.__name__,f': Something went wrong! -> \n Basename: {basename}, relative path: {rel_pth}, filename: {filename}')
        file = None
    return file


def read_in_signal_dataset(obj, filename=None):
    '''
    read signal dataset
    '''
    rf.rea

    return



#----------------------------------------------------
def read_in_signal_dataset(obj, basename=None, filename=None, search_key='end sig', search_key2='metadata'):
    # decide wether to use date from pars or from df

    ### read specs
    #print(' +++SIG path: ', self.filepath)
    if not basename:
        basename = obj.cwd
    filepath = os.path.join(basename, filename)
    print('Filename for find_line: ', filepath)
    line_specs_end = find_line(filepath, search_key, s_key2=search_key2)
    if line_specs_end:
        specs   = read_metadata(filepath, line_specs_end) # returns dict
        skprws  = line_specs_end + 3
        print('-->Specs: ', specs)
    else:
        specs   = None
        skprws  = None

    ### read data
    #data = None
    df = pd.read_csv(filepath, skiprows=skprws, header=[0])
    df['Date'] = pd.to_datetime(df['Date'])
    #df = df.set_index('Date') # in louter
    return specs, df#specs, data


def find_line(pth_to_file, search_key, s_key2=None, num_end=30):
    '''
        get line in csv, with search key
        or any specified search text
    '''
    with open(pth_to_file, 'r') as f:
        for num, line in enumerate(f):#,1):
            if (search_key in line.lower()) & (s_key2 in line.lower()):
                return num
            if num > num_end:
                return None

def read_metadata(pth, ln):
    d = {}
    with open(pth, 'r') as f:
        rf =csv.reader(f)
        for num, line in enumerate(rf):
            print('num: ', num, 'ln: ', ln)
            if num < ln:
                if line[0].strip()[0].isalpha():
                    print('--line: ', line)
                    key, value = line[0], line[1]
                    d[key.strip()] = value.strip()
    return d

'''
def read_in_parameter_file(obj, pth_par_file):

    #read full parameter file for Simulation

    if not os.path.isfile(pth_par_file):
        flnm = obj.basepath + '/par/' + pth_par_file

    return parameters
'''
