'''
functions/methods for reading files
'''

def read_json_file(filename):
    '''
    read json file from full_path
    '''
    with open(filename, 'r') as f:
        file = json.load(f)
    return file




#----------------------------------------------------

def read_in_parameter_file(obj, pth_par_file):
    '''
    read full parameter file for Simulation
    '''
    if not os.path.isfile(pth_par_file):
        flnm = obj.basepath + '/par/' + pth_par_file

    return parameters


def read_in_signal_dataset(self):
    '''
    read signal dataset
    '''
    self.prms.filepath_sig_input

    return
