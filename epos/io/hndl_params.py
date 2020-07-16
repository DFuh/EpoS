

def read_parameters(self):
    '''
    read parameter file ***enable multiple files ???
    and convert to namedtuple

    returns
    ---------
    namedtuple
    '''
    #with open(self.basepath+"/par_v001.json") as f:
    with open(self.basepath + '/' + self.par_flnm) as f:
        par_dict = json.load(f)
        NT = namedtuple('Params',list(par_dict.keys()))
        nt = NT(**par_dict)

    print('Parameters as dict: ', par_dict)
    return nt
