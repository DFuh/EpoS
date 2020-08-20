import json
from collections import namedtuple

def read_par_file(basepath, flnm, nm, to_namedtuple=False):
    '''
    read parameter file ***enable multiple files ???
    and convert to namedtuple

    returns
    ---------
    namedtuple
    '''
    #with open(self.basepath+"/par_v001.json") as f:
    with open(basepath + '/par/' + flnm) as f:
        par_dict = json.load(f)

    if to_namedtuple:
        NT = namedtuple(nm,list(par_dict.keys()))
        nt = NT(**par_dict)
        return nt
    else:
        return par_dict
