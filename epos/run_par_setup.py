'''
run poar setup in order to create scenario files
'''
import os
import sys
from epos.par import setup_params as sepa


def run_scen_setup(flnm_spp=None, flnm_sig_par=None, dir_super_pars=None):
    home_pth = os.path.expanduser('~')
    curr_dir = os.path.dirname(__file__)

    # if flnm_spp is None:
    #    flnm_spp = input('Please insert valid filename of parameters')


    sepa.ScenarioSetup(flnm_spp, flnm_sig_par, home_pth, curr_dir,
                        dir_super_pars)
    return

if __name__ =='__main__':
    fl_nm = sys.argv[1] if len(sys.argv) >1 else None
    add_inp = sys.argv[2] if len(sys.argv) >2 else None

    if not fl_nm:
        inp = input('Please insert valid filename of parameters: /n (Enter or >skip< for default)')
        if (inp == '') or (inp.lower()=='skip'):
            print('Using default super_parameters: par_v01.json')
            fl_nm = 'par_v01.json' # From config or similar file?

    if '.json' not in fl_nm:
        fl_nm = fl_nm+'.json'

    if os.path.isfile(fl_nm):
        print('Using super_parameters: ', fl_nm)
        run_scen_setup(filename_super_parameters, add_inp)
    else:
        print('no valid filename --> skip')
