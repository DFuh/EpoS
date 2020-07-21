'''
runfile
'''
from main.classes import EpoS

if __name__ == '__main__':
    flnm_par = 'par_v001.json'
    flnm_sig = 'sig_par_v001.json'
    Run = EpoS(flnm_par, flnm_sig)
    print('sig instances: ', Run.sig_instances)
    print('nm_lst_sig: ', list(Run.sig_instances.keys()))
    print('simu_inst: ', Run.simu_instances)
    print('todd: ', Run.tdd)
