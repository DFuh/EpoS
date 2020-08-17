'''
runfile
'''

import sys
import os
import time
from timeit import default_timer as timer
from multiprocessing import Process, Pool, cpu_count

from main.classes import EpoS
# TODO: comments !

def run_simu(simu_inst):
    print('Run simulation: ', simu_inst.tag)
    #print('process-name: ', current_process().name)
    #print('parent process id:', os.getppid())
    #print('process id:', os.getpid())
    simu_inst.run()
    return


def main(Simu):

    start = timer()
    noc = cpu_count()
    print(f'starting computations on {noc} cores')

    values = Simu.simu_instances #(2, 4, 6, 8)

    #with Pool() as pool:
    #    res = pool.map(run_simu, values)
        #print(res)
    process_lst = []

    l_inst = len(Simu.simu_instances)

    for i in range(l_inst):
        inst = Simu.simu_instances[i]
        p = Process(name=inst.tag, target= run_simu, args=(inst,))
        p.start()
        process_lst.append(p)
    for prc in process_lst:
        prc.join()

    end = timer()
    print(f'elapsed time: {end - start}')
    return

if __name__ == '__main__':
    print('+++ start EpoS +++' )
    args = sys.argv
    la = len(args)
    if  la <2:
        # use default parameter-files
        print('_ using default parameter-files')
        flnm_par = 'par_v001.json'
        flnm_sig = 'sig_par_v001.json'
    elif la<3:
        flnm_par = args[1]
        print('_ using default signal parameter-file')
        flnm_sig = 'sig_par_v001.json'
    else:
        # use specified
        flnm_par = args[1]
        flnm_sig = args[2]
    print(f'Parameter files: \n \t{flnm_par}\n \t{flnm_sig}')

    Sim = EpoS(flnm_par, flnm_sig)

    '''
    print('sig instances: ', Sim.sig_instances)
    print('nm_lst_sig: ', list(Sim.sig_instances.keys()))
    print('simu_inst: ', Sim.simu_instances)
    print('todd: ', Sim.tdd)

    main(Sim)
    '''
    print(' --- end --- ')
