'''
runfile
'''

import sys
import os
import time
from timeit import default_timer as timer
import multiprocessing as mp #import Process, Pool, cpu_count

from main.classes import EpoS
# TODO: comments !
# TODO implement proper logging https://docs.python.org/3/howto/logging.html
# TODO: implement processbar: https://stackoverflow.com/questions/3160699/python-progress-bar#26761413


def run_simu(simu_inst):
    print('Run simulation: ', simu_inst.tag)
    #print('process-name: ', current_process().name)
    #print('parent process id:', os.getppid())
    #print('process id:', os.getpid())

    simu_inst.run()

    print('+'*15)
    print('Simu.cln_log_dict:')
    print(simu_inst.cln_log_dict)
    print('+'*15)
    simu_inst.log_simu_time()
    return


def main(Simu):

    start = timer()
    noc = mp.cpu_count()
    print(f'starting computations on {noc-1} cores')

    #values = Simu.simu_instances #(2, 4, 6, 8)

    #with Pool() as pool:
    #    res = pool.map(run_simu, values)
        #print(res)
    #process_lst = []
    arg_lst = list()
    #l_inst = len(Simu.simu_instances)

    #for i in range(l_inst):
    #    inst = Simu.simu_instances[i]
    #for inst in Simu.simu_instances:
#        arg_lst.append(inst)
        #p = Process(name=inst.tag, target= run_simu, args=(inst,))#
        # TODO: insert Process for storage clc here !
        #p.start()
        #process_lst.append(p)
    #for prc in process_lst:
    #    prc.join()
    p = mp.Pool(noc-1)
    res = p.map(run_simu, Simu.simu_instances) #arg_lst)
    p.join()

    end = timer()
    print(f'elapsed time: {end - start}')
    return

if __name__ == '__main__':
    print('+++ start EpoS +++' )
    args = sys.argv
    la = len(args)
    if  la <2: #TODO: put it in a function
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
    '''
    main(Sim)

    print(' --- end --- ')
