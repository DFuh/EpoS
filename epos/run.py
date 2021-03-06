'''
runfile
'''

import sys
import os
import time
import logging
import glob

from timeit import default_timer as timer
import multiprocessing as mp #import Process, Pool, cpu_count

#from main.classes import EpoS
from main.simulation import ElSim
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
    #print(simu_inst.cln_log_dict)
    print('+'*15)
    #simu_inst.log_simu_time()
    return


def main(pth_in, nms, *argvs):
    print('pth_in: ', pth_in)
    print('nms: ', nms)
    print('argvs: ', argvs)

    la = len(argvs)
    if la == 2:
        pth = argvs[1]
    elif la >1:
        pth = argvs[1]
        flnms = argvs[2:]
    else:
        pth = pth_in
        nm_lst = nms

    start = timer()
    noc = mp.cpu_count()
    logging.info(f'starting computations on {noc-1} cores')

    logging.info('Initialize Simu instances')

    #TODO: what, if pth_lst ?
    if not nm_lst: # if no list of files: use all in dir
        lst = glob.glob(pth + '/Scen*.json')
    else:
        lst = []
        for nm in nm_lst:
            lst.append(os.path.join(pth,nm))

    # make list of simu-instances
    inst_lst = []
    for flnm in lst:
        print('flnm: ', flnm)
        inst_lst.append(ElSim(flnm))

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
    res = p.map(run_simu, inst_lst) #arg_lst)
    #p.join()
    p.close()
    end = timer()
    print(f'elapsed time: {end - start}')
    return

if __name__ == '__main__':
    logging.basicConfig(filename='example_df.log',level=logging.INFO)
    logging.info('+++ start EpoS +++' )

    '''
    print('sig instances: ', Sim.sig_instances)
    print('nm_lst_sig: ', list(Sim.sig_instances.keys()))
    print('simu_inst: ', Sim.simu_instances)
    print('todd: ', Sim.tdd)
    '''
    args = sys.argv
    pth='data/scen/20200930'
    flnm= []#['Scen__PEM_0.6_1_sig_05_WEAoff_2000__.json']

    main(pth, flnm, args)

    print(' --- end --- ')
