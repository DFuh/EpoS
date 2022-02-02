'''
runfile
'''

import sys
import os
import time
#import logging
import glob

import datetime

from timeit import default_timer as timer
import multiprocessing as mp #import Process, Pool, cpu_count


#from main.classes import EpoS
from epos.main.simulation import ElSim
import epos.auxf.faux as fx
# TODO: comments !
# TODO implement proper logging https://docs.python.org/3/howto/logging.html
# TODO: implement processbar: https://stackoverflow.com/questions/3160699/python-progress-bar#26761413


def run_simu(simu_inst):
    nmstr = simu_inst.name
    #slogger.info('Run simulation: %s', nmstr)
    #print('process-name: ', current_process().name)
    #print('parent process id:', os.getppid())
    #print('process id:', os.getpid())

    simu_inst.setup_sim()
    simu_inst.run()

    print('+'*15)
    print('Simu.cln_log_dict:')
    #print(simu_inst.cln_log_dict)
    print('+'*15)
    #simu_inst.log_simu_time()
    return


def main(pth_in, nms, *argvs, cwd=None):
    print('pth_in: ', pth_in)
    print('nms: ', nms)
    print('argvs: ', argvs)

    if not cwd:
        epos_path = os.path.dirname(__file__)
        cwd = os.getcwd()
    else:
        epos_path = cwd
    print('epos_path: ', epos_path)
    #print('cwd: ', os.getcwd())
    print('__file__', __file__, os.path.dirname(__file__))

    now = datetime.datetime.now()
    slogger, logger_nm = fx.ini_logging(name=str(now)+'slog',pth=epos_path+'/logfiles')
    #print('Base -> logger_nm: ', logger_nm)
    #slogger = lgg.getLogger(logger_nm)
    #logging.basicConfig(filename=str(now)+'xmpl.log',level=logging.DEBUG)
    slogger.info('+++ start EpoS +++' )

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


    slogger.info('Initialize Simu instances')

    #TODO: what, if pth_lst ?
    if not nm_lst: # if no list of files: use all in dir
        lst = glob.glob(epos_path + '/'+ pth + '/Scen*.json')
    else:
        lst = []
        for nm in nm_lst:
            lst.append(os.path.join(pth,nm))

    # make list of simu-instances
    inst_lst = []
    for flnm in lst:
        print('flnm: ', flnm)

        # Ini simulation instance
        inst_lst.append(ElSim(flnm))

    l_sim_inst = len(inst_lst) # Number of Simulation instances
    if not inst_lst:
        slogger.info('No Simulations initialized. Check Scenario-Files...')
    else:
        noc = mp.cpu_count()
        noc_max = 20 #inst_lst[0].no_ac
        # if noc >noc_max:
        slogger.info(f'Available cores: {noc}')
        slogger.info(f'Number of initialized simulations: {l_sim_inst}')
        ret = input('How many cores shall be used? -> ')
        try:
            noc_act = int(ret)
        except:
            if noc_max > noc:
                noc_act = noc-1
            else:
                noc_act = noc_max
        # noc=noc_max

        slogger.info(f'Starting computations on {noc_act} cores')
        ### ini logging
        #logger_nm = fx.ini_logging(self)
        #logger = logging.getLogger(logger_nm)
        #logger.propagate = False

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




    print(' --- end --- ')
    slogger.info(' +++ End EpoS +++')
    end = timer()
    print(f'elapsed time: {end - start}')
    slogger.info(f'elapsed time: {end - start}')
    return

if __name__ == '__main__':
    '''
    now = datetime.datetime.now()
    lgg, logger_nm = fx.ini_logging(name=str(now)+'slog',pth='logfiles')
    print('Base -> logger_nm: ', logger_nm)
    slogger = lgg.getLogger(logger_nm)
    #logging.basicConfig(filename=str(now)+'xmpl.log',level=logging.DEBUG)
    slogger.info('+++ start EpoS +++' )


    print('sig instances: ', Sim.sig_instances)
    print('nm_lst_sig: ', list(Sim.sig_instances.keys()))
    print('simu_inst: ', Sim.simu_instances)
    print('todd: ', Sim.tdd)

    args = sys.argv
    pth='data/scen/dftest/20201117'
    flnm= []#['Scen__PEM_0.6_1_sig_05_WEAoff_2000__.json']
    '''
    args = sys.argv
    if len(args) <2:
        print('...use default path...')
        pth = 'data/scen/test/20211020'#'data/scen/test/20210526' #'data/scen/dftest/20201117'
    if len(args) >1:
        pth = os.path.join('data/scen/', args[1])
        print('pth: ', pth)
    if len(args) > 2:
        nms = args[2]
    else:
        nms = []
    cwd = os.getcwd()
    #print('args: ', args)
    #pth='data/scen/dftest/20201117'
    #flnm= []#['Scen__PEM_0.6_1_sig_05_WEAoff_2000__.json']
    main(pth, nms, cwd=cwd)#, flnm, args)
    '''
    print(' --- end --- ')
    slogger.info(' +++ End EpoS +++')
    '''
