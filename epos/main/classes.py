"""
Created on Mon Jul 13 08:46:47 2020

@author: dafu_res
"""

import os
import itertools
import json
import glob
from collections import namedtuple
import multiprocessing as mup
import datetime

import io.hndl_dates as hdldt
import io.hndl_params as hdlpar
import io.hndl_output as hdlout

#TODO: paths only relative ->> add absolute/ basepath

class EpoS():
    '''
    primary class for simulation ctrl
    - read parameters
    - make signal instances
    - make simulation instances

    - run simulations *** enable multithreading/ multiprocessing
    '''

    def __init__(self,):
        self.basepath = os.getcwd()
        self.tdd = hdldt.todaysdate()

        #read_setup_file()
        self.parameters = hdlpar.read_parameters(self) # json-file with simulation ctrl parameters
        self.sig_instances = self.ini_sig_instances()
        #self.simu_instances = self.ini_simu_instances()

    ###
    def ini_sig_instances(self,):
        '''
        initialize instances of signal input
        - read given paths
        - check, if file or directory
        -- file -> to list
        -- directory --> all files in dir to list *** only one subdirectory considered
        -- non of both: skip
        '''
        sig_inst = []
        for pth in self.parameters.sig_pth:
            if os.path.isfile(pth):
                sig_inst.append( Sig(pth) )
            elif os.path.isdir(pth):
                for fl in glob.glob(pth + '/*.csv'):
                    sig_inst.append( Sig(fl) )
            else:
                print('Skipped: ', pth)
        return sig_inst


    ###
    def ini_simu_instances(self, ):
        '''
        make list of instances for simulation(s)

        Parameters
        ----------
         : TYPE
            DESCRIPTION.

        Returns
        -------
        instances : TYPE list
            DESCRIPTION.

        '''

        s = [self.parameters.tec_lst,
             self.parameters.scl_lst,
             self.parameters.nominal_power_ee,
             self.parameters.nominal_power_el,
             ]
        instances = []
        for l in list(itertools.product(*s)):
            #print(l)
            #l ->
            instances.append( SimuInst( tec = l[0],
                                       scl = l[1],
                                       nom_pwr_ee = l[2],
                                       nom_pwr_e = l[2],
                                       ) )
        return instances

    def ctrl_simu_run(self,):
        '''
        control all runs of simulations
        -- initialize multiprocessing
        --- get number of available cores

        --

        Returns
        -------
        None.

        '''

        noc = mup.cpu_count()-1 # number of cores (keeping 1 core for os-functionality)
        if not len(self.simu_inst) > noc:
            pass # possible to run all simulations in parallel
        else:
            #need to split list of instances according to available cpu
            pass

        self.run_simu()

        return

    def run_simu(self,):

        return
