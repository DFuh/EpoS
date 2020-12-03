'''
contains basic simualtion class
'''
import sys
import os
import datetime
#import logging
import uuid

import time
from collections import namedtuple

import aux.readingfiles as rf
import aux.handlingfiles as hf
import aux.handlingdata as hd
import aux.faux as fx

from main import louter


class ElSim():

    def __init__(self, scenario_filename, full_simu=True):
        ### auxilliary parameters
        #logging.basicConfig(filename='example_df.log',level=logging.INFO)

        # date
        self.tdd            = datetime.datetime.now()
        self.today_ymd      = str(self.tdd.year) +str(self.tdd.strftime("%m")) +str(self.tdd.strftime("%d"))
        self.today_ymdhs    = self.today_ymd +str(self.tdd.strftime("%H")) +str(self.tdd.strftime("%M"))
        self.cwd = os.getcwd()



        ### Parameters
        self.prms = rf.read_json_file(filename=scenario_filename) # Parameters as dict

        ### name and tag
        self.name = self.prms['scen_name'].replace('Scen','Sim')
        self.tag = uuid.uuid1()

        if full_simu: # TODO: what for???
            ### input data
            self.metadata_sig, self.data_sig = rf.read_in_signal_dataset(self,
                                                                        rel_flpth=self.prms['relpth_sig_data'],
                                                                        search_key=self.prms['searchkey_sig_metadata'])
            #print(f'Data head of simulation {self.name}:', self.data_sig.head())
            # check properties of df ?
            #hf.ini_logfile(self,)

            ### ini output data
            self.df0, self.df0_keys, self.lst_pths_out = hd.ini_data_output(self,)

        ### ini logging
        lgg, logger_nm = fx.ini_logging(self, pth='logfiles')
        print('Simu -> logger_nm: ', logger_nm)
        logger = lgg.getLogger(logger_nm)
        self.logger = logger

        #logger.info('Today_ymdhs: %s', self.today_ymdhs)




        #self.calc_modules = fx.ini_clc_versions(self,) # returns tuple
        logger.info('Initialized Simulation: %s', self.name)
        #self.parameters_tec = fx.ini_tec_params(self,)

    #### TODO: to be copmnsidered (?)
    # -- power control (PID) ? -> Pin <-> Pact
    # --> grid technology: are very short peaks less harmful to grid stability?

    def setup_sim(self, ):
        self.logger.info('Setup parameters')
        ### Setup Params
        self.clc_m = fx.ini_clc_versions(self)
        #par = obj.parameters_tec # namedtuple does not work with Pool.map()
        par = self.prms['parameters_tec_el']

        #TODO: ini el-properties (pwr, etc)
        dct_par_plnt = par['plant']                 # Plant parameters as dict
        dct_par_cell = par['cell']                  # Cell PArameters as dict
        dct_par_op = par['operation']             # Operation parameters as dict
        dct_par_elchem = par['electrochemistry']    # Electrochemistry parameters as dict

        self.pplnt = hd.dct_to_nt(dct_par_plnt, subkey='value') # Plant parameters as namedtuple
        self.pcll = hd.dct_to_nt(dct_par_cell, subkey='value')  # Cell parameters as namedtuple
        self.pop = hd.dct_to_nt(dct_par_op, subkey='value')     # Operation parameters as namedtuple
        self.pec = hd.dct_to_nt(dct_par_elchem, subkey='value') # Electrochemistry parameters as namedtuple

        # option: dacite.from_dict to build dataclasss from dict (faster than namedtuple?)

        hd.ini_auxvals(self,)

        ### Setup auxilliary values (dataclass)

        return

    def run(self,):
        #logging.info
        self.logger.info('Run Simulation: %s', self.name)
        t0 = time.time()

        ### Start calculations
        louter.mainloop(self, )

        t1 = time.time()
        #logging.info
        self.logger.info('End Simulation: %s', self.name)
        dt = t1-t0 # time in seconds
        # TODO: update output data

        return None
