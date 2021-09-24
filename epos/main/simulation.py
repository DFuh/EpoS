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
from pathlib import Path

import epos.auxf.readingfiles as rf

import epos.auxf.handlingfiles as hf
import epos.auxf.handlingdata as hd
import epos.auxf.faux as fx
from epos.clc import ctrl

from epos.main import louter


class ElSim():

    def __init__(self, scenario_filename, full_simu=True):
        ### auxilliary parameters
        #logging.basicConfig(filename='example_df.log',level=logging.INFO)

        # date
        self.tdd            = datetime.datetime.now()
        self.today_ymd      = self.tdd.strftime("%Y%m%d")
        self.today_ymdhs    = self.tdd.strftime("%Y%m%d%H%M")
        #self.cwd = os.getcwd()
        self.cwd = Path(__file__).parents[1]


        ### Parameters
        self.prms = rf.read_json_file(filename=scenario_filename) # Parameters as dict

        ### name and tag
        self.name = self.prms['scen_name'].replace('Scen','Sim')
        self.tag = uuid.uuid1()

        ### ini logging
        logpth = os.path.join(self.cwd, 'logfiles')
        #lgg, logger_nm = fx.ini_logging(self, pth=logpth)
        #print('Simu -> logger_nm: ', logger_nm)
        #logger = lgg.getLogger(logger_nm)
        logger, logger_nm = fx.ini_logging(self, pth=logpth, notest=full_simu)
        self.logger = logger
        # TODO: pth to logfile hardcoded

        if full_simu: # No sig needed for testing
            ### Read Power input data
            self.metadata_sig, self.data_sig = rf.read_in_dataset(self,
                                        rel_flpth=self.prms['relpth_sig_data'],
                                        search_key=self.prms['searchkey_sig_metadata'])
            ### Read H2-Demand data
            self.metadata_H2dmnd, self.data_H2dmnd = rf.read_in_dataset(self,
                                        rel_flpth=self.prms['relpth_H2dmnd_data'],
                                        search_key=self.prms['searchkey_H2dmnd_metadata'])
            #print(f'Data head of simulation {self.name}:', self.data_sig.head())
            # check properties of df ?
            #hf.ini_logfile(self,)

            ### ini output data
            self.df0, self.df0_keys, self.lst_pths_out = hd.ini_data_output(self,)
            self.prms['pathlist_outputfiles'] = self.lst_pths_out
            ### Store actual Params
            hf.store_simu_params(self, )

        #logger.info('Today_ymdhs: %s', self.today_ymdhs)

        #self.calc_modules = fx.ini_clc_versions(self,) # returns tuple
        logger.info('Initialized Simulation: %s', self.name)
        #self.parameters_tec = fx.ini_tec_params(self,)

    #### TODO: to be copmnsidered (?)
    # -- power control (PID) ? -> Pin <-> Pact
    # --> grid technology: are very short peaks less harmful to grid stability?

    def setup_sim(self, test=False, testmode=None):
        self.logger.info('Setup parameters')
        ### Setup Params
        self.clc_m = fx.ini_clc_versions(self)
        #par = obj.parameters_tec # namedtuple does not work with Pool.map()
        par = self.prms['parameters_tec_el']

        #TODO: ini el-properties (pwr, etc)
        #dct_par_plnt = par['plant']                 # Plant parameters as dict
        #dct_par_cell = par['cell']                  # Cell PArameters as dict
        #dct_par_op = par['operation']             # Operation parameters as dict
        #dct_par_elchem = par['electrochemistry']    # Electrochemistry parameters as dict

        self.pplnt = hd.dct_to_nt(par['plant'], subkey='value') # Plant parameters as namedtuple
        self.pcll = hd.dct_to_nt(par['cell'], subkey='value')  # Cell parameters as namedtuple
        self.pop = hd.dct_to_nt(par['operation'], subkey='value')     # Operation parameters as namedtuple
        self.pec = hd.dct_to_nt(par['electrochemistry'], subkey='value') # Electrochemistry parameters as namedtuple
        self.bop = hd.dct_to_nt(par['periphery'], subkey='value') # Periphery parameters as namedtuple
        self.p = hd.dct_to_nt(par['operation']['nominal_electrode_pressure'],
                                    subkey='value')
        ''' CHECK below !!!'''
        if par['functions']:
            self.fnct = hd.dct_to_nt(par['functions'],
                                        subkey='value')

        #print('self.pop: ', self.pop)
        #print('self.pcll:', self.pcll)
        # option: dacite.from_dict to build dataclasss from dict (faster than namedtuple?)

        # in test-mode, read in refvals
        if test:
            self.refvals = hd.setup_refvals_nt(par['refvals'], testmode)
            print('self.refvals: ', self.refvals)

        ### Setup auxilliary values (dataclass)
        hd.ini_auxvals(self, par)

        ### ??? Setup PID cntrl ->> ??? here ???
        self.pid_ctrl = ctrl.PID_controller()
        self.pid_ctrl.reset()
        self.pid_ctrl.tune_man(*self.bop.pid_parameters)

        return

    def run(self,):
        #logging.info
        self.logger.info('Run Simulation: %s', self.name)
        t0 = time.time()

        ### EL calculations
        louter.mainloop(self, )
        t1 = time.time()

        ### Process data
        df_lst, lst_meda = hf.final_fl_to_df(self.lst_pths_out)

        ### Run Storage Model (optional)
        if self.prms['storage_clc_iso']:
            fin_df = self.clc_m.strg.clc_strg_state_iso(self, df_lst, )
            df_lst = hd.slice_df_by_years(fin_df)
        t2 = time.time()
        dt0 = t1-t0 # time in seconds
        dt1 = t2-t0 # time in seconds

        ### Rewrite files
        # TODO: Ensure correct order !
        hf.update_metadata(lst_meda, dt0, dt1)
        hf.rewrite_output_files(self, self.lst_pths_out, lst_meda, df_lst)

        #logging.info
        self.logger.info('End Simulation: %s', self.name)

        # TODO: update output data

        return None
