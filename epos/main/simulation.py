'''
contains basic simualtion class
'''
import sys
import os
import datetime
import pandas as pd
#import logging
import uuid

import time
from collections import namedtuple
from pathlib import Path

import epos.auxf.readingfiles as rf

import epos.auxf.handlingfiles as hf
# import epos.auxf.handlingparams as hp
import epos.auxf.handlingdata as hd
import epos.auxf.faux as fx
from epos.clc import ctrl
from epos.clc import teco_matbal_v20 as tecomb

from epos.main import louter


class ElSim():

    def __init__(self, scenario_filename=None, home_dir=None, curr_dir=None,
                        full_simu=True, scn_setup=False,scn_dct=None, logpath=None):
        ### auxilliary parameters
        #logging.basicConfig(filename='example_df.log',level=logging.INFO)

        # date
        self.tdd            = datetime.datetime.now()
        self.today_ymd      = self.tdd.strftime("%Y%m%d")
        self.today_ymdhs    = self.tdd.strftime("%Y%m%d%H%M")
        #self.cwd = os.getcwd()
        # self.cwd = Path(__file__).parents[1]
        self.dir_home = home_dir
        self.dir_epos = curr_dir

        self.scn_setup = scn_setup
        self.full_simu = full_simu

        # self.hndl_prms = hp
        ### Parameters and Mode
        if self.scn_setup:
            self.par_thrm_out=False
            self.wrt_output_data = False
            self.prms = scn_dct
            logpth = os.path.dirname(scn_dct['scen_filepath']) # Set path of Scenario file as logpath
            # print('logpth: ', logpth)
            self.full_simu=False

            ### name
            self.name = self.prms['scen_name'].replace('Scen','Sim')

            # print('prms: ', self.prms.keys())
            pth_data_bumptest = hf.mk_abspath(self,
                                pth=scn_dct.get("relpth_bumptest_data", None))
            (self.metadata_input,
            self.data_input) = rf.read_in_dataset(self, pth_data_bumptest,
                                                search_key="end Sig - metadata")
            self.data_input['Power'] = self.data_input['Power'] * 60000

            # print('Data-Input: ', self.data_input.head(5))

        elif scenario_filename is not None:
            if not scn_dct:
                self.prms = rf.read_json_file(abspth_to_fl=scenario_filename) # Parameters as dict
            else:
                self.prms = scn_dct
            self.wrt_output_data = True
            self.par_thrm_out = False
            # if self.prms['fctr_scl_sig']: # already implemented in louter !!!
            #    self.data_input['Power'] = self.data_input['Power'] * self.prms['fctr_scl_sig']

            ### name and tag
            self.name = self.prms['scen_name'].replace('Scen','Sim')
            self.tag = uuid.uuid1()
            self.no_ac = self.prms.get('number_of_cores_max',False) # Limit for utilized cores

            ### output path
            try:
                self.pth_data_out = os.path.join(self.prms['pth_data_out'],
                                                    self.prms['reldir_data_output'],
                                                    self.today_ymd,
                                                    self.prms['scen_name'])
                self.flpth_out_basic = os.path.join(self.pth_data_out, str(self.tag))
            except:

                self.pth_data_out = None
                self.flpth_out_basic = None
                print('self.prms[pth_data_out]: ', self.prms['pth_data_out'])
                print('self.prms[reldir_data_output]: ',self.prms['reldir_data_output'])
                print('self.prms[scen_name]:', self.prms['scen_name'])
                print('--- Pth_out in simulatio_init: ', self.pth_data_out)

            ### mk_dir
            if self.pth_data_out is not None:
                if self.full_simu:
                    hf.mk_dir(self.pth_data_out)
                if logpath is not None:
                    logpth = logpath
                else:
                    logpth = self.pth_data_out
            else:
                logpth=os.getcwd() #
        else:
            self.logger.info(' --- No scenario file available --- ')





        ### ini logging
        # logpth = os.path.join(self.cwd, 'logfiles')
        #lgg, logger_nm = fx.ini_logging(self, pth=logpth)
        #print('Simu -> logger_nm: ', logger_nm)
        #logger = lgg.getLogger(logger_nm)
        logger, logger_nm = fx.ini_logging(self, pth=logpth,
                                            notest=self.full_simu)
        self.logger = logger
        # TODO: pth to logfile hardcoded

        if self.full_simu: # No sig needed for testing
            ### Read Power input data
            self.metadata_input, self.data_input = hd.read_data(self)
            # print('--> (data, simu): ', self.data_input.head(4))
            # self.metadata_sig, self.data_sig = rf.read_in_dataset(self,
            #                            rel_flpth=self.prms['relpth_sig_data'],
            #                            search_key=self.prms['searchkey_sig_metadata'])
            ### Read H2-Demand data
            # if self.prms['relpth_H2dmnd_data'] is not None:
            #    self.metadata_H2dmnd, self.data_H2dmnd = rf.read_in_dataset(self,
            #                            rel_flpth=self.prms['relpth_H2dmnd_data'],
            #                            search_key=self.prms['searchkey_H2dmnd_metadata'])

            #print(f'Data head of simulation {self.name}:', self.data_sig.head())
            # check properties of df ?
            #hf.ini_logfile(self,)

            ### handle sig metadata /(MOVE TO ???)
            if self.prms['metadata_sig'] is None:
                self.prms['metadata_sig'] = {}
                self.prms['metadata_sig']['how'] = 'extracted_from_df'
                self.data_input['Date'] = pd.to_datetime(self.data_input['Date'])
                sig_prop = hd.get_properties_df(self.data_input)
                self.prms['metadata_sig']['start_date'] = sig_prop['start_date']
                self.prms['metadata_sig']['end_date'] = sig_prop['end_date']
                self.prms['metadata_sig']['years'] = list(sig_prop['years'])
                self.prms['metadata_sig']['time_incr'] = int(sig_prop['time_incr'])
                self.prms['metadata_sig']['generator_technology'] = '-'

            ### ini output data
            self.df0, self.df0_keys, self.lst_pths_out = hd.ini_data_output(self,)
            self.prms['pathlist_outputfiles'] = self.lst_pths_out
            ### Store actual Params
            hf.store_simu_params(self, )

        #logger.info('Today_ymdhs: %s', self.today_ymdhs)

        #self.calc_modules = fx.ini_clc_versions(self,) # returns tuple
        logger.info('Initialized Simulation: %s \n', self.name)
        #self.parameters_tec = fx.ini_tec_params(self,)

    #### TODO: to be copmnsidered (?)
    # -- power control (PID) ? -> Pin <-> Pact
    # --> grid technology: are very short peaks less harmful to grid stability?

    def setup_sim(self, test=False, testmode=None):
        self.logger.info('Setup parameters')
        self.logger.info('Mode: %s \n', str(testmode))
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
        if not test:
            self.pstrg = hd.dct_to_nt(self.prms['parameters_strg']['strg_0'], subkey='value') # Storage parameters as namedtuple

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
        self.pid_ctrl = ctrl.PID_controller(setpoint=self.pcll.temperature_nominal)
        self.pid_ctrl.reset()
        self.pid_ctrl.tune_man(*self.bop.pid_parameters)

        return

    def run(self,):
        #logging.info
        self.logger.info('Run Simulation: %s \n', self.name)
        t0 = time.time()

        # if False:
            ### EL calculations
        louter.mainloop(self, )
        t1 = time.time()


        ### Process data
        # print('self.lst_pths_out: ',self.lst_pths_out)
        lst_df, lst_meda = hf.final_fl_to_df(self.lst_pths_out)

        ### add data


        ### Run Storage Model (optional)

        if self.prms['storage_clc_iso']:
            self.logger.info('Run Storage-Simu ... ')
            #try:
                # strg_v19
                # fin_df = self.clc_m.strg.clc_strg_state_iso(self, lst_df, )
                # lst_df = hd.slice_df_by_years(fin_df)
                #strg_v20
            lst_df = self.clc_m.strg.clc_strg_state_iso(self, lst_df, )
            #except:
            #    self.logger.info('Running Storage-Simu failed ... ')
        t2 = time.time()
        dt0 = t1-t0 # time in seconds
        dt1 = t2-t0 # time in seconds

        ### Rewrite files

        # TODO: update scenario/ parameters !!!

        if True:
            ### extract data (decrease size of files to use)
            self.logger.info('Extract data ... ')
            # matbal_df_lst, yr_lst, extr_df_lst, extr_meda_lst
            (lst_matbal_dfs,
            lst_yrs,
            lst_extr_df,
            lst_extr_meda,
            lst_extr_pths   )= hd.extract_data(self, lst_df, lst_meda,
                                    pths_orig_data=self.lst_pths_out,
                                    extr_keys=self.prms['keys_data_extraction'],
                                    extr_nms=self.prms['nms_data_extraction'],
                                    basic_pth=self.flpth_out_basic)
            ### Write extracted data to csv
            self.logger.info('Write extracted data ... ')
            if any(isinstance(i, list) for i in lst_extr_pths):
                for k,extr_pths in enumerate(lst_extr_pths):
                    hf.rewrite_output_files(self, extr_pths,
                                            lst_extr_meda[k], lst_extr_df[k])
            else:
                hf.rewrite_output_files(self, lst_extr_pths,
                                        lst_extr_meda, lst_extr_df)
                                        # update meda dict -> insert note on extraction
            # update lst of pth out
            # extract dfs

            # mk matbal file for elTeco
            self.logger.info('Make Matbal df ...')
            tecomb.make_matbal_df(self, lst_matbal_dfs, lst_meda, lst_yrs)

        ### ReWrite Original (full) data to csv
        # TODO: Ensure correct order !

        hf.update_metadata(lst_meda, dt0, dt1)
        self.logger.info('Rewrite output files ... ')
        hf.rewrite_output_files(self, self.lst_pths_out, lst_meda, lst_df)

        #logging.info
        self.logger.info('\n ---- End Simulation: %s', self.name)

        # TODO: update output data

        return None
