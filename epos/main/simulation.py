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
                        full_simu=True, scn_setup=False, scn_dct=None, logpath=None,
                        name=None):
        """
        Main scenario class.

        Parameters
        ----------
        scenario_filename : str
            Path to a json file with all scenario parameters.
            (standard way of creating simu instances)

        home_dir : str
            +++ still in use??? +++

        curr_dir : str
            +++ ? +++

        full_simu : boolean
            (if True), Run full setup.
            Including
                - initialisation of parameter and data storing
                - reading input (sig) datasets
                - storing simulation-parameters

        scn_setup : boolean
            If True, obj. only used for scenario-setup
            Only read bumptest-data (no full sig-datasets)

        scn_dct : dict
            scenario-file as dict (e.g. for direct/external application)

        logpath : str
            path to directory for storing logfiles

        name : str
            +++ may be deleted +++
            Name for ini of logging, if no simu-instance is given



        Attributes
        ----------

        tdd : Datetime Object
            Current date and time
        today_ymd : str
            Year month and day from current datetime
        today_ymdhs :
            Year, month day, hour and minute from current datetime
        dir_home : str
            Path to dir ??? +++ may be deleted?? +++
        dir_epos : str
            P
        scn_setup : Bool
            Select (set True) in order to invoke ScenarioSetup
        full_simu : Bool
            Select (if True) full setup including data reading
        pth_data_out : str
            pth to data-output directory
        name : str
            Identification name (just for external run -> logging-name)
        tag : str
            uuid or string for identifying simulation files
        logger : object +++ ? +++
            Logging-object
        metadata_input : dict
            Metadata of signal dataset, containing basic properties of input-signal
        data_input : pd.DataFrame
            DataFrame containing input (sig) data
        df0 : pd.DataFrame
            DataFrame with basic structure for initializing output df
        df0_keys : iterable
            Column-keys of df0
        lst_pths_out : list
            List with pths to output files

        Methods
        -------
        setup_scenario

        setup_sim

        run


        Returns
        -------
        df : TYPE
            DESCRIPTION.



        """
        '''


    Parameters
    ----------
    dct : TYPE
        DESCRIPTION.
    sort : TYPE, optional
        DESCRIPTION. The default is False.
    val_key : TYPE, optional
        DESCRIPTION. The default is None.



    '''
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
        self.pth_data_out = None

        self.prms = None
        self.wrt_output_data = False
        self.par_thrm_out = False
        self.flpth_out_basic = None

        #self.name = None
        self.tag = None
        #self.logger = None

        self.metadata_input = None
        self.data_input = None

        self.df0 = None
        self.df0_keys = None
        self.lst_pths_out = None

        # self.hndl_prms = hp
        ### Parameters and Mode
        if self.scn_setup:
            self.setup_scenario(scn_dct)
            self.name = 'ScenarioSetup'
            no_scen = False

        elif (scenario_filename is not None) or (scn_dct is not None):
            if not scn_dct:
                self.prms = rf.read_json_file(abspth_to_fl=scenario_filename) # Parameters as dict
            else:
                self.prms = scn_dct
            if self.prms.get('parameters',False):
                self.prms = self.prms['parameters']
            self.wrt_output_data = True
            self.par_thrm_out = False
            # if self.prms['fctr_scl_sig']: # already implemented in louter !!!
            #    self.data_input['Power'] = self.data_input['Power'] * self.prms['fctr_scl_sig']

            ### name and tag
            self.name = self.prms.get('scen_name',name).replace('Scen','Sim')
            self.tag = uuid.uuid1()
            self.no_ac = self.prms.get('number_of_cores_max',False) # Limit for utilized cores

            ### output path
            if self.full_simu:
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
            else:
                self.pth_data_out = None
            ### mk_dir
            no_scen = False
        else:
            self.name = 'ElSim'
            no_scen = True



        if self.pth_data_out is not None:
            if self.full_simu:
                hf.mk_dir(self.pth_data_out)
            #if logpath is not None:
            #    logpth = logpath
            #else:
                logpth = self.pth_data_out
        elif logpath is not None:
            logpth = logpath
        else:
            logpth=os.getcwd() #
        print('logpath: ', logpth)


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
        if no_scen:
            self.logger.info(' --- No scenario file available --- ')
        #self.parameters_tec = fx.ini_tec_params(self,)

    #### TODO: to be copmnsidered (?)
    # -- power control (PID) ? -> Pin <-> Pact
    # --> grid technology: are very short peaks less harmful to grid stability?

    def setup_scenario(self, scn_dct):
        '''
        Set up simulation-scenario
        -> read and prepare bumptest-data

        '''
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
        return

    def setup_sim(self, test=False, testmode=None):
        '''
        Run setup of simulation-instance (scenario)
        - ini logging
        - prepare parameter-sets as named-tuples
        - ini auxilliary values
        - setup PID-ctrl

        Parameters
        ----------
        test : bool
            Just excludes storage-paramteres for testing issues
            The default is False.
        testmode : str
            Just for logging +++ may be deleted +++
            The default is None.

        '''

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
        '''
        Run simulation
        -> invoke mainloop (main simulation)
        -> run storage simu (if selected)
        -> extract data / make materialbalance (for elTeco)
        -> (re)write outputfiles
            - simu data output
            - matbal data
        '''
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
            (lst_matbal_dfs, # Yearly matbal data | list of dfs
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
            tecomb.make_matbal_df(self, lst_matbal_dfs, lst_meda, lst_yrs) # Just data output


            ### Edit 202206 -> MAke high-resolution Matbal
            self.logger.info('Make hires-matbal df ...')
            try:
                tecomb.make_hires_matbal(self, lst_df) # Extraction and output
                self.logger.info('Successfully made hires-matbal df ...')
            except:
                self.logger.info('Applying hires-matbal df failed ...')
        ### ReWrite Original (full) data to csv
        # TODO: Ensure correct order !

        hf.update_metadata(lst_meda, dt0, dt1)
        self.logger.info('Rewrite output files ... ')
        hf.rewrite_output_files(self, self.lst_pths_out, lst_meda, lst_df)

        #logging.info
        self.logger.info('\n ---- End Simulation: %s', self.name)

        # TODO: update output data

        return None
