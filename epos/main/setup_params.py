'''
new structural approach:
- setup all parameters and dump into file

'''

'''
what about: ....

- https://pypi.org/project/dataclasses-json/ # not necessary, when using dataclass.asdict ?

'''

# TODO: setup for all PIDs
#TODO: setup for cell-area and stack size

import os
import datetime
from pathlib import Path

import epos.auxf.handlingfiles as hf
import epos.auxf.handlingparams as hp
import epos.auxf.handlingdata as hd
import epos.auxf.readingfiles as rf
import epos.auxf.writingfiles as wf
import epos.auxf.faux as fx

class ScenarioSetup():

    '''
    Class for Scenario Setup
    creates parameter-sets and stores respective files for each simu-scenario

    Attributes
    ----------
    
    '''

    def __init__(self, flnm_supprs, flnm_sigprs,
                        home_dir, curr_dir, slct_dir_supprs=None):

        self.tdd            = datetime.datetime.now()
        self.today_ymd      = str(self.tdd.year) +str(self.tdd.strftime("%m")) + str(self.tdd.strftime("%d"))
        self.today_ymdhs    = self.today_ymd + str(self.tdd.hour) + str(self.tdd.minute)

        # self.cwd = Path(__file__).parents[1] #os.getcwd()
        '''
        TODO : make paths flexible (for storing scenario files)
        '''
        self.dir_home = home_dir
        self.dir_epos = curr_dir
        print('self.dir_home: ', self.dir_home)
        print('self.dir_epos: ', self.dir_epos)
        ### read parameters
        # -> superior Pars
        if not slct_dir_supprs:
            self.abspth_par_sup = os.path.join(self.dir_epos,'par',flnm_supprs)
        else:
            self.abspth_par_sup = os.path.join(slct_dir_supprs,flnm_supprs)
        # print('abspath suppar: ', self.abspth_par_sup)

        ### read super-parameters
        self.par_sup = rf.read_json_file(self.abspth_par_sup) # return dict

        ### Create basic paths
        self.dct_pths = hf.set_pths(self)
        print('dct_pths: ', self.dct_pths)
        self.pth_scen_out = hf.mk_abspath(self, abspth=self.pth_scen_out,
                                            add_pth=self.today_ymd)
        print('pth_scen_out: ', self.pth_scen_out)
        if not os.path.exists(self.pth_scen_out):
            os.makedirs(self.pth_scen_out)

        # ini logging
        self.logger, logger_nm = fx.ini_logging(self, name=self.today_ymdhs+'ScenarioSetup',
                                                pth=self.pth_scen_out)
        #self.logger = logger
        self.logger.info('Ini Setup Log')

        ### Read sig-parameters
        self.logger.info('Read sig parameters ...')
        self.par_sig = rf.read_json_file(pth=self.pth_prms_sig,
                                        flnm=self.par_sup.get('filename_sig_par',None))

        if (self.par_sup is not None) and (self.par_sig is not None):
            # tec-par-path must be specified in full_dict -ini

            # print('self.sig_par',self.sig_par)

            # select scenarios
            ### select simulation-scenarios
            # clc versions
            # sig
            # tec
            # ee_power
            # el_power_fraction
            # scl  ---> scaling or absolute power val ?
            # print('pth_scen_out (0): ',self.pth_scen_out)

            #self.pth_data_out = hf.mk_abspath(self, abspth=self.pth_data_out)
            # print('pth_scen_out: ',self.pth_scen_out)

            self.logger.info('Select Scenarios')
            self.scen_dict = hp.select_scenarios(self)
            # print('Scen Dict: \n', self.scen_dict)

            self.metadata_sig_dicts = hd.check_sig_data(self,)

            hp.mk_scen_filenames_and_paths(self, sffx='.json') # Updates Scen_dict



            self.fin_scen_lst_o_dict = hp.list_all_final_scen_dicts(self)

            hp.store_scenario_files(self)
        else:
            self.logger.info('No valid parameters -> abort')
            self.logger.info('Basic Parameters: %s', str(self.par_sup))
            self.logger.info('Sig Parameters: %s', str(self.par_sig))
        self.logger.info(' -- Done --')
        # -> set all paths

        # TODO: check properties of sig-df ?



        # -> tec pars
        # -> par_data_output

        ### units ?

        ### scenario_name

    def make_scenario_files():
        '''
        extracting necessary parameters for simulation scenarios to store them in separate files
        '''


        return
