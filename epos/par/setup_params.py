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


import aux.handlingfiles as hf
import aux.handlingparams as hp
import aux.handlingdata as hd
import aux.readingfiles as rf
import aux.writingfiles as wf
import aux.faux as fx

class ScenarioSetup():

    def __init__(self, *flnms):
        for i,flnm in enumerate(flnms):
            if (len(flnms) >0) & i==0:
                filename_super_parameters = flnm
            if (len(flnms) >1) & i==1:
                filename_sig_parameters = flnm

        self.tdd            = datetime.datetime.now()
        self.today_ymd      = str(self.tdd.year) +str(self.tdd.strftime("%m")) + str(self.tdd.strftime("%d"))
        self.today_ymdhs    = self.today_ymd + str(self.tdd.hour) + str(self.tdd.minute)

        self.cwd = os.getcwd()

        lgg, logger_nm = fx.ini_logging(name='ScenarioSetup')
        #print('ScenarioSetup -> logger_nm: ', logger_nm)
        logger = lgg.getLogger(logger_nm)
        self.logger = logger

        ### read parameters
        # -> superior Pars

        self.relpth_sup_par = hf.mk_relpath(filename=filename_super_parameters, rel_parents= 'par')
        print('self.relpth_sup_par', self.relpth_sup_par)
        self.sup_par = rf.read_json_file(basename =self.cwd,
                                        rel_pth=self.relpth_sup_par) # return dict

        # -> sig pars
        if not 'filepath_sig_parameters' in locals():
            self.relpth_sig_par = self.sup_par['filepath_sig_parameters']
        else:
            self.relpth_sig_par = filepath_sig_parameters
        self.sig_par = rf.read_json_file(basename=self.cwd,
                                        rel_pth=self.relpth_sig_par)
        # tec-par-path must be specified in full_dict -ini

        print('self.sig_par',self.sig_par)

        # select scenarios
        ### select simulation-scenarios
        # clc versions
        # sig
        # tec
        # ee_power
        # el_power_fraction
        # scl  ---> scaling or absolute power val ?
        self.scen_dict = hp.select_scenarios(self)

        self.metadata_sig_dicts = hd.check_sig_data(self,)

        #self.name = mk_scenario_name()

        #self.pth_scen_files =
        #self.scen_filename = None

        ########## ------- ###################
        #TODO: clc rated power of plant/ stacks
        ## -> in hp.mk_full_scenario_dict() | l.: 66
        # read tec-params
        #self.tec_par = rf.read_json_file(rel_pth=self.metadata_sig_dicts['relpth_tec_parameters'])
        # clc rated power,


        hp.mk_scen_filenames_and_paths(self, sffx='.json')
        self.fin_scen_lst_o_dict = hp.list_all_final_scen_dicts(self)

        hp.store_scenario_files(self)


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
