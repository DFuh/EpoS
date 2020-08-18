"""
Created on Mon Jul 13 08:46:47 2020

@author: dafu_res
"""

import os
import numpy as np
import itertools
import json
import glob
from collections import namedtuple
import multiprocessing as mup
from timeit import default_timer as timer
import datetime
import uuid

import pandas as pd

import tools.hndl_dates as hdldt
import tools.hndl_params as hdlpar
import tools.hndl_files as hdlfls
from tools.sig import aux_v01 as sigauxf

#TODO: paths only relative ->> add absolute/ basepath
#TODO: combine hndl_files - methods !
#TODO: change module name!

class EpoS():
    '''
    primary class for simulation ctrl
    - read parameters
    - make signal instances
    - make simulation instances

    - run simulations *** enable multithreading/ multiprocessing
    '''

    def __init__(self, par_flnm, sig_par_flnm):
        self.basepath = os.getcwd()

        self.flnm_par_bsc = par_flnm
        self.flnm_par_sig = sig_par_flnm
        self.tdd = hdldt.todaysdate()
        self.tday = str(self.tdd.year)+str(self.tdd.month)+str(self.tdd.day)

        #read_setup_file()
        self.parameters = hdlpar.read_par_file(self.basepath, self.flnm_par_bsc, 'Params') # json-file with simulation ctrl parameters
        #self.pth_bsc_in = self.parameters.path_data_input
        '''
        if self.parameters.path_data_output:
            self.out_pth = self.parameters.output_path
        else:
            self.out_pth = None
        '''
        # TODO: create separate file for sig paths and names ?
        self.sig_instances = self.ini_sig_instances() # returns dict !
        self.simu_instances = self.ini_simu_instances()

    ###
    def ini_sig_instances(self,):
        '''
        initialize instances of signal input
        - read given paths
        - check, if file or directory
        -- file -> to list || not list, but dict ??
        -- directory --> all files in dir to list *** only one subdirectory considered
        -- non of both: skip

        TDOD:
        -> REDUNDANT code blocks
        '''
        nt_sig = hdlpar.read_par_file(self.basepath, self.flnm_par_sig, 'Sig')

        #sig_inst = []
        #sig_nms = []
        sig_dict = {}
        for sig in nt_sig:

            if os.path.isfile(sig['path']):
                #sig_inst.append( Sig(name=sig.name, file=sig.path) )
                sig_dict[sig['name']]= Sig(name=sig['name'], file=sig['path'], ref_pth=sig['ref_path'])
                #sig_nms.append(sig['name'])

            elif os.path.isfile(self.basepath+sig['path']):
                #sig_inst.append( Sig(name=sig.name, file=sig.path) )
                #print('sig[name]:', sig['name'])
                sig_dict[sig['name']] = Sig(name=sig['name'], file=self.basepath+sig['path'], ref_pth=sig['ref_path'])
                #sig_nms.append(sig['name'])

            elif os.path.isdir(sig['path']):
                pth = sig['path']
                for num,fl in enumerate(glob.glob(pth + '/*.csv')):
                    #sig_inst.append( Sig(name=nm, file=fl) )
                    fnm = os.path.splitext(os.path.split(fl)[1])[0]
                    nm = sig['name'] + '_' + str(num)+ '_' + fnm
                    sig_dict[nm]= Sig(name=nm, file=fl, ref_pth=sig['ref_path'])
                    #sig_nms.append(nm)

            elif os.path.isdir(self.basepath+sig['path']):
                pth = self.basepath+sig['path']
                for num,fl in enumerate(glob.glob(pth + '/*.csv')):
                    #sig_inst.append( Sig(name=nm, file=fl) )
                    fnm = os.path.splitext(os.path.split(fl)[1])[0]
                    nm = sig['name'] + '_' + str(num)+ '_' + fnm
                    sig_dict[nm]= Sig(name=nm, file=fl, ref_pth=sig['ref_path'])
                    #sig_nms.append(nm)
            else:
                print('Not found! ->> Skipped: ', sig['name'], sig['path'])

        return sig_dict#sig_inst


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
        #TODO: consider case of multiple tec-par-files

        s = [self.parameters.tec_lst,
             self.parameters.scl_lst,
             self.parameters.nominal_power_ee,
             self.parameters.nominal_power_el,
             list(self.sig_instances.keys())
             ]
        print(__name__, '---> s: ', s)
        pssbl_sim = []
        instances = []

        print('all possible simualtions:')
        for num,li in enumerate( list( itertools.product(*s))):
            #print(l)
            #l ->
            print('['+str(num)+']---> l: ', li)
            pssbl_sim.append(li)


        npt = input('Please select simulations to run: [input comma-separated integers or all -> enter]')
        try:
            idxlst = []
            for el in npt.split(','):
                idxlst.append(int(el))
        except:
            idxlst = []
        print('idxlst: ', idxlst)
        for num,l in enumerate( pssbl_sim):
            if num not in idxlst:
                print('Skip: ['+str(num)+']---> l: ', l)
            else:
                instances.append( Simu( b_path = self.basepath,
                                        s_pars = self.parameters,
                                        #bsc_in_pth = self.pth_bsc_in,
                                        tec = l[0],
                                           scl = l[1],
                                           nom_pwr_ee = l[2],
                                           nom_pwr_el = l[3],
                                           sig = self.sig_instances[l[4]],
                                           tday = self.tday,
                                           todd = self.tdd,
                                           #out_pth = self.out_pth
                                           ) )
        return instances
        '''
    def ctrl_simu_run(self,):
        ''''''
        control all runs of simulations
        -- initialize multiprocessing
        --- get number of available cores

        --

        Returns
        -------
        None.
        '''
        '''

        noc = mup.cpu_count()-1 # number of cores (keeping 1 core for os-functionality)
        if len(self.simu_inst) <= noc:
            pass # possible to run all simulations in parallel
        else:
            #need to split list of instances according to available cpu
            pass

        self.run_simu()

        return

    def run_simu(self, multiprcss):

        if multiprcss:
            pass
        else:
            for sim in self.simu_instances:
                pass
        return
        '''
################################################################################

class Simu(EpoS):
    '''
    simulation
    '''

    def __init__(self,  b_path = None, s_pars = None, tec=None, sig=None, scl=None, nom_pwr_el=None, nom_pwr_ee=None, tday=None, todd = None):

        self.s_parameters = s_pars
        self.tag = uuid.uuid1() # unique identifier for simulation

        #self.name = None
        self.basepath = b_path # ???
        self.pth_bsc_in = self.s_parameters.bsc_path_data_input #bsc_in_pth
        self.pth_bsc_out = self.s_parameters.bsc_path_data_output
        self.tec = tec # string
        self.sig = sig # instance?
        self.scl = scl

        self.nominal_power_ee = nom_pwr_ee
        self.nominal_power_el = nom_pwr_el

        self.name = self.create_name(tec, nom_pwr_el, scl, sig.name, nom_pwr_ee,  )

        self.log_dict = self.__dict__.copy() #specs of data
        self.cln_log_dict = self.log_dict.copy()
        #del self.cln_log_dict['tec_parameters'] # avoid duplicating tec_parameters
        del self.cln_log_dict['s_parameters']
        self.today = tday # e.g. 20200817
        ### data input
        self.filepath_sig_input = sig.filepath

        ### data output
        self.path_data_out = hdlfls.mk_output_path(self, sig, tday=tday, name=self.name)
        '''
        self.output_path, self.output_filepath = hdlfls.handleFiles.pth_mirror(sig_filepath = self.sig_input_filepath, # filepath of sig
                                                              #self.basepath,
                                                                                bsc_pth = self.bsc_in_pth, # default: data/in
                                                                                bsc_dir=None, # e.g. '/in'
                                                                                now= self.now,
                                                                                fl_prfx='',
                                                                                mk_newdir=True,
                                                                                out_pth = out_pth)
        '''

        hdlfls.mk_output_dir(self.path_data_out,)
        ### read tec - specific parameters
        #print(' ---- >>> Parameters: ', self.s_parameters)
        self.tec_parameters = self.read_tec_params() # returns named_tuple
        self.log_filename = str(self.tag) +'_'+self.today
        hdlfls.ini_logfile(self) # watch out: self -> simu!!
        # TODO: create logfile



        # timerange
        '''
        #TODO: make datetimeobjects
        # datetime.strptime('2019-01-01', '%Y-%m-%d')
        # datetime.strptime('2019-01-01 00:02:02', '%Y-%m-%d %H:%M:%S')
        # USE pandas !!!
        # ->> pd.to_datetime('2019-01-01 00:02:02')
        if len(self.s_parameters.yrs_to_clc)>0:
            self.sim_yrs = self.s_parameters.yrs_to_clc # check, if datasets applicable! (below)
        else:
            if sig.starttime > self.s_parameters.starttime:
                self.starttime = sig.starttime
            else:
                self.starttime = self.s_parameters.starttime
            if sig.stoptime < self.s_parameters.stoptime:
                self.stoptime = sig.stoptime
            else:
                self.stoptime = self.s_parameters.stoptime

            self.sim_yrs = np.arange(self.starttime.year, self.stoptime.year+1)
        '''

    def create_name(self,*args):
        if len(args)>0:
            ostr = ''
            for nm in args:
                if isinstance(nm, float):
                    nm = str(nm).replace('.', '-')
                ostr += str(nm) +'_'
        else:
            ostr = None
        return ostr

    def read_tec_params(self):
        '''
        read electrolyser-technology specific parameters
        '''
        #print(self.basepath + '/par' )#+ self.s_parameters.parameter_files[self.tec+'_par_file')
        #print(self.s_parameters.parameter_files[self.tec+'_par_file'])
        with open(self.basepath + '/par' + self.s_parameters.parameter_files[self.tec+'_par_file'][0])as tp_f:
            par_dict = json.load(tp_f)
            NT = namedtuple('Params',list(par_dict.keys()))
            nt = NT(**par_dict)

        # tec
        # scl
        # type ?
        # vers = _ _ _ _ - _ _ _ _

        return nt

    def run(self):
        print('-'*10)
        self.time_simu_start = timer()
        print('Running Simulation: ', self.tag)
        print('sig:', self.sig.name)
        print('tec: ', self.tec)

        self.time_simu_end = timer()
        self.time_duration_simu = self.time_simu_end - self.time_simu_start
        self.log_simu_time()
        print('-'*10)

        return

    def mk_meda_dict(self):
        x=None
        return # as dict

    def log_simu_time(self):
        with open(self.path_data_out + '/eposLog_' + self.log_filename + '.txt', 'r') as f:
            flines = f.readlines()

        flines[4] = 'elapsed time for simu: \t \t '+ str(self.time_duration_simu)+'\n'
        with open(self.path_data_out + '/eposLog_' + self.log_filename + '.txt', 'w') as f:
            f.writelines(flines)
        return

class Sig():
    '''
    instances for input-signals
    '''

    def __init__(self, name=None, file=None, ref_pth=None):
        self.name = name
        self.filepath = file
        self.ref_pth = ref_pth
        specs, self.df = self.read_sig()
        self.normalized = False #

        '''
        self.nominal_power = specs.nominal_power
        self.starttime = specs.starttime    # compare to df.date !!
        self.stoptime = specs.stoptime      # compare to df.date !!
        self.timerange = None #?
        '''

    def read_sig(self,):


        # decide wether to use date from pars or from df

        ### read specs
        line_specs_end = sigauxf.find_line(path?, 'end ?')
        if line_specs_end:
            specs= 
            skprws=
        else:
            specs = None
            skprws = None

        ### read data
        #data = None
        df = pd.read_csv(self.filepath, skiprows=skprws)
        return specs, data
