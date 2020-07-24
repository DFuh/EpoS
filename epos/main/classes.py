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
import uuid

import pandas as pd

import tools.hndl_dates as hdldt
import tools.hndl_params as hdlpar
import tools.hndl_files as hdlfls

#TODO: paths only relative ->> add absolute/ basepath
#TODO: combine hndl_files - methods !

class EpoS():
    '''
    primary class for simulation ctrl
    - read parameters
    - make signal instances
    - make simulation instances

    - run simulations *** enable multithreading/ multiprocessing
    '''

    def __init__(self, par_flnm, sig_flnm):
        self.basepath = os.getcwd()

        self.par_flnm = par_flnm
        self.sig_flnm = sig_flnm
        self.tdd = hdldt.todaysdate()

        #read_setup_file()
        self.parameters = hdlpar.read_par_file(self.basepath, self.par_flnm, 'Params') # json-file with simulation ctrl parameters
        self.bsc_in_pth = self.parameters.basic_input_path
        if len(self.parameters.output_path) >0:
            self.out_pth = self.parameters.output_path
        else:
            self.out_pth = None

        # TODO: create separate file for sig paths and names ?
        self.sig_instances = self.ini_sig_instances() # returns dict !
        self.simu_instances = self.ini_simu_instances()

    ###
    def ini_sig_instances(self,):
        '''
        initialize instances of signal input
        - read given paths
        - check, if file or directory
        -- file -> to list
        -- directory --> all files in dir to list *** only one subdirectory considered
        -- non of both: skip

        TDOD:
        -> REDUNDANT code blocks
        '''
        nt_sig = hdlpar.read_par_file(self.basepath, self.sig_flnm, 'Sig')

        #sig_inst = []
        #sig_nms = []
        sig_dict = {}
        for sig in nt_sig:

            if os.path.isfile(sig['path']):
                #sig_inst.append( Sig(name=sig.name, file=sig.path) )
                sig_dict[sig['name']]= Sig(name=sig['name'], file=sig['path'])
                #sig_nms.append(sig['name'])

            elif os.path.isfile(self.basepath+sig['path']):
                #sig_inst.append( Sig(name=sig.name, file=sig.path) )
                #print('sig[name]:', sig['name'])
                sig_dict[sig['name']] = Sig(name=sig['name'], file=self.basepath+sig['path'])
                #sig_nms.append(sig['name'])

            elif os.path.isdir(sig['path']):
                pth = sig['path']
                for num,fl in enumerate(glob.glob(pth + '/*.csv')):
                    #sig_inst.append( Sig(name=nm, file=fl) )
                    fnm = os.path.splitext(os.path.split(fl)[1])[0]
                    nm = sig['name'] + '_' + str(num)+ '_' + fnm
                    sig_dict[nm]= Sig(name=nm, file=fl)
                    #sig_nms.append(nm)

            elif os.path.isdir(self.basepath+sig['path']):
                pth = self.basepath+sig['path']
                for num,fl in enumerate(glob.glob(pth + '/*.csv')):
                    #sig_inst.append( Sig(name=nm, file=fl) )
                    fnm = os.path.splitext(os.path.split(fl)[1])[0]
                    nm = sig['name'] + '_' + str(num)+ '_' + fnm
                    sig_dict[nm]= Sig(name=nm, file=fl)
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
                instances.append( Simu( bsc_in_pth = self.bsc_in_pth,
                                        tec = l[0],
                                           scl = l[1],
                                           nom_pwr_ee = l[2],
                                           nom_pwr_el = l[3],
                                           sig = self.sig_instances[l[4]],
                                           now = self.tdd,
                                           out_pth = self.out_pth
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

    def run_simu(self, multiprcss):

        if multiprcss:
            pass
        else:
            for sim in self.simu_instances:
                pass
        return

################################################################################

class Simu():
    '''
    simulation
    '''

    def __init__(self,  bsc_in_pth=None, tec=None, sig=None, scl=None, nom_pwr_el=None, nom_pwr_ee=None, now=None, out_pth = None):

        self.tag = uuid.uuid1() # unique identifier for simulation

        self.name = None
        #self.basepath = basepath
        self.bsc_in_pth = bsc_in_pth
        self.tec = tec # string
        self.sig = sig # instance?
        self.scl = scl
        self.nominal_power_ee = nom_pwr_ee
        self.nominal_power_el = nom_pwr_el

        self.now = now
        ### data input
        self.input_filepath = sig.filepath

        ### data output
        self.output_path, self.output_filepath = hdlfls.handleFiles.pth_mirror(self.input_filepath,
                                                                                #self.basepath,
                                                                                self.bsc_in_pth,
                                                                                now= self.now,
                                                                                fl_prfx='',
                                                                                mk_newdir=True,
                                                                                out_pth = out_pth)

        ### read tec - specific parameters
        self.parameters = None

        # timerange
        #TODO: make datetimeobjects
        # datetime.strptime('2019-01-01', '%Y-%m-%d')
        # datetime.strptime('2019-01-01 00:02:02', '%Y-%m-%d %H:%M:%S')
        # USE pandas !!!
        # ->> pd.to_datetime('2019-01-01 00:02:02')
        if len(self.yrs_to_clc)>0:
            self.sim_yrs = self.yrs_to_clc # check, if datasets applicable! (below)
        if sig.starttime > starttime:
            self.starttime = sig.starttime
        else:
            self.starttime = starttime
        if sig.stoptime < stoptime:
            self.stoptime = sig.stoptime
        else:
            self.stoptime

        self.sim_yrs = np.arange(starttime.year, stoptime.year+1)



    def read_el_params():
        '''
        read electrolyser-technology specific parameters
        '''

        # tec
        # scl
        # type ?
        # vers = _ _ _ _ - _ _ _ _

        return

    def run_simu(self):

        return

class Sig():
    '''
    instnances for input-signals
    '''

    def __init__(self, name=None, file=None):
        self.name = name
        self.filepath = file
        specs, self.df = self.read_sig()
        '''
        self.nominal_power = specs.nominal_power
        self.starttime = specs.starttime    # compare to df.date !!
        self.stoptime = specs.stoptime      # compare to df.date !!
        self.timerange = None #?
        '''

    def read_sig(self,):
        pd.read_csv(self.filepath)
        ### read specs
        specs= None

        ### read data
        data = None
        return specs, data
