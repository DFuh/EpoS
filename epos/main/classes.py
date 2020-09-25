"""
Created on Mon Jul 13 08:46:47 2020

@author: dafu_res
"""

import sys
import os
import numpy as np
import itertools
import json
import glob
from collections import namedtuple
from importlib import import_module as impm

import multiprocessing as mup
from timeit import default_timer as timer
import datetime
import uuid

import pandas as pd

import tools.hndl_dates as hdldt
import tools.hndl_params as hdlpar
import tools.hndl_files as hdlfls
from tools.sig import aux_v01 as sigauxf

from main import louter

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
        print('tdd:',self.tdd)
        self.tday = str(self.tdd.year)+str(self.tdd.month)+str(self.tdd.day)

        #read_setup_file()
        self.parameters = hdlpar.read_par_file(self.basepath, self.flnm_par_bsc, 'Params', to_namedtuple=True) # json-file with simulation ctrl parameters
        #self.pth_bsc_in = self.parameters.path_data_input
        '''
        if self.parameters.path_data_output:
            self.out_pth = self.parameters.output_path
        else:
            self.out_pth = None
        '''
        print('par, clc_vers.: ', self.parameters.version_clc_files_tec)
        print('par, clc_vers.plr: ', self.parameters.version_clc_files_tec['plr'])
        # TODO: create separate file for sig paths and names ?
        #self.select_simus() !!! #

        self.sig_dict = self.ini_sig_dict() # returns dict !
        self.simu_instances = self.ini_simu_instances()




    ###
    def select_simus():
        '''
        select simulations to run
        just read paths and Params
        -> select simulations
        -> afterwards: initialize simu and sig instances
        '''

        return


    def ini_sig_dict(self,):
        '''
        initialize instances of signal input
        - read given paths
        - check, if file or directorykeys=
        -- file -> to list || not list, but dict ??
        -- directory --> all files in dir to list *** only one subdirectory considered
        -- non of both: skip

        TDOD:
        -> REDUNDANT code blocks
        '''
        nt_sig = hdlpar.read_par_file(self.basepath, self.flnm_par_sig, 'Sig', to_namedtuple=True)

        #sig_inst = []
        #sig_nms = []
        sig_dict = {}
        for sig in nt_sig:
            nosig = False
            if os.path.isfile(sig['path']):
                #sig_inst.append( Sig(name=sig.name, file=sig.path) )
                #sig_dict[sig['name']]= Sig(name=sig['name'], file=sig['path'], ref_pth=sig['ref_path'])
                name = sig['name']
                file = sig['path']
                #ref_pth = sig['ref_pth']
                #sig_nms.append(sig['name'])

            elif os.path.isfile(self.basepath+sig['path']):
                #sig_inst.append( Sig(name=sig.name, file=sig.path) )
                #print('sig[name]:', sig['name'])
                name = sig['name']
                file = self.basepath + sig['path']
                #ref_pth = sig['ref_pth']
                #sig_dict[sig['name']] = Sig(name=sig['name'], file=self.basepath+sig['path'], ref_pth=sig['ref_path'])
                #sig_nms.append(sig['name'])

            elif os.path.isdir(sig['path']):
                pth = sig['path']
                for num,fl in enumerate(glob.glob(pth + '/*.csv')):
                    #sig_inst.append( Sig(name=nm, file=fl) )
                    fnm = os.path.splitext(os.path.split(fl)[1])[0]
                    name = sig['name'] + '_' + str(num)+ '_' + fnm
                    file = fl
                    #ref_pth = sig['ref_pth']
                    #sig_dict[nm]= Sig(name=nm, file=fl, ref_pth=sig['ref_path'])
                    #sig_nms.append(nm)keys=

            elif os.path.isdir(self.basepath+sig['path']):
                pth = self.basepath+sig['path']
                for num,fl in enumerate(glob.glob(pth + '/*.csv')):
                    #sig_inst.append( Sig(name=nm, file=fl) )
                    fnm = os.path.splitext(os.path.split(fl)[1])[0]
                    name = sig['name'] + '_' + str(num)+ '_' + fnm
                    file = fl
                    #ref_pth =
                    #sig_dict[nm]= Sig(name=nm, file=fl, ref_pth=sig['ref_path'])
                    #sig_nms.append(nm)
            else:
                print('Not found! ->> Skipped: ', sig['name'], sig['path'])
                nosig = True

            if not nosig:
                ref_pth = sig['ref_path']
                col_nm_p = sig['clmn_nm_p']
                sig_dict[name] = {'name': name, 'filepath': file, 'ref_pth': ref_pth, 'clmn_nm_p': col_nm_p}
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
             list(self.sig_dict.keys()),
             [None] #space holder for clc_ver
             #self.parameters.clc_versions)
             ]
        print(__name__, '---> s: ', s)
        pssbl_sim = []
        instances = []
        #PossSim = namedtuple('PossSim', ['tec', 'scl', 'rpow_ee', 'rpow_el', 'sig'])
        #s_nt = PossSim(tec=self.parameters.tec_lst)
        #ClcV = namedtuple('ClcV', 'nm plr flws dgr pwr thrm')
        #clc_vers =  # namedtuple
        k_lst = ['tec', 'scl', 'rpow_ee', 'rpow_el', 'sig', 'clc_ver']
        v_klst = ['plr', 'flws', 'dgr' ,'pwr','thrm']
        #sim_ = dict.fromkeys(k_lst)# dict
        #sim_['clc_ver'] = dict.fromkeys(self.parameters.version_clc_files_tec.keys())
        sim_ = {}
        iter_lst = list( itertools.product(*s))

        print('all possible simualtions:')
        for num, i_lst in enumerate( iter_lst):

            #print(l)
            ## add possible clc version:
            #k = 0
            #for clc_mod,val in sim_['clc_ver'].items():
            v_lst = []
            for key in v_klst:
                vers = self.parameters.version_clc_files_tec[key][i_lst[0]]
                v_lst.append(vers)
                #if len(vers) >1:
            #        #for ele in vers:
                        #v_lst.append(ele)
                #else:
                #    v_lst.append(vers)
            #[ele for ele in self.parameters.version_clc_files_tec[key][i_lst[0]] for key in v_klst]
            vi_lst = []
            for ele in list(itertools.product(*v_lst)):
                if ele not in vi_lst:
                    vi_lst.append(ele)

            for k,v_i in enumerate(vi_lst):
                nm = 'simu_'+str(num)+str(k)
                sim_[nm] = dict.fromkeys(k_lst)
                for key, val in zip(k_lst, i_lst[:-1]):
                    sim_[nm][key] = val
                sim_[nm]['clc_ver'] = dict.fromkeys(v_klst)
                for n,key in enumerate(sim_[nm]['clc_ver'].keys()):
                    sim_[nm]['clc_ver'][key] = v_i[n]
            '''
            print(l0)
            print(z)
            for k in range(len(v_lst)):
                nm = 'simu_'+str(num)+str(k)
                sim_[nm] = dict.fromkeys(k_lst)
                for key, val in zip(k_lst, i_lst[:-1]):
                    sim_[nm][key] = val
                sim_[nm]['clc_ver'] = dict.fromkeys(v_lst)
                vers_lst = []
                for key in sim_[nm]['clc_ver'].keys():
                    vers = self.parameters.version_clc_files_tec[key][i_lst[0]]
                    vers_lst.append(vers)
                vers_iter_lst = itertools.product(*vers_lst)
                for  elem in vers_iter_lst:
                    sim_[nm]['clc_ver'][clc_mod] = vers_i #self.parameters.version_clc_files_tec[clc_mod][li[0]]
                    #else:
                    #    sim_[nm]['clc_ver'][key] = vers


                #dgr_ver = self.parameters.version_calculation_files['dgr'][li[0]]
                #plr_ver = self.parameters.version_calculation_files['dgr'][li[0]]
                #k +=1
            '''
            #l ->
        num = 0
        for key, val in sim_.items():
            print(f'[{num}]', key, val)
            num +=1
        #    print('['+str(num)+']---> l: ', li)
        #    pssbl_sim.append(li)


        npt = input('Please select simulations to run: [input comma-separated integers or all -> enter]')
        try:
            idxlst = []
            for el in npt.split(','):
                idxlst.append(int(el))
        except:
            idxlst = []
        print('idxlst: ', idxlst)
        for num,pss in enumerate( sim_):
            if num not in idxlst:
                print('Skip: ['+str(num)+']---> l: ', pss)
        print(sim_)
        for num,ps in enumerate(sim_):
            if num in idxlst:
                print(sim_[ps])#, ps['tec'])
                instances.append( Simu( b_path = self.basepath,
                                        s_pars = self.parameters,
                                        #TODO: delete clc_versions from par ???
                                        v_pars = sim_[ps]['clc_ver'],
                                        #bsc_in_pth = self.pth_bsc_in,
                                        tec = sim_[ps]['tec'],
                                           scl = sim_[ps]['scl'],
                                           nom_pwr_ee = sim_[ps]['rpow_ee'],
                                           nom_pwr_el = sim_[ps]['rpow_ee'],
                                           sig_dct = self.sig_dict[sim_[ps]['sig']],
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

class Simu():#EpoS): #?? inheritance necessary?
    '''
    simulation
    '''

    def __init__(self,  b_path = None, s_pars = None, v_pars=None,
                tec=None, sig_dct=None, scl=None, nom_pwr_el=None, nom_pwr_ee=None, tday=None, todd = None):
        print('sys.path: ', sys.path)
        self.s_parameters   = s_pars
        self.v_pars         = v_pars
        self.basepath       = b_path # ???
        print('v_pars: ', v_pars)
        ### module import
        #nm = '.plr_'+v_pars['plr']#+'.py'
        plr = impm('.plr_'+v_pars['plr'], package='clc')

        self.tag = uuid.uuid1() # unique identifier for simulation

        #self.name = None

        self.pth_bsc_in = self.s_parameters.bsc_path_data_input #bsc_in_pth
        self.pth_bsc_out = self.s_parameters.bsc_path_data_output
        self.tec = tec # string
        self.sig_dct = sig_dct # dict
        self.scl = scl

        ### data input
        self.sig = self.ini_sig_instance()
        self.filepath_sig_input = self.sig.filepath
        #self.nm_sig_pcol = self.parameters.clmn_nm_P # column name of signal df for input power


        ############ years
        # list of years in self.sig.years

        #self.years = self.sig.years #[]

        self.nominal_power_ee = nom_pwr_ee
        self.nominal_power_el = nom_pwr_el

        self.name = self.create_name(tec, nom_pwr_el, scl, self.sig.name, nom_pwr_ee,  )

        self.log_dict = self.__dict__.copy() #specs of data
        self.cln_log_dict = self.log_dict.copy()
        #del self.cln_log_dict['tec_parameters'] # avoid duplicating tec_parameters
        del self.cln_log_dict['s_parameters']
        print('self.cln_log_dict: ', self.cln_log_dict)
        #del self.cln_log_dict['version_clc_files_tec']
        #self.cln_log_dict['clc_versions'] = v_pars
        self.filename_log = None
        self.today = tday # e.g. 20200817
        self.tdd = todd


        ### data output
        self.path_data_out = hdlfls.mk_output_path(self, self.sig, tday=tday, name=self.name)
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
        self.parameters_output = hdlpar.read_par_file(self.basepath,
                                                        self.s_parameters.output_parameters[self.tec+'_output_par_file'][0],
                                                         'doutParams',
                                                         to_namedtuple=False)
        ### read tec - specific parameters
        #print(' ---- >>> Parameters: ', self.s_parameters)
        self.tec_parameters = hdlpar.read_par_file(self.basepath,
                                                    self.s_parameters.parameter_files[self.tec+'_par_file'][0],
                                                    'tecParams',
                                                    to_namedtuple=True) #self.read_tec_params() # returns named_tuple

        hdlfls.ini_logfile(self) # watch out: self -> simu!!
        #hdlfls.ini_df_data_output(self)
        self.df0, self.lst_svals = hdlfls.ini_df_data_output(self)  # returns full df (all variables), output df (selected variables to be stroed), key_lst of values to be stored
        #self.df_out = df0.copy()                    # df for storing results

        hdlfls.ini_output_files(self)
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

    def ini_sig_instance(self):
        sig_inst = Sig(name=self.sig_dct['name'], file = self.sig_dct['filepath'],
                        ref_pth=self.sig_dct['ref_pth'], p_col=self.sig_dct['clmn_nm_p'])
        return sig_inst

    def create_name(self,*args):
        if len(args)>0:
            ostr = ''
            for nm in args:
                if isinstance(nm, float): #replace decimal number for filename
                    nm = str(nm).replace('.', '-')
                ostr += str(nm) +'_'
        else:
            ostr = None
        return ostr


    def read_params(self):
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
        print('|---> Running Simulation: ', self.tag)
        print('sig:', self.sig.name)
        print('tec: ', self.tec)

        louter.mainloop(self)

        self.time_simu_end = timer()
        self.time_duration_simu = self.time_simu_end - self.time_simu_start
        self.log_simu_time()
        print('-'*10)

        return

    def mk_meda_dict(self):
        x=None
        return # as dict

    def log_simu_time(self):
        with open(self.path_data_out + '/'+ self.filename_log + '.txt', 'r') as f:
            flines = f.readlines()

        flines[4] = 'elapsed time for simu: \t \t '+ str(self.time_duration_simu)+'\n'
        with open(self.path_data_out + '/'+ self.filename_log + '.txt', 'w') as f:
            f.writelines(flines)
        return

class simu_year(Simu):

    def __init__():
        super.__init__(self, *args, **kwargs)

class Sig():
    '''
    instances for input-signals
    '''

    def __init__(self, name=None, file=None, ref_pth=None, p_col=None):
        self.name       = name
        self.filepath   = file
        self.ref_pth    = ref_pth
        self.props, self.df  = self.read_sig()
        self.normalized = False #
        self.p_col=p_col
        check_props = self.get_sig_properties()
        print('Props: ', self.props)
        self.rated_power = self.props['rated_power']
        self.starttime  = self.props['start_date']    # compare to df.date !!
        self.stoptime   = self.props['end_date']      # compare to df.date !!
        self.timerange  = None #?
        self.timediff   = int(self.props['time_incr'])# time diff between valus // in s |
        if not 'years' in self.props:
            self.years = check_props['p_years']


    def read_sig(self,):


        # decide wether to use date from pars or from df

        ### read specs
        #print(' +++SIG path: ', self.filepath)


        line_specs_end = sigauxf.find_line(self.filepath, 'end Sig - metadata')
        if line_specs_end:
            specs   = sigauxf.read_metadata(self.filepath, line_specs_end) # returns dict
            skprws  = line_specs_end + 3
            print('-->Specs: ', specs)
        else:
            specs   = None
            skprws  = None

        ### read data
        #data = None
        df = pd.read_csv(self.filepath, skiprows=skprws, header=[0])
        df['Date'] = pd.to_datetime(df['Date'])
        #df = df.set_index('Date') # in louter
        return specs, df#specs, data


    def get_sig_properties(self, ):
        # get startdate
        df_c    = self.df.copy()

        sd      = df_c.Date.min().strftime("%Y-%m-%d %H:%M:%S") # startdate
        ed      = df_c.Date.max().strftime("%Y-%m-%d %H:%M:%S") # enddate
        df_c['tdiff'] = df_c.Date.diff()

        years = df_c.Date.dt.year.drop_duplicates().to_list() # list of years in sig-df

        # get enddate
        # get list of years
        return {'p_startdate': sd, 'p_enddate': ed, 'p_timediff': df_c.tdiff, 'p_years': years}


    def unique_cols(cl):
        a = cl.to_numpy() # df.values (pandas<0.24)
        return (a[0] == a).all() #0)

    def check_sig_properties():
        '''
        check, if metadata consistent
        '''
        for item, value in prop_d:
            if value:
                pass
        if self.startdate != sd:                # Check consitency of if  startdate
            print('Non-consistent startdate')
            print('Metadata: ', self.startdate)
            print('df-data: ', sd)

        if self.enddate != ed:
            print('Non-consistent enddate')
            print('Metadata: ', self.enddate)
            print('df-data: ', ed)
        # get startdate

        # get enddate

        # time difference
        td = df_c.tdiff[1:] # timedifference
        if not unique_cols(td):
            print('Unhomogenous Time Range')
        # get list of years
        return
