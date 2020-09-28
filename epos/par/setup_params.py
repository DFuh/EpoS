'''
new structural approach:
- setup all parameters and dump into file

'''

'''
what about: ....

- https://pypi.org/project/dataclasses-json/ # not necessary, when using dataclass.asdict ?

'''

import aux.handlingfiles as hf
import aux.readingfiles as rf
import aux.writingfiles as wf

class ScenarioSetup():

    def __init__(self, *pths):
        for i,pth in enumerate(pths):
            if (len(pths) >0) & i==0:
                filepath_super_parameters = pth
            if (len(pths) >1) & i==1:
                filepath_sig_parameters = pth

        self.tdd            = datetime.datetime.now()
        self.today_ymd      = str(self.tdd.year) +str(self.tdd.month) + str(stef.tdd.day)
        self.today_ymdhs    = self.today_yms + str(self.tdd.hour) + str(self.tdd.minute)

        self.cwd = os.getcwd()

        ### read parameters
        # -> superior Pars
        self.sup_par = rf.read_json_file(filepath_super_parameters) # return dict

        # -> sig pars
        if not 'filepath_sig_parameters' in locals():
            filepath_sig_parameters = self.sup_par['filepath_sig_parameters']
        self.sig_par = rf.read_json_file(filepath_sig_parameters)

        # select scenarios
        ### select simulation-scenarios
        # clc versions
        # sig
        # tec
        # ee_power
        # el_power_fraction
        # scl  ---> scaling or absolute power val ?
        self.scen_dict = self.select_scenarios()

        #self.name = mk_scenario_name()

        self.scen_path = self.cwd + self.sup_par['basic_path_scenario_files'] + '/'+self.today_ymd
        self.scen_filename =



        wf.write_to_json(self.scen_filename, scen_dict)
        # -> set all paths




        # -> tec pars
        # -> par_data_output

        ### units ?

        ### scenario_name

    def make_scenario_files():
        '''
        extracting necessary parameters for simulation scenarios to store them in separate files
        '''
        

        return

    def select_scenarios(self,):
    #def ini_simu_instances(self, ):
        '''
        make list of possible simulation(scenarios)

        Parameters
        ----------
         : TYPE
            DESCRIPTION.

        Returns
        -------
        instances : TYPE list
            DESCRIPTION.

        '''


        s = [self.sup_par['tec_el'],
             self.sup_par['scl_el'],
             self.sup_par['nominal_power_ee'],
             self.sup_par['nominal_power_el'],
             list(self.sig_par.keys()),
             [None] #space holder for clc_ver
             #self.parameters.clc_versions)
             ]
        print(__name__, '---> s: ', s)
        pssbl_sim = []
        instances = []

        k_lst = ['tec', 'scl', 'rpow_ee', 'rpow_el', 'sig', 'clc_ver']
        v_klst = ['plr', 'flws', 'dgr' ,'pwr','thrm']

        sim_dict = {}
        iter_lst = list( itertools.product(*s))

        print('Find all possible simualtions:')
        for num, i_lst in enumerate( iter_lst):

            v_lst = []
            for key in v_klst:
                # extract possible clc versions by key (e.g. plr) and tec (i_lst[0])
                vers = self.parameters['version_clc_files_tec'][key][i_lst[0]]
                v_lst.append(vers)

            # combine allpossible versions-variations
            vi_lst = []
            for ele in list(itertools.product(*v_lst)):
                if ele not in vi_lst:
                    vi_lst.append(ele)

            # make simu-dict entry based on k_lst (names -> keys)
            for k,v_i in enumerate(vi_lst):
                nm = 'simu_'+str(num)+str(k)
                sim_dict[nm] = dict.fromkeys(k_lst)
                for key, val in zip(k_lst, i_lst[:-1]):
                    sim_dict[nm][key] = val
                sim_dict[nm]['clc_ver'] = dict.fromkeys(v_klst)
                # additional sub-dict for clc versions
                for n,key in enumerate(sim_[nm]['clc_ver'].keys()):
                    sim_[nm]['clc_ver'][key] = v_i[n]

        num = 0
        l = [3,8]
        title = '___'
        for k in k_lst:
            title += '_'+k+'_'*(l[1]-len(k)-1)
        for key, val in sim_dict.items():
            strnum = str(num)
            number = '['+strnum + '_'*(l[0]-len(strnum)+'] '

            for skey, sval in key,val:
                strkey = '['+skey + '_'*(l[1]-len(skey)+'] '
                strval = '['+skey + '_'*(l[1]-len(skey)+'] '
                print(number, strkey, strval)
            num +=1
        #    print('['+str(num)+']---> l: ', li)
        #    pssbl_sim.append(li)

        simu_dict_out = {}
        npt = input('Please select simulations to run: [input comma-separated integers or all -> enter]')
        if npt=='':
            pass_all = True
            idxlst = []
        else:
            pass_all = False
            try:
                idxlst = []
                for el in npt.split(','):
                    idxlst.append(int(el))
            except:
                idxlst = []
            print('idxlst: ', idxlst)

        for num,pss in enumerate( sim_):
            if (num in idxlst) or pass_all:
                simu_dict_out.update(pss)
            else:
                print('Skip: ['+str(num)+']---> l: ', pss)
        return simu_dict_out
