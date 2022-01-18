'''
handle parameters / scenario_settings
'''
import os
import glob
import itertools

import numpy as np
import pandas as pd
import scipy.signal as scys

import epos.auxf.handlingfiles as hf
# import epos.auxf.handlingdata as hd
import epos.auxf.readingfiles as rf
import epos.auxf.writingfiles as wf

import epos.clc.bsc as bc
from epos.main.simulation import ElSim
from epos.main import louter

def list_all_final_scen_dicts(obj):
    print('Metadata_sig_dicts: ', obj.metadata_sig_dicts)
    dict_lst = []
    for key, dct in obj.scen_dict.items():
        full_dct = mk_full_scenario_dict(obj,dct, obj.metadata_sig_dicts[key])
        dict_lst.append(full_dct)
    return dict_lst


def mk_full_scenario_dict(obj, dct_in, sig_mtd):
    '''
    make full scenario dict
    '''
    # TODO write this properly !!!

    fin_dct = {}
    # basic par
    #nm_keys = ['tec_el', 'scl_el', 'rpow_el','sig', 'tec_el', 'rpow_ee']
    #nm_lst =[]
    #for key,val in dct_in.items():
    #    if key in nm_keys:
    #        nm_lst.append(val)

    #  print('=====>>>> dict_in: \n', dct_in)
    #  print('=====>>>> sup_par: \n', obj.sup_par)

    #fin_dct['scen_name'] = create_name(nm_lst)
    fin_dct['creation_date'] = obj.today_ymdhs
    fin_dct['scen_name'] = dct_in['scen_name'] # name of scenario
    fin_dct['scen_filename'] = dct_in['filename'] # name of scenario with file-suffix
    fin_dct['scen_dir'] = dct_in['scen_pth'] # path to directory of scen file
    fin_dct['scen_filepath'] = dct_in['flpth'] # full path to scen file
    fin_dct['external_dir_data_output'] = obj.sup_par['external_dir_data_output'] # path to external data dir
    fin_dct['nms_data_extraction'] = obj.sup_par['nms_data_extraction'] #
    fin_dct['keys_data_extraction'] = obj.sup_par['keys_data_extraction'] #
    fin_dct['mode_dgr'] = obj.sup_par['degradation_mode'] #

    # avoid duplicates in final dict:
    del dct_in['scen_name']
    del dct_in['filename']
    del dct_in['flpth']

    fin_dct['bsc_par'] = dct_in
    fin_dct['time_incr_clc'] = obj.sup_par['time_increment_clc']
    # paths
    fin_dct['pth_sup_par']      = obj.relpth_sup_par
    fin_dct['relpth_sig_par']   = obj.relpth_sig_par
    fin_dct.update({key:val for (key, val) in obj.sig_par[dct_in['input']].items()})
    # fin_dct['relpth_sig_data']  = obj.sig_par[dct_in['input']]['path_sig'] # sig name in dct_in -> returning path from sig-par-file
    # fin_dct['refpth_input_data']  = obj.sig_par[dct_in['input']]['ref_path'] # reference path for output-dir-structure
    # fin_dct['relpth_H2dmnd_data']  = obj.sig_par[dct_in['input']]['path_H2dmnd'] # sig name in dct_in -> returning path from sig-par-file
    # fin_dct['relpth_c_electr_data']  = obj.sig_par[dct_in['input']]['path_c_electr'] # sig name in dct_in -> returning path from sig-par-file
    # fin_dct['relpth_f_emiss_data']  = obj.sig_par[dct_in['input']]['path_f_emiss'] # sig name in dct_in -> returning path from sig-par-file

    # fin_dct['nm_col_sig']      = obj.sig_par[dct_in['input']]['clmn_nm_p']
    # fin_dct['nm_col_H2dmnd']      = obj.sig_par[dct_in['input']]['clmn_nm_H2dmnd']
    # fin_dct['searchkey_sig_metadata'] = obj.sig_par[dct_in['input']]['searchkey_sig_metadata']
    # fin_dct['searchkey_H2dmnd_metadata'] = obj.sig_par[dct_in['input']]['searchkey_H2dmnd_metadata']
    fin_dct['thrml_scaling'] = obj.sup_par['thrml_scaling']
    fin_dct['metadata_sig'] = sig_mtd

    fin_dct['storage_clc_iso'] = obj.sup_par['storage_clc_iso']
    fin_dct['storage_clc_dyn'] = obj.sup_par['storage_clc_dyn']



    fin_dct['flpth_logfile']    = None # Needs to be generated in Simu-inst
    fin_dct['flpth_data_out']   = []
    fin_dct['bsc_pth_data_out'] = obj.sup_par['basic_path_data_output']
    # TODO: check following lines (20210127: inserted basename)
    fin_dct['refpth_out_data']  = hf.mirror_output_path(basename=obj.cwd, ref_pth=obj.sig_par[dct_in['input']]['refpth_input_data'],
                                                        bsc_pth_out=obj.sup_par['basic_path_data_output'],
                                                        tday=None, name=None)


    fin_dct['relpth_tec_parameters'] = dct_in['clc_ver']['tec_par']
    fin_dct['relpth_strg_parameters'] = dct_in['clc_ver']['strg_par']

    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    ### ini tec parameters
    flcntnt = rf.read_json_file(rel_pth=fin_dct['relpth_tec_parameters'])
    strgpar = rf.read_json_file(rel_pth=fin_dct['relpth_strg_parameters'])
    #if not flcntnt:
    #raise Exception('Could not read jsonfile; relpath: ', fin_dct['relpth_tec_parameters'])
    fin_dct['parameters_tec_el'] = flcntnt
    # fin_dct['parameters_tec_el']['tec'] = bsc_par['tec_el']
    fin_dct['parameters_strg'] = strgpar
    #obj.prms = flcntnt



    ### clc power vals
    upd_dct = bc.clc_pwr_vls(obj, fin_dct,fin_dct['parameters_tec_el'])
    fin_dct['parameters_tec_el'].update(upd_dct) # update tec_params


    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#


    fin_dct['select_stored_values'] = obj.sup_par['select_stored_values']
    fin_dct['output_parameters'] = rf.read_json_file(basename=obj.cwd,
                                        rel_pth=obj.sup_par['output_parameters'][dct_in['tec_el']][0])['varkeys']
    #  print('=== >>> fin_dct[output_parameters]: ', fin_dct['output_parameters'])
    # print('Fin Dict: ', fin_dct)

    if True:
        # ===========================================
        # // // // // // // // // // // // // // // /
        # ===========================================
        ### tune PID-cntrl
        ### create instance for calc ()
        sim = ElSim(scn_setup=True, scn_dct=fin_dct)
        sim.setup_sim(test=True, testmode='scn_setup')

        # sim.df0 = hd.mk_df_data_output(obj, dates)
        full_lst = []
        for k,v in fin_dct['output_parameters'].items():
            full_lst.append(k)
        default_data = [0]*len(fin_dct['output_parameters'])
        sim.df0 = pd.DataFrame( data=[default_data], columns=full_lst) # full df
        sim.df0['T_st'] = 298
        #print(sim.df0.head(5))

        # mssflw_max = fin_dct['parameters_tec_el']['periphery']['massflow_coolant_max']['value']
        # sim.mssflw_max = mssflw_max

        # Ku = np.linspace(mssflw_max/5*0.1, mssflw_max/5, 20)
        Ku = 0.1#mssflw_max/5*0.1
        T_tar = fin_dct['parameters_tec_el']['cell']['temperature_nominal']['value']
        # what, if Temp > T_crit?
        sim.pid_ctrl.reset()

        Ku_lim = 0.2
        ratio = 0
        mr = 0
        k=0
        while (ratio< 0.5) and (mr < 0.75) and (Ku<Ku_lim):
            print('=======================================================')
            sim.pid_ctrl.tune_man(Ku, 0,0)
            print('PID-tuning: (hd): ', sim.pid_ctrl.get_tuning_par())
            print(f'PID-tuning... mr= {mr}  || ratio= {ratio}')
            #df = pd.DataFrame()
            df = louter.mainloop(sim, )
            print(df.tail(5))
            pks, _ = scys.find_peaks(df.T_st.to_numpy(),height=T_tar, threshold=0.5)
            ratio = len(pks)/(len(df)/2)
            #if pks:
            try:
                mr = max(pks)/len(df)
            except:
                mr = 0
            # binc = np.bincount(np.diff(pks))

            # k +=1
            if Ku<5:
                Ku += 0.1
            elif Ku>10:
                Ku += 1
            else:
                Ku += 0.5
        print(f'Found Ku (osci): Ku = {Ku}')
        # Ziegler-Nichols, classic
        # Kp = Ku*0.6
        # Ki = 1.2*Ku/20
        # Kd = 0.105*Ku*20
        # Ziegler-Nichols, no overshoot
        Kp = Ku*0.2
        Ki = 0.4*Ku/20
        Kd = 0.6666666*Ku*20
        fin_dct['parameters_tec_el']['periphery']['pid_parameters']["value"]=[
                Kp, Ki, Kd]

    ####### Scenario name (edit 20211209) #########
    lst_attrs = [fin_dct['bsc_par']['tec_el'],
                    round(fin_dct['parameters_tec_el']['plant']['power_of_plant_act']['value']),
                    fin_dct['bsc_par']['tec_ee'],
                    'par'+fin_dct['relpth_tec_parameters'].lower().split('par_'+fin_dct['bsc_par']['tec_el'].lower())[-1].replace('.json','')
                    ]
    fin_dct['scen_name'] = 'Scen_'+ create_name(lst_attrs)
    fin_dct['scen_filename'] = fin_dct['scen_name']+'.json'
    return fin_dct

def create_name(nmlst):
    if len(nmlst)>0:
        ostr = ''
        for nm in nmlst:
            if isinstance(nm, float): #replace decimal number for filename
                nm = str(nm).replace('.', '-')
            ostr += str(nm) +'_'
    else:
        ostr = None
    return ostr


######################################################
def store_scenario_files(obj):
    for dct in obj.fin_scen_lst_o_dict:
        #pth = os.path.join(obj.cwd, obj.pth_scen_files)
        pth = os.path.join(obj.cwd,dct['scen_dir'])
        print('Path for scen storage: ', pth)
        parents_lst = []
        while not os.path.isdir(pth):
            #print('Make new dir: ', pth)
            try:
                os.mkdir(pth)
                print('Make new dir: ', pth)
            except:
                splt = os.path.split(pth)
                parents_lst.append(splt[1])
                pth = splt[0]
        for prnt in parents_lst:
            npth = os.path.join(pth,prnt)
            cpth = [npth,]
            # TODO check below! edit_DF_20210210: [0]
            while not os.path.exists(cpth[-1]): # ???
                print('cpth 0: ', cpth[-1])
                cpth.append(os.path.dirname(cpth[-1]))
            print('cpth 1: ', cpth)
            for pthi in cpth[::-1][:-1]:
                if not os.path.exists(pthi): # Redundant /// UGLY !
                    os.mkdir(pthi)
                    print('Make new dir: ', pthi)
        #flpth = os.path.join(pth, dct['filename'])
        flpth = os.path.join(pth, dct['scen_filename'])
        print('Filepath for storing: ', flpth)
        #print('dct: ', dct)
        wf.write_to_json(flpth, dct)
    return

def mk_scen_filenames_and_paths(obj, version='00', prfx='Scen', sffx='.json'):
    '''
    partly duplicate of create_name() !
    '''

    key_lst=['tec_el', 'scl_el', 'rpow_el', 'input', 'tec_ee', 'rpow_ee']
    nm_lst = []

    pth = os.path.join(obj.sup_par['basic_path_scenario_files'], obj.today_ymd)
    flpth_lst = []
    for fl in glob.glob(os.path.join(obj.cwd,pth)+'/*.json'):
        flpth_lst.append(fl)
    print('flpth_lst: ', flpth_lst)
    for dct in obj.scen_dict: #fin_scen_lst_o_dict:
        name = ''
        print('Dct: ', dct)
        for key, val in obj.scen_dict[dct].items():#.items(): #['bsc_par'].items():
            print('Key: ', key)
            if key in key_lst:
                #nm_lst.append(val)
                name += str(val) + '_'
        fin_nm = prfx+'_'+name + version +'_'
        #nm_lst.append(fin_nm)
        #dct['scen_name'] = fin_nm
        #obj.scen_dict[dct]['scen_name'] = fin_nm
        #dct['filename'] = fin_nm+sffx
        flnm = fin_nm+sffx
        #obj.scen_dict[dct]['filename'] = fin_nm+sffx


        #pth = os.path.join(obj.cwd, obj.pth_scen_files)

        flpth = os.path.join(pth, flnm)
        fllflpth = os.path.join(obj.cwd, flpth)

        #if os.path.isfile(obj.scen_dict[dct]['flpth']):

        ### check for duplicates
        # TODO: check, if identical dict exists !

        print('flpth_lst in loop: ', flpth_lst)
        while fllflpth in flpth_lst:
            s = flnm.split('_')
            num = int(s[-2])+1
            if num <10:
                ns = str(0)+str(num)
            else:
                ns = str(num)
            s[-2] = ns
            flnm = '_'.join(s)
            #dct[]'filename'] = nns
            #dct['scen_name'] = nns.replace('.json','')
            #dct.update({"filename": nns})
            #dct.update({"scen_name": nns.replace('.json','')})
            fin_nm = flnm.replace('.json','')

            fllflpth = os.path.join(obj.cwd,pth, flnm) #dct['filename'])
            print('flnm: ', flnm)
        flpth_lst.append(fllflpth)

        obj.scen_dict[dct]['scen_name'] = fin_nm
        obj.scen_dict[dct]['filename'] = flnm
        obj.scen_dict[dct]['scen_pth']= pth
        flpth = os.path.join(pth, flnm)
        obj.scen_dict[dct]['flpth'] = flpth

    return

def select_scenarios(obj,):
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


    s = [obj.sup_par['tec_el'],
         obj.sup_par['scaling_el'],
         obj.sup_par['nominal_power_el'],
         list(obj.sig_par.keys()),
         [None], # space holder for tec_ee
         obj.sup_par['nominal_power_ee'],
         [None], #space holder for clc_ver
         [None], # par versions (tec_el) #self.parameters.clc_versions)
         [None] # par versions (storage)
         ]
    print(__name__, '---> s: ', s)
    # pssbl_sim = []
    # instances = []

    # TODO: remove/ replace  hardcoded lines below
    k_lst = ['tec_el', 'scl_el', 'rpow_el', 'input', 'tec_ee', 'rpow_ee', 'clc_ver'] # positional consistency with s
    # v_klst = ['plr', 'flws', 'dgr' ,'pwr','thrm', 'aux', 'strg', 'tec_par']
    v_klst = list(obj.sup_par['version_clc_files_tec'].keys())+['tec_par', 'strg_par']
    # print('obj.sup_par[version_clc_files_tec]: ', list(obj.sup_par['version_clc_files_tec'].keys()))
    sim_dict = {}
    iter_lst = list( itertools.product(*s))

    # print(iter_lst)

    print('Possible simulations:')
    for num, i_lst in enumerate( iter_lst):

        v_lst = []
        for key in v_klst:
            # extract possible clc versions by key (e.g. plr) and tec (i_lst[0])
            if key == 'tec_par':
                vers = obj.sup_par['parameter_files'][i_lst[0]]
            elif key == 'strg_par':
                vers = obj.sup_par['parameter_files']['strg']
            else:
                vers = obj.sup_par['version_clc_files_tec'][key][i_lst[0]]
            v_lst.append(vers)

        # combine allpossible versions-variations
        vi_lst = []
        for ele in list(itertools.product(*v_lst)):
            if ele not in vi_lst:
                vi_lst.append(ele)

        # print(vi_lst)
        # make simu-dict entry based on k_lst (names -> keys)
        for k,v_i in enumerate(vi_lst):
            nm = 'simu_'+str(num)+str(k)
            sim_dict[nm] = dict.fromkeys(k_lst)
            for key, val in zip(k_lst, i_lst[:-1]):
                if not key == 'tec_ee':
                    sim_dict[nm][key] = val
                    if key == 'input':
                        sig_nm = val
                else:
                    sim_dict[nm]['tec_ee'] = obj.sig_par[sig_nm]['gen_tec']
            sim_dict[nm]['clc_ver'] = dict.fromkeys(v_klst)
            # additional sub-dict for clc versions
            for n,key in enumerate(sim_dict[nm]['clc_ver'].keys()):
                sim_dict[nm]['clc_ver'][key] = v_i[n]
    # print(sim_dict)
    num = 0
    l = [3,8,5,5,9,8,8] # TODO: find proper solution
    title = '___'
    for n,k in enumerate(k_lst):
        title += '_'+str(k) + '_'*(l[1]-len(k)-1) + '_'*((n*3))
    print('------'+title)
    for key, val in sim_dict.items():
        strnum = str(num)
        spcs = '_' * (l[0]-len(strnum))
        number = '[' +strnum + spcs+'] '

        prnt_lst =[]
        numi=0
        for skey,sval in val.items():
            #strkey = '_'+skey + '_'*(l[numi]-len(skey))+'_'
            strval = '_'+str(sval) + '_'*(l[numi]-len(str(sval)))+'_'
            prnt_lst.append(strval)
            numi +=1
        print(number, prnt_lst)
        num +=1

    #    print('['+str(num)+']---> l: ', li)
    #    pssbl_sim.append(li)

    simu_dict_out = {}#dict.fromkeys(list(sim_dict.keys()))
    npt = input('Please select simulations to run: [input comma-separated integers or all -> enter]')
    if npt == '':
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

    for num,pss in enumerate( sim_dict):
        if (num in idxlst) or pass_all:
            print('pss: ',pss)
            simu_dict_out[pss] = sim_dict[pss]
        else:
            print('Skip: ['+str(num)+']---> l: ', pss)
    print('Dict with selected simu-pars: ', simu_dict_out)
    return simu_dict_out
