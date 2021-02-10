'''
handle parameters / scenario_settings
'''
import os
import glob
import itertools

import epos.aux.handlingfiles as hf
import epos.aux.readingfiles as rf
import epos.aux.writingfiles as wf

import epos.clc.bsc as bc

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
    fin_dct = {}
    # basic par
    #nm_keys = ['tec_el', 'scl_el', 'rpow_el','sig', 'tec_el', 'rpow_ee']
    #nm_lst =[]
    #for key,val in dct_in.items():
    #    if key in nm_keys:
    #        nm_lst.append(val)

    #fin_dct['scen_name'] = create_name(nm_lst)
    fin_dct['scen_name'] = dct_in['scen_name'] # name of scenario
    fin_dct['scen_filename'] = dct_in['filename'] # name of scenario with file-suffix
    fin_dct['scen_dir'] = dct_in['scen_pth'] # path to directory of scen file
    fin_dct['scen_filepath'] = dct_in['flpth'] # full path to scen file
    # avoid duplicates in final dict:
    del dct_in['scen_name']
    del dct_in['filename']
    del dct_in['flpth']

    fin_dct['bsc_par'] = dct_in
    fin_dct['time_incr_clc'] = obj.sup_par['time_increment_clc']
    # paths
    fin_dct['pth_sup_par']      = obj.relpth_sup_par
    fin_dct['relpth_sig_par']   = obj.relpth_sig_par
    fin_dct['relpth_sig_data']  = obj.sig_par[dct_in['sig']]['path'] # sig name in dct_in -> returning path from sig-par-file
    fin_dct['refpth_sig_data']  = obj.sig_par[dct_in['sig']]['ref_path'] # reference path for output-dir-structure
    fin_dct['nm_pcol_sig']      = obj.sig_par[dct_in['sig']]['clmn_nm_p']
    fin_dct['searchkey_sig_metadata'] = obj.sig_par[dct_in['sig']]['searchkey_sig_metadata']
    fin_dct['metadata_sig'] = sig_mtd




    fin_dct['flpth_logfile']    = None # Needs to be generated in Simu-inst
    fin_dct['flpth_data_out']   = []
    fin_dct['bsc_pth_data_out'] = obj.sup_par['basic_path_data_output']
    # TODO: check following lines (20210127: inserted basename)
    fin_dct['refpth_out_data']  = hf.mirror_output_path(basename=obj.cwd, ref_pth=obj.sig_par[dct_in['sig']]['ref_path'],
                                                        bsc_pth_out=obj.sup_par['basic_path_data_output'],
                                                        tday=None, name=None)


    fin_dct['relpth_tec_parameters'] = dct_in['clc_ver']['tec_par']


    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    ### ini tec parameters
    flcntnt = rf.read_json_file(rel_pth=fin_dct['relpth_tec_parameters'])
    #if not flcntnt:
    #raise Exception('Could not read jsonfile; relpath: ', fin_dct['relpth_tec_parameters'])
    fin_dct['parameters_tec_el'] = flcntnt
    upd_dct = bc.clc_pwr_vals(obj, fin_dct['bsc_par'],fin_dct['parameters_tec_el'])
    fin_dct['parameters_tec_el'].update(upd_dct) # update tec_params
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#


    fin_dct['select_stored_values'] = obj.sup_par['select_stored_values']
    fin_dct['output_parameters'] = rf.read_json_file(basename=obj.cwd,
                                        rel_pth=obj.sup_par['output_parameters'][dct_in['tec_el']][0])['varkeys']
    #print('Fin Dict: ', fin_dct)

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
            while not os.path.exists(cpth[0]): # ???
                cpth.append(os.path.split(cpth))
            for pthi in cpth[::-1]:
                os.mkdir(pthi)
                print('Make new dir: ', pth)
        #flpth = os.path.join(pth, dct['filename'])
        flpth = os.path.join(pth, dct['scen_filename'])
        print('Filepath for storing: ', flpth)
        wf.write_to_json(flpth, dct)
    return

def mk_scen_filenames_and_paths(obj, version='00', prfx='Scen', sffx='.json'):
    '''
    partly duplicate of create_name() !
    '''

    key_lst=['tec_el', 'scl_el', 'rpow_el', 'sig', 'tec_ee', 'rpow_ee']
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
         [None] # par versions (tec_el) #self.parameters.clc_versions)
         ]
    print(__name__, '---> s: ', s)
    pssbl_sim = []
    instances = []

    k_lst = ['tec_el', 'scl_el', 'rpow_el', 'sig', 'tec_ee', 'rpow_ee', 'clc_ver'] # positional consistency with s
    v_klst = ['plr', 'flws', 'dgr' ,'pwr','thrm', 'tec_par']

    sim_dict = {}
    iter_lst = list( itertools.product(*s))

    print('Possible simulations:')
    for num, i_lst in enumerate( iter_lst):

        v_lst = []
        for key in v_klst:
            # extract possible clc versions by key (e.g. plr) and tec (i_lst[0])
            if key == 'tec_par':
                vers = obj.sup_par['parameter_files'][i_lst[0]]
            else:
                vers = obj.sup_par['version_clc_files_tec'][key][i_lst[0]]
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
                if not key == 'tec_ee':
                    sim_dict[nm][key] = val
                    if key == 'sig':
                        sig_nm = val
                else:
                    sim_dict[nm]['tec_ee'] = obj.sig_par[sig_nm]['gen_tec']
            sim_dict[nm]['clc_ver'] = dict.fromkeys(v_klst)
            # additional sub-dict for clc versions
            for n,key in enumerate(sim_dict[nm]['clc_ver'].keys()):
                sim_dict[nm]['clc_ver'][key] = v_i[n]

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
