'''
handle parameters / scenario_settings
'''
import os
import itertools

import aux.handlingfiles as hf
import aux.readingfiles as rf
import aux.writingfiles as wf

def list_all_final_scen_dicts(obj):
    dict_lst = []
    for key, dct in obj.scen_dict.items():
        full_dct = mk_full_scenario_dict(obj,dct)
        dict_lst.append(full_dct)
    return dict_lst


def mk_full_scenario_dict(obj, dct_in):
    '''
    make full scenario dict
    '''
    fin_dct = {}
    # basic par
    nm_keys = ['tec_el', 'scl_el', 'rpow_el','sig', 'tec_el', 'rpow_ee']
    nm_lst =[]
    for key,val in dct_in.items():
        if key in nm_keys:
            nm_lst.append(val)

    fin_dct['scen_name'] = create_name(nm_lst)
    fin_dct['bsc_par'] = dct_in
    # paths
    fin_dct['pth_sup_par']      = obj.relpth_sup_par
    fin_dct['relpth_sig_par']   = obj.relpth_sig_par
    fin_dct['relpth_sig_data']  = obj.sig_par[dct_in['sig']]['path'] # sig name in dct_in -> returning path from sig-par-file
    fin_dct['refpth_sig_data']  = obj.sig_par[dct_in['sig']]['ref_path'] # reference path for output-dir-structure
    fin_dct['nm_pcol_sig']      = obj.sig_par[dct_in['sig']]['clmn_nm_p']
    fin_dct['key_end_metadata'] = obj.sig_par[dct_in['sig']]['key_end_metadata']


    fin_dct['flpth_logfile']    = None # Needs to be generated in Simu-inst
    fin_dct['flpth_data_out']   = []
    fin_dct['bsc_pth_data_out'] = obj.sup_par['basic_path_data_output']
    fin_dct['refpth_out_data']  = hf.mirror_output_path(ref_pth=obj.sig_par[dct_in['sig']]['ref_path'],
                                                        bsc_pth_out=obj.sup_par['basic_path_data_output'],
                                                        tday=None, name=None)

    fin_dct['relpth_tec_parameters'] = dct_in['clc_ver']['tec_par']
    fin_dct['tec_el_parameters'] = rf.read_json_file(rel_pth=fin_dct['relpth_tec_parameters'])

    print('Fin Dict: ', fin_dct)

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
        pth = os.path.join(obj.cwd, obj.pth_scen_files)
        print('Path for scen storage: ', pth)
        if not os.path.isdir(pth):
            print('Make new dir: ', pth)
            os.mkdir(pth)
        flpth = os.path.join(pth, dct['filename'])
        print('Filepath for storing: ', flpth)
        wf.write_to_json(flpth, dct)
    return

def mk_scen_filename(obj, version='', prfx='Scen', sffx='.json'):
    '''
    partly duplicate of create_name() !
    '''
    key_lst=['tec_el', 'scl_el', 'rpow_el', 'sig', 'tec_ee', 'rpow_ee']
    nm_lst = []
    for dct in obj.fin_scen_lst_o_dict:
        name = ''
        for key, val in dct['bsc_par'].items():
            print('Key: ', key)
            if key in key_lst:
                #nm_lst.append(val)
                name += '_' + str(val)
        fin_nm = prfx+'_'+name+'_'+ version +'_'+ sffx
        nm_lst.append(fin_nm)
        dct['filename'] = fin_nm
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
