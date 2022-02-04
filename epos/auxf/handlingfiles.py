'''
handling of files
'''
import os
import glob
import numpy as np
import pandas as pd
import epos.auxf.readingfiles as rf
import epos.auxf.writingfiles as wr


def mk_dir(full_path=None, add_suffix=None, no_duplicate=False):
    #mk_output_dir(full_path, add_suffix=None, no_duplicate=False):#, name=None, tday=None):
    '''
    make new output directory (if not existing)
    if tday is not none -> create subdir with tday
    if name is not none -> create subdir with name

    if suffix is not None -> suffix is added to dirnm (only for simu-specific dirnms)
    '''
    #if name:
    #    full_path += '/' + name
    #if tday:
    #    full_path += '/' + tday
    if add_suffix is not None:
        full_path = full_path + add_suffix
    if not os.path.exists(full_path):
        print('+++ make new dir: ', full_path)
        os.makedirs(full_path)
    else:
        if no_duplicate:
            full_path+'_'
            while os.path.exists(full_path):
                try:
                    num = int(fullpath[-1])+1

                    l = len(str(num))
                except:
                    num = 1
                    l=0
                full_path = full_path[:-l] + str(num)
            os.makedirs(full_path)
    return

def set_pths(obj, ):
    '''
    setup paths to relevant files


    -----
    relpth -> relative to usr-dir
    abspth -> absolute path on os
    epospth -> relative to .../EpoS/epos/

    -----
    Retruns
        dct_pths: dict with pths
    '''

    dct_pths = {}
    # Set paths

    dct_pths['pth_scen_out']    = obj.par_sup.get('relative_pth_scenario_files',
                                                    'epos/data/scen')
    dct_pths['pth_prms_sig']    = obj.par_sup.get('relative_pth_sig_parameters',
                                                    'epos/data/in')
    dct_pths['pth_data_out']    = obj.par_sup.get('relative_pth_data_out',
                                                    'epos/data/out')

    # make abspaths and check
    for key,val in dct_pths.items():
        ret = mk_abspath(obj, val, check_pth=True)
        setattr(obj,key,ret)
        val = ret
    # print('pth dict: ', dct_pths)
    return dct_pths

def mk_abspath(obj,pth='', abspth=None, flnm='',
                add_pth=None,
                check_flpth=False, check_pth=False):
    '''
    check, if epos-path or external ans join accordingly
    '''
    if add_pth is not None:
        pth += '/'+add_pth
    if flnm !='':
        pth+='/'+flnm
    if pth[0] == '/':
        pth = pth[1:]
    if abspth is None:
        ret = os.path.commonpath(['/epos', os.path.join('/',pth)])
        if pth[0] == '/':
            pth = pth[1:]
        if ret != '/':
            abspth = os.path.join(os.path.dirname(obj.dir_epos),pth)
        else:
            abspth = os.path.join(obj.dir_home, pth)
    else:
        abspth = os.path.join(abspth, pth)
    # print('abspath: ', abspth)
    if (check_flpth) and (flnm !='') and (not os.path.isfile(abspth)):
        if not os.path.isdir(abspth):
            # obj.logger.warning('Invalid path: %s', abspth)
            abspth=None
    if (check_pth) and (not os.path.isdir(abspth)):
        if not os.path.isfile(abspth):
            # obj.logger.warning('Invalid path: %s', abspth)
            abspth=None
    return abspth



def mk_relpath(filename=None, rel_parents= None):
    if rel_parents:
        return os.path.join(rel_parents,filename)
    else:
        return None

def mk_path(basename=None, rel_parents=None, filename=None):
    if not basename:
        #basename = os.getcwd()
        basename = Path(__file__).parents[1]
    if not filename:
        return basename
    if not rel_parents:
        os.path.join(basename, filename)
    else:
        return os.path.join(basename, rel_parents, filename)

def setup_paths(obj, ):
    '''
    setup all important paths
    '''

    pth_parameter_input = None
    pth_scenario_file = None
    bsc_pth_data_input = concat_pths(obj.cwd, obj.sup_par['basic_path_data_input'])
    pth_data_input = concat_pths(obj.cwd, obj.sig_par['path'])
    pth_data_output = None

    return

def concat_pths(*pths):

    '''
    make absolute paths from multiple components
    '''

    #while not os.path.isdir():
    res_pth = os.path.join(*pths)
    if not os.path.isdir(res_pth):
        res_pth = os.path.abspath(pths[-1])
    if not os.path.isdir(res_pth):
        res_pth = os.path.realpath(pths[-1])
    '''
    while not os.path.isdir(res_pth):
        npth = input(f'...non-valid path from {pths [-1]}...-> please insert valid path: ')
        res_pth = os.path.join(npth)
    '''
    return res_pth

def mirror_output_path(basename=None, ref_pth=None, bsc_pth_out='data/out', filepath=None, tday=None, name=None):
    pth = pth_mirror(basename=basename, ref_pth=ref_pth,
                        filepath = filepath,
                        bsc_pth_out = bsc_pth_out
                        )
    if tday:
        pth += '/' + str(tday)
    if name:
        pth += '/' + str(name)
    return pth


def pth_mirror(basename=None, ref_pth='data/in', filepath=None, bsc_pth_out=None):
    '''
    input:
    sig_filepath    | string | full input filepath
    bsc_pth         | string | ?
    ref_dir_pth     | string | reference path (subpath used for mirroring)

    create directory according to input file-location

    sig_filepath: string | relative path from
    bsc_pth: string | (rel) path of data

    #TODO: Code structure // redundant variable assignment
    '''
    print('basename:', basename)
    print('filepath:', filepath)
    print('bsc_pth_out:', bsc_pth_out)
    if not filepath:
        out_pth = ref_pth
        subpath=''
    elif filepath:
        # split filepath
        filepath_head, filename = os.path.split(filepath)
        subpath = filepath_head.split(ref_pth)[1]
        if subpath[0] == '/':
            subpath = subpath[1:]
    #print('pathmirror: filepath: ', filepath)
    #print('pathmirror: filepath: ', filepath_head)
    print('subpath:', subpath)
    # TODO: check the following lines !

    # out_path = os.path.join(basename, bsc_pth_out, subpath)#+'/'+tday)
    out_path = os.path.join(bsc_pth_out, subpath)#+'/'+tday)
    print('out_path (hf pth_mirror):', out_path)
    '''
    skip=False
    while (not os.path.isdir(out_path)) and (not skip):
        outpath = input(f'{out_path} not valid...\n Please insert valid path or skip(skp) ')
        if (outpath.lower() == 'skp') or (outpath.lower() == ''):
            skip = True
            out_path=None
    '''
    return out_path#, full_filepath

# ---------------------------------------------

def ini_logfile(obj, ):
    ### name of logfile

    ### path of logfile

    ### information on simulation
    # elapsed time
    # errors

    return


def check_for_duplicates():
    lst = glob.glob(flpth)
    print('existing df: ', lst)
    if lst:
        num = int(os.path.splitext(os.path.basename(pth))[0][-2:])
    return num

def get_line(filepth, search_text='end Simu - metadata', num_end=100):
    '''
    get line in csv, where ist says 'end info' (default)
    or any specified search text
    '''
    with open(filepth, 'r') as f:
        for num, line in enumerate(f,1):
            if search_text in line:
                return num
            if num > num_end:
                return None

def update_csv_file(obj, key, val, flpth=None):
    '''
    update existing json file
    '''
    if not flpth:
        flpth = os.path.join(obj.cwd, obj.path_data_out)

    #https://stackoverflow.com/questions/46126082/how-to-update-rows-in-a-csv-file

    return


def mk_output_file(obj, yr, n, l, flpth, df, dates):
    '''
    initialize output file:
    --- metadata ---
    specs of simu
    --- end metadata ---

    --- data ---
    header of data
    data of simu

    '''

    #simu_inst.auxcnt_yrs = # Auxilliary counter variable for simulation-years or periods

    metadata                = {}                    # Metadata of simu | dict |
    metadata['tag']         = obj.tag               # Unique identifier of simu
    metadata['ini_time']    = obj.tdd               # Time of initialization
    metadata['elapsed_time_simu']= None           # Placeholder f time passed during simulation
    metadata['elapsed_time_tot']= None           # Placeholder f time passed in total
    metadata['sigtag']      = None  # Unique identifier of signal dataset
    metadata['df_cnt']      = str(n+1)+'/'+str(l)   # Counter of output files
    metadata['name']        = obj.name              # Name of Simu
    metadata['tec_el']      = obj.prms['bsc_par']['tec_el']                     # Electrolysis technology
    metadata['el_pwr_nom']  = obj.prms['parameters_tec_el']['plant']['power_of_plant_act'].get('value',False)
    #                            obj.prms['parameters_tec_el']['plant']['power_of_plant_act'].get('value',False))#['nominal']
    metadata['el_pwr_ol']  = obj.prms['parameters_tec_el']['plant']['power_of_plant_overload']['value']#['max']
    metadata['tec_gen']     = obj.prms['metadata_sig']['generator_technology']  # Generator technology (signal; power-source)
    #metadata['years']       = []
    metadata['year']        = yr                    # Current year of data


    metadata['start_date']  = dates[0]                    #
    metadata['end_date']    = dates[2]
    #ts.strftime('%Y-%m-%d %H:%M:%S')
    #d0 =
    #df0 =
    #data        = None# Data-template (header, first col) of simu
    l0          = 10
    nm0         = 'Simu'
    headline    = ['-' * l0*2, wr.txt_symline(text=nm0)]
    #footline = ['-' * l0*2]

    metadata_hl     = wr.txt_symline(text='begin Simu - metadata')
    data_hl         = wr.txt_symline('begin Simu - data')
    data_headline   = [metadata_hl, data_hl]

    metadata_fl     = wr.txt_symline('end Simu - metadata')
    data_fl         = wr.txt_symline('end Simu - data')
    data_footline   = [metadata_fl, None]

    #flpth = simu_inst.path_data_out+f'/dataframe_{yr}_{n}.csv'
    #flpth = obj.flpth_data_out
    #print('flpth: ',flpth)
    lognm = os.path.basename(os.path.dirname(flpth))+'/'+os.path.basename(flpth)
    obj.logger.info('Setup output file: \n --> %s', lognm)
    wr.write_to_csv(flpth, datasets=[metadata, df],
                    headline=headline,
                    footline=None,
                    data_headline=data_headline,
                    data_footline=data_footline )

    return

def rewrite_output_files(obj, fllst, md_lst, df_lst):
    for fl, md, df in zip(fllst, md_lst, df_lst):
        # print('rewrite file: ', fl)
        ### following line creates output-file
        l0          = 10
        nm0         = 'Simu'
        headline    = ['-' * l0*2, wr.txt_symline(text=nm0)]
        #footline = ['-' * l0*2]

        metadata_hl     = wr.txt_symline(text='begin Simu - metadata')
        data_hl         = wr.txt_symline('begin Simu - data')
        data_headline   = [metadata_hl, data_hl]

        metadata_fl     = wr.txt_symline('end Simu - metadata')
        data_fl         = wr.txt_symline('end Simu - data')
        data_footline   = [metadata_fl, None]

        #flpth = simu_inst.path_data_out+f'/dataframe_{yr}_{n}.csv'
        #flpth = obj.flpth_data_out
        #print('flpth: ',flpth)
        lognm = os.path.basename(os.path.dirname(fl))+'/'+os.path.basename(fl)
        obj.logger.info('ReWrite output file: \n -> %s |', lognm)
        wr.write_to_csv(fl, datasets=[md, df],
                        headline=headline,
                        footline=None,
                        data_headline=data_headline,
                        data_footline=data_footline,
                        df_index=True )
    return

def final_fl_to_df(fllst):
    '''
    read final .csv-files in order to append data
    return full df (of one year)
    '''
    df_lst = []
    md_dct_lst = []
    for fl in fllst:
        line_metad = get_line(fl, 'begin Simu - data')
        line_units = get_line(fl, '>units<')
        # print(f'line_metad: {line_metad}, line_units: {line_units}')
        skprws = list(np.arange(line_metad))
        if line_units is not None:
            skprws = skprws+[line_units-1]
        # print(f'skprws: {skprws}')
        df_rd = pd.read_csv(fl, skiprows=skprws, nrows=30)
        # print(df_rd.head(3))
        df_lst.append(pd.read_csv(fl, skiprows=skprws))

        md_dct_lst.append(rf.read_metadata(fl,line_metad))
    return df_lst, md_dct_lst

def update_metadata(dct_lst, dt0, dt1):
    for dcti in dct_lst:
        dcti['elapsed_time_simu'] = dt0
        dcti['elapsed_time_tot'] = dt1
    return dct_lst

def store_simu_params(self, ):
    '''
    Init dict for storing actual parameters (to avoid confusion, in case of changed Scenario files)
    '''
    str_par_dct = {}
    str_par_dct['tag_sim'] = str(self.tag)
    str_par_dct['tag_sig'] = self.prms['metadata_sig'].get('tag',None)
    str_par_dct['parameters'] = self.prms
    pth, flnm0 = os.path.split(self.lst_pths_out[0])

    # for key,val in str_par_dct.items(): (error was json cant print timestamp)
    #     if isinstance(val, pd.Timestamp):
    #         str_par_dct[key] = val.strftime
    testdct = {}
    for key,val in str_par_dct['parameters']['metadata_sig'].items():
        #if isinstance(val, int):
        # print(f'key: {key}, val = -, type = {type(val)}')
        testdct[key]=val
        #print(testdct)
        flnm = flnm0.split('results')[0]+'parameters.json'
        wr.write_to_json(os.path.join(pth, flnm), testdct)
    #flnm = flnm0.replace('results', 'parameters').replace('.csv','.json')
    flnm = flnm0.split('results')[0]+'parameters.json'
    wr.write_to_json(os.path.join(pth, flnm), str_par_dct)
    return

def lst_files_in_dir(pth, bpth=None, suffix='.csv'):
    pth0 = bpth # os.getcwd()
    print('bpth: ', pth0)
    pth1 = 'data/scen'
    # pth1 = 'Scen_PEM_34170_syn_bump_v22_60e3_par_nopow_vrs002_'
    flpth = os.path.join(pth0,pth1,pth+'/*'+suffix)

    #spth = flpth+'/*'+suffix
    print('Searching files : ', flpth)
    fllst0 = glob.glob(flpth)

    return fllst0

def select_file_from_filelist(fllst, ):
    for i,file in enumerate(fllst):
        print(f'[{i}] -> ', file)
    if len(fllst)>0:
        nvi = False
        skip=False
    else:
        nvi=True
        skip=True
    while not nvi:
        slct0 = input('Please insert index to file: ')
        if slct0 == '':
            print('---Skip---')
            nvi=True
            skip=True
        else:
            skip=False
            try:
                slct = int(slct0)
                nvi=True
            except:
                print('No valid input (must be int)')
                nvi=False
            #print('slct = ', slct)
    if not skip:
        fl = fllst[slct]
    else:
        fl = None
    return fl

def select_single_item_from_list(lst, ):
    for i,item in enumerate(lst):
        print(f'[{i}] -> ', item)
    nvi = False
    while not nvi:
        slct0 = input('Please insert index to item: ')
        if slct0 == '':
            print('---Skip---')
            nvi=True
            skip=True
        else:
            skip=False
            try:
                slct = int(slct0)
                nvi=True
            except:
                print('No valid input (must be int)')
                nvi=False
    if not skip:
        out = lst[slct]
    else:
        out = None
    return out

################################################################################
def mk_dir_recursively(pth):
    '''
    make dir(s) in path, if not existing
    '''
    print('Path for scen storage: ', pth)
    parents_lst = []
    while not os.path.isdir(pth):
        print('Make new dir: ', pth)
        inp = input('Do it? (Enter -> yes / any other key -> skip)')
        if inp == '':
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
    return
