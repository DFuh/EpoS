'''
handling of files
'''
import os
import aux.readingfiles as rf
import aux.writingfiles as wr

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


def mk_relpath(filename=None, rel_parents= None):
    if rel_parents:
        return os.path.join(rel_parents,filename)
    else:
        return None

def mk_path(basename=None, rel_parents=None, filename=None):
    if not basename:
        basename = os.getcwd()
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

def mirror_output_path(ref_pth=None, bsc_pth_out='data/out', filepath=None, tday=None, name=None):
    pth = pth_mirror(ref_pth=ref_pth,
                        filepath = filepath,
                        bsc_pth_out = bsc_pth_out
                        )
    if tday:
        pth += '/' + str(tday)
    if name:
        pth += '/' + str(name)
    return pth


def pth_mirror(ref_pth='data/in', filepath=None, bsc_pth_out=None):
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
    if not filepath:
        out_pth = ref_pth
        subpath=''
    elif filepath:
        # split filepath
        filepath_head, filename = os.path.split(filepath)
        subpath = filepath_head.split(ref_pth)[1]
    #print('pathmirror: filepath: ', filepath)
    #print('pathmirror: filepath: ', filepath_head)

    out_path = os.path.join(bsc_pth_out + subpath)#+'/'+tday)
    while not os.path.isdir(out_path):
        out_path = input(f'{out_pth} non valid...\n Please insert valid path: ')

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
    metadata['sigtag']      = None  # Unique identifier of signal dataset
    metadata['df_cnt']      = str(n+1)+'/'+str(l)   # Counter of output files
    metadata['name']        = obj.name              # Name of Simu
    metadata['tec_el']      = obj.prms['bsc_par']['tec_el']                     # Electrolysis technology
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
    wr.write_to_csv(flpth, datasets=[metadata, df],
                    headline=headline,
                    footline=None,
                    data_headline=data_headline,
                    data_footline=data_footline )

    return
