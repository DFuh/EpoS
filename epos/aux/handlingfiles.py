'''
handling of files
'''
import os
import aux.readingfiles as rf


def mk_relpath(filename=None, rel_parents= None):
    if rel_parents:
        return os.path.join(rel_parents,filename)
    else:
        return None

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
