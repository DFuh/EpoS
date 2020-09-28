'''
handling of files
'''
import os
import readingfiles as rf



def setup_paths(obj, ):
    '''
    setup all important paths
    '''

    pth_parameter_input =
    pth_scenario_file = 
    bsc_pth_data_input = concat_pths(obj.cwd, obj.sup_par['basic_path_data_input'])
    pth_data_input = concat_pths(obj.cwd, obj.sig_par['path'])
    pth_data_output =

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
