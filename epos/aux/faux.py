'''
...sonstiges...
'''
import os
import sys
import logging

from collections import namedtuple

from importlib import import_module as impm

def ini_logging(*obj, name=None, pth=None):
    print('obj in ini_logging: ', obj)
    if not obj:
        nm = name
    else:
        nm = str(obj[0].tdd)+'_'+obj[0].name
    logging.root.setLevel(logging.DEBUG)
    log = logging.getLogger('lggr_'+nm)
    formatter =logging.Formatter('[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s')#,

    if not pth:
        fh = logging.FileHandler(filename=nm+'.log')
    else:
        flnm = os.path.join(pth,nm+'.log')
        fh = logging.FileHandler(filename=flnm)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    log.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(formatter)
    log.addHandler(sh)

    return logging, 'lggr_'+nm



def ini_clc_versions(obj):

    #import

    ver = obj.prms['bsc_par']['clc_ver']
    tec = obj.prms['bsc_par']['tec_el'].lower()

    plr_clc     = impm('clc.' +tec+ '_plr_' + ver['plr'])
    flws_clc    = impm('clc.' +tec+ '_flws_' + ver['flws'])
    dgr_clc     = impm('clc.' +tec+ '_dgr_' + ver['dgr'])
    pwr_clc     = impm('clc.' +tec+ '_pwr_' + ver['pwr'])
    thrm_clc    = impm('clc.' +tec+ '_thrm_' + ver['thrm'])

    # namedtuple causing pickling error
    NT = namedtuple('NT', 'plr flws dgr pwr thrm')
    nt = NT(plr_clc, flws_clc, dgr_clc, pwr_clc, thrm_clc)
    #return nt
    return nt #plr_clc, flws_clc, dgr_clc, pwr_clc, thrm_clc


'''
def ini_tec_params(obj):
    par_dct = obj.prms['parameters_tec_el']
    Par = namedtuple('Par', par_dct)
    par = Par(**par_dct)
    return par
'''
