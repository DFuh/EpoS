'''
...sonstiges...
'''
import os
import sys
import logging

from collections import namedtuple

from importlib import import_module as impm

def ini_logging(*obj, name=None, pth=None, notest=True):
    print('obj in ini_logging: ', obj)
    if not obj:
        nm = name
    else:
        nm = str(obj[0].tdd)+'_'+obj[0].name
    if notest:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    #logging.root.setLevel(logging.DEBUG)
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

    return log, 'lggr_'+nm # old: return logging, 'lggr_'+nm



def ini_clc_versions(obj, prm_dct=None):

    #import
    if prm_dct:
        ver = prm_dct['clc_ver']
        tec = prm_dct['tec_el'].lower()
    else:
        ver = obj.prms['bsc_par']['clc_ver']
        tec = obj.prms['bsc_par']['tec_el'].lower()
    print('ver = ', ver)
    plr_clc     = impm('epos.clc.' +tec+ '_plr_' + ver['plr'])
    flws_clc    = impm('epos.clc.' +tec+ '_flws_' + ver['flws'])
    dgr_clc     = impm('epos.clc.' +tec+ '_dgr_' + ver['dgr'])
    #pwr_clc     = impm('epos.clc.' +tec+ '_pwr_' + ver['pwr'])
    pwr_clc     = impm('epos.clc.gnrl_pwr_' + ver['pwr'])
    # thrm_clc    = impm('epos.clc.' +tec+ '_thrm_' + ver['thrm'])
    thrm_clc    = impm('epos.clc.gnrl_thrm_' + ver['thrm'])
    strg_clc    = impm('epos.clc.strg_' + ver['strg'])
    aux_clc    = impm('epos.clc.auxf.' +tec+ '_aux_' + ver['aux'])

    # namedtuple causing pickling error
    NT = namedtuple('NT', 'plr flws dgr pwr thrm strg aux')
    nt = NT(plr_clc, flws_clc, dgr_clc, pwr_clc, thrm_clc, strg_clc, aux_clc)
    #return nt
    return nt #plr_clc, flws_clc, dgr_clc, pwr_clc, thrm_clc

def dyn_aux_import(flnm, nm):
    '''
    flnm => __file__
    nm => __name__
    '''
    sffx = os.path.basename(flnm).replace('.py','')
    prfx = nm.split(sffx)[0] +'auxf.'
    sffx = sffx.replace('_v', '_aux_v')
    mod = impm(prfx+sffx)
    return mod
'''
def ini_tec_params(obj):
    par_dct = obj.prms['parameters_tec_el']
    Par = namedtuple('Par', par_dct)
    par = Par(**par_dct)
    return par
'''
