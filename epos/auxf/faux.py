'''
...sonstiges...
'''
import os
import sys
import logging
import pandas as pd

from collections import namedtuple

from importlib import import_module as impm

def ini_logging(*obj, name=None, pth=None, notest=True):
    print('obj in ini_logging: ', obj)
    if not obj:
        nm = name
    else:
        nm = str(obj[0].today_ymdhs)+'_'+obj[0].name
    if notest:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    #logging.root.setLevel(logging.DEBUG)
    log = logging.getLogger('lggr_'+nm)
    formatter =logging.Formatter('[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s %(name)s \n- %(message)s')#,

    if not pth:
        pth=''
        # fh = logging.FileHandler(filename=nm+'.log')
    # else:
    lgpth = os.path.join(pth, 'logfiles')
    if not os.path.exists(lgpth):
        print('Make logpth: ',lgpth)
        os.makedirs(lgpth)
    flnm = os.path.join(lgpth,nm+'.log')
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
    # print('prm_dct: ', prm_dct)
    #import
    if prm_dct:
        ver = prm_dct['bsc_par']['clc_ver']
        tec = prm_dct['bsc_par']['tec_el'].lower()
        prms = prm_dct
    else:
        ver = obj.prms['bsc_par']['clc_ver']
        tec = obj.prms['bsc_par']['tec_el'].lower()
        prms = obj.prms
    print('ver = ', ver)
    plr_clc     = impm('epos.clc.' +tec+ '_plr_' + ver['plr'])
    flws_clc    = impm('epos.clc.' +tec+ '_flws_' + ver['flws'])
    dgr_clc     = impm('epos.clc.' +tec+ '_dgr_' + ver['dgr'])
    # dgr_clc     = impm('epos.clc.gnrl_dgr_' + ver['dgr'])
    #pwr_clc     = impm('epos.clc.' +tec+ '_pwr_' + ver['pwr'])
    pwr_clc     = impm('epos.clc.gnrl_pwr_' + ver['pwr'])
    # thrm_clc    = impm('epos.clc.' +tec+ '_thrm_' + ver['thrm'])
    thrm_clc    = impm('epos.clc.gnrl_thrm_' + ver['thrm'])
    strg_clc    = impm('epos.clc.strg_' + ver['strg'])
    aux_clc    = impm('epos.clc.auxf.' +tec+ '_aux_' + ver['aux'])

    #if obj.full_simu:
    ### Setup Degradation mode
    mode_dgr = prms['mode_dgr'].get(prms['bsc_par']['tec_el'].upper(), False)

    if (mode_dgr.lower() == 'lfun' and
            prms['bsc_par']['tec_el'].lower() == 'pem'):
        obj.logger.info('Degradation mode: %s', 'lfun')
        dgr_clc.voltage_increase = dgr_clc.voltage_increase_lfun
    elif (mode_dgr.lower() == 'lin' or
            prms['mode_dgr'][prms['bsc_par']['tec_el'].upper()].lower() == 'default'):
        obj.logger.info('Degradation mode: %s', 'default')
        dgr_clc.voltage_increase = dgr_clc.voltage_increase_lin
    else:
        obj.logger.info('No valid Degradation mode -> not applied')

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
def parameter_log_thrm(obj, flpth_out, par_a, par_b, idx=None):
    # par_a, par_b = pars
    header_a = 'C_cw,n_H2,n_O2,n_H2O_cns_in,n_H2O_resid_out,n_c,n_st,A_c,C_t,R_t,U_HAx,cp_mH2,cp_mO2,cp_mH2O,exp_f'.split(',')
    header_b = 'T_a,T_cwi,C_cw,n_H2,n_O2,n_H2O_cns_in,n_H2O_resid_out,U,i_cell,eta_e,n_c,n_st,A_c,C_t,R_t,U_HAx,dQ_heat,cp_mH2,cp_mO2,cp_mH2O,exp_f'.split(',')
    fl_a = flpth_out.replace('results', 'par_thrm_a')
    fl_b = flpth_out.replace('results', 'par_thrm_b')
    # print('flpth (par_a, thrm): ', fl_a)
    if os.path.exists(fl_a):
        hdr = False
    else:
        hdr=True
    pd.DataFrame(data=par_a, columns=header_a, index=idx).to_csv(fl_a, mode='a', header=hdr)
    pd.DataFrame(data=par_b, columns=header_b, index=idx).to_csv(fl_b, mode='a', header=hdr)

    return
