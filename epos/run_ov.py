'''
+++ OLD VERSION +++

run only claculation of polarisation curve
options:
-save as csv
-plot
'''
print(__name__)

import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

from epos.main.simulation import ElSim
import epos.auxf.readingfiles as rf


def run_plr(scn_pth, T_spcs, i_spcs, p_spcs, mode=None, plot=True):

    sim = ElSim(scn_pth, full_simu=False)
    sim.setup_sim()
    #sim.logger.info(f'Scen_name: {sim.name}')
    sim.logger.info('Run plr test... mode= %s', mode)

    ### import parameters and clc_module
    #prms = rf.read_json_file(re_pth=scn_pth) # |dict
    #ver = sim.prms['bsc_par']['clc_ver']
    #tec = sim.prms['bsc_par']['tec_el'].lower()
    #sim.plr_clc     = impm('clc.' +tec+ '_plr_' + ver['plr'])
    #sim.flws_clc     = impm('clc.' +tec+ '_flws_' + ver['flws'])
    p_in = setup_p(p_spcs) # Returns list of namedtuples

    [i_s, i_e, inum] = i_spcs
    i_arr = np.linspace(i_s, i_e, inum)

    T_in = T_spcs

    tl = ['mlt', 'multi']
    if (not mode) or ('sngl' in mode):
        res_ca, res_an, res_cell, res_Urev, res_Utn = plr_sngl(sim, T_in, i_arr, p_in)
    elif any(mode in ele for ele in tl):
        i_arr, res = run_plr_mlt()
    else:
        print(f'No valid mode input (sys.argv = {mode})...')

    if plot:
        plt_bsc_res(res_cell, T_in, i_arr, p_in)
        clc_plt_eff(res_cell, res_Urev, res_Utn, i_arr)
        clc_plt_Pcell(res_cell, i_arr)
    return

def clc_plt_eff(res_arr, resurev, resutn, i_in):
    '''
    calc. theor. cell-efficiency

    only for single pressure and temp.
    '''
    u_in = res_arr[0,-1,:]
    urev = resurev[0,-1,:]
    utn = resutn[0,-1,:]
    eta_LHV = np.zeros(len(u_in))
    eta_HHV = np.zeros(len(u_in))

    for k in range(len(u_in)):
        #Vdot = i * A_cell * n_cells / (2*F)
        #e_spec_H2 = 3 # kWh/m³ (STP) (33.3 kWh/kg)
        #E_H2 = Vdot * e_spec_H2
        #E_in = i * A_cell * n_cells * u_cell
        eta_LHV[k] = urev[k] / u_in[k]
        eta_HHV[k] = utn[k] / u_in[k]

    plt.figure(111)
    plt.plot(i_in, eta_LHV, label='eta_LHV')
    plt.plot(i_in, eta_HHV, label='eta_HHV')
    plt.legend()
    plt.show()
    return

def clc_plt_Pcell(res_arr,i):
    u_in = res_arr[0,-1,:]
    P_cell = np.zeros(len(u_in))
    P_cell = u_in*i
    print('P_cell max: ', P_cell.max())
    plt.figure(222)
    plt.plot(i, P_cell/P_cell.max(), label='P_cell')#'$P_{cell}^{nrm}$')
    #plt.plot(i_in, eta_HHV, label='eta_HHV')
    plt.legend()
    plt.show()
    return

def plr_sngl(sim, T_in, i_in, p_in):
    '''
    process plr-function for single parameter set
    '''

    li = len(i_in)
    lT = len(T_in)
    lp = len(p_in)

    res_cell = np.zeros((lp, lT, li))
    res_ca = np.zeros((lp, lT, li))
    res_an = np.zeros((lp, lT, li))
    res_Urev = np.zeros((lp, lT, li))
    res_Utn = np.zeros((lp, lT, li))
    res_Pcell = np.zeros((lp, lT, li))

    for j in range(lp):

        for k in range(lT):

            for m in range(li):
                sim.clc_m.aux.clc_auxvals(sim, T_in[k])
                results_plr= sim.clc_m.plr.voltage_cell(sim, sim.pec, T_in[k], i_in[m], p_in[j], pp=None, ini=True)
                res_Urev[j,k,m] = sim.clc_m.plr.cv_rev(sim, sim.pec, T_in[k], p_in[j])[1]
                res_Utn[j,k,m] = sim.clc_m.plr.cv_rev(sim, sim.pec, T_in[k], p_in[j])[2]
                res_cell[j,k,m]= results_plr[-1]
                res_ca[j,k,m] = results_plr[0]
                res_an[j,k,m] = results_plr[1]
                #res_Pcell[j,k,m] = results_plr[-1]*i_in[m]
    return res_ca, res_an, res_cell, res_Urev, res_Utn

def plr_mlt():
    '''
    process plr-function for multiple parameter sets
    '''

    #for ...:
    #    run_plr_sng()

    return

def setup_p(p_spcs):
    lp_ca = len(p_spcs[0]) # length of pressure array
    lp_an = len(p_spcs[1]) # length of pressure array
    if lp_an > lp_ca:
        lp = lp_an
    else:
        lp = lp_ca

    if lp != lp_ca:
        p_ca_sngl = True
    elif lp != lp_an:
        p_ca_sngl = True
    else:
        p_ca_sngl = False
        p_an_sngl = False

    pnt = namedtuple('pnt', 'anode, cathode')

    p0_in = []
    for j in range(lp):
        if p_ca_sngl:
            p0_ca = p_spcs[0][0]
            p0_an = p_spcs[1][j]
        elif p_an_sngl:
            p0_ca = p_spcs[0][j]
            p0_an = p_spcs[1][0]
        else:
            p0_ca = p_spcs[0][j]
            p0_an = p_spcs[1][j]
        p0_in.append(pnt(p0_an, p0_ca)) # Append namedtuple to list
        #p0_in.append([p0_ca, p0_an])
    return p0_in

def plt_bsc_res(results, T_in, i_in, p_in):
    '''
    plot temp and pressure behaviour in two plots
    '''
    li = len(i_in)
    lT = len(T_in)
    lp = len(p_in)

    temp_res = results[0,:,:].reshape(lT,li)
    press_res = results[:,0,:].reshape(lp,li)

    lbs = [T_in, p_in]
    for n,(res, l) in enumerate(zip([temp_res, press_res], [lT,lp])):
        plt.figure(n)
        if n == 0:
            lb = 1
        if n == 1:
            lb =0
        for k in range(l):
            plt.plot(i_in, res[k,:], label=str(lbs[lb][0])+'_'+str(lbs[n][k]))
        plt.grid()
        plt.legend()
    plt.show()
    return

def plt_Pcell():
    return

def plt_bsc_res3d(res):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ll=len(r_lst)
    li=len(i_arr)
    x = np.array([[k]*li for k in range(1,ll+1)])#[np.arange(1,11,1)*10] #10 * np.outer(np.cos(u), np.sin(v))
    print('x.shape: ', x.shape)
    y = i_arr
    print('y.shape: ', y.shape)
    z = z_arr #par_arr[:,0,:]#.reshape(10,100)
    print('z.shape: ', z.shape)
    ax.plot_wireframe(x,y,z)#
    ax.set_xticks(np.arange(1,ll+1,1)) # plot all
    ax.set_xticklabels(prnms, rotation='vertical', fontsize=10) # plot all

    #ax.set_ylabel('Current density $ /~ A/cm²$')
    #ax.set_zlabel('Sensitivity coefficient')
    plt.show()
    return

if __name__ == '__main__':

    if (len(sys.argv) >1):
        if ('ael' in sys.argv[1].lower()):
            tec = 'ael'
        elif ('pem' in sys.argv[1].lower()):
            tec = 'pem'
    else:
        print('-- run default --')
        tec = 'pem'
    print('tec: ', tec)

    if len(sys.argv)>2:
        mod = sys.argv[2]
    else:
        mod=None

    if tec == 'pem':
        #scn_pth = 'data/scen/test/20210225/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'# Scenario name
        #scn_pth = 'data/scen/test/20210325/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'
        #scn_pth = 'data/scen/test/20211005/?
        scn_pth = 'data/scen/test/20211125/Scen_PEM_0.6_1_sig_03_syn_bump_v22_60e3_2000_00_.json'
        i_stp = 3.0*1e4 # Stop Current density // in A/m²
    elif tec == 'ael':
        #scn_pth = 'data/scen/test/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'
        scn_pth = 'data/scen/test/20210916/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'
        i_stp = 0.5*1e4 # Stop Current density // in A/m²

    T_spcs = [333,353] # fixed temperature(range)
    p0_ca = [101325, 101325*2, 101325*10] # fixed pressure at cathode
    p0_an = [101325, 101325*2, 101325*10] # fixed pressure at anode
    p_spcs = [p0_ca, p0_an]

    i_stt = 10 # start Current density
    i_num = 1000
    i_spcs = [i_stt, i_stp, i_num]


    run_plr(scn_pth, T_spcs, i_spcs, p_spcs, mode=mod)
