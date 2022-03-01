'''
run only claculation of polarisation curve
options:
-save as csv
-plot
'''
print(__name__)

import sys
import os
import glob
import numpy as np
try:
    import matplotlib.pyplot as plt
    plots=True
except:
    plots=False
from collections import namedtuple

from epos.main.simulation import ElSim
import epos.auxf.handlingfiles as hf
import epos.auxf.readingfiles as rf


def run_plr(scn_pth, T_spcs, i_spcs, p_spcs, mode=None, plot=True,
            scn_dct=None, logpath=None):

    if plots:
        sim = ElSim(scn_pth, full_simu=False, scn_dct=scn_dct, logpath=logpath)
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
            plt_bsc_res(res_ca, T_in, i_arr, p_in,nfig=2201)
            plt_bsc_res(res_cell, T_in, i_arr, p_in)
            clc_plt_eff(res_cell, res_Urev, res_Utn, i_arr)
            clc_plt_Pcell(res_cell, i_arr)
        else:
            print(' --- Skip plr-plots ---')
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
        sim.p.anode, sim.p.cathode = p_in[j]
        for k in range(lT):
            n_in = (0,0,0,0)
            c_in = (0,0,0,0)
            for m in range(li):
                sim.clc_m.aux.clc_auxvals(sim, T_in[k])

                sim.av.fctr_n_c_A_abs = sim.pplnt.number_of_cells_in_plant_act * sim.pcll.active_cell_area
                sim.av.fctr_n_c_A_st = sim.pplnt.number_of_cells_in_stack_act * sim.pcll.active_cell_area
                m_ely = sim.av.rho_ely*sim.bop.volumetricflow_ely_nominal * sim.pplnt.number_of_cells_in_stack_act
                flws_o = sim.clc_m.flws.materialbalance(sim,T_in[k],  i_in[m],
                                    m_ely/sim.av.fctr_n_c_A_st,
                                    sim.p, c_in, n_in,
                                    stf=sim.av.fctr_n_c_A_abs,
                                    ntd=None, sns=False, m=m)
                n_in = flws_o.n_H2_out_an, flws_o.n_H2_out_ca, flws_o.n_O2_out_an, flws_o.n_O2_out_ca
                c_in = flws_o.c_H2_out_an, flws_o.c_H2_out_ca, flws_o.c_O2_out_an, flws_o.c_O2_out_ca
                
                pp = sim.p.pp_H2_mem_ca, sim.p.pp_O2_mem_an, 0 #pp_H2O
                # pp = pp = sim.clc_m.flws.partial_pressure(sim, sim.pec, T_in[k], sim.p)
                print('pp: ', pp)
                #results_plr= sim.clc_m.plr.voltage_cell(sim, sim.pec, T_in[k], i_in[m], p_in[j], pp=pp, ini=True)
                results_plr= sim.clc_m.plr.voltage_cell(sim, sim.pec, T_in[k], i_in[m], sim.p, pp=pp, ini=True)
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

def plt_bsc_res(results, T_in, i_in, p_in, nfig=None):
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
        if nfig is not None:
            plt.figure(nfig)
        else:
            plt.figure(n)
        if n == 0:
            lb = 1
        if n == 1:
            lb =0
        for k in range(l):
            plt.plot(i_in, res[k,:], label=str(lbs[lb][0])+'_'+str(lbs[n][k]))
        plt.grid()
        plt.xlabel('Current Density / $A/m²$')
        plt.ylabel('Cell Voltage / $V$')
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

    # print(os.path.abspath(os.path.join(pth, os.pardir)))
    pth0 = os.getcwd()
    print('cwd: ', pth0)
    skip = False
    if (len(sys.argv) >1):
        fllst0 = hf.lst_files_in_dir(sys.argv[1], bpth=pth0, suffix='.json')
        scn_pth = hf.select_file_from_filelist(fllst0)
    if scn_pth is None:
        skip = True
    ### Select tec
    tec = hf.select_single_item_from_list(['ael','pem'])

    if (tec == 'pem') and not skip:
        #scn_pth = 'data/scen/test/20210225/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'# Scenario name
        #scn_pth = 'data/scen/test/20210325/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'
        #scn_pth = 'data/scen/test/20211005/?
        # scn_pth = 'data/scen/test/20211125/Scen_PEM_0.6_1_sig_03_syn_bump_v22_60e3_2000_00_.json'
        i_stp = 4.0*1e4 # Stop Current density // in A/m²
    elif (tec == 'ael') and not skip:
        #scn_pth = 'data/scen/test/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'
        # scn_pth = 'data/scen/test/20210916/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'
        i_stp = 0.5*1e4 # Stop Current density // in A/m²
    else:
        print('No valid input ... skip ...')
        skip = True

    T_spcs = [333,353] # fixed temperature(range)
    p0_ca = [101325, 101325*2, 101325*10] # fixed pressure at cathode
    p0_an = [101325, 101325*2, 101325*10] # fixed pressure at anode
    p_spcs = [p0_ca, p0_an]

    i_stt = 10 # start Current density
    i_num = 1000
    i_spcs = [i_stt, i_stp, i_num]

    if not skip:
        mod=None
        run_plr(scn_pth, T_spcs, i_spcs, p_spcs, mode=mod)
