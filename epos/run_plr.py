'''
run only claculation of polarisation curve
options:
-save as csv
-plot
'''
import sys
import numpy as np
import matplotlib.pyplot as plt

from main.simulation import ElSim
import aux.readingfiles as rf


def run_plr(scn_pth, T_spcs, i_spcs, p_spcs, mode=None, plot=True):

    sim = ElSim(scn_pth, full_simu=False)
    sim.setup_sim()

    ### import parameters and clc_module
    #prms = rf.read_json_file(re_pth=scn_pth) # |dict
    #ver = sim.prms['bsc_par']['clc_ver']
    #tec = sim.prms['bsc_par']['tec_el'].lower()
    #sim.plr_clc     = impm('clc.' +tec+ '_plr_' + ver['plr'])
    #sim.flws_clc     = impm('clc.' +tec+ '_flws_' + ver['flws'])


    tl = ['mlt', 'multi']
    if (not mode) or ('sngl' in mode):
        res = run_plr_sngl(sim, T_spcs[0], i_spcs)
    elif any(mode in ele for ele in tl):
        i_arr, res = run_plr_mlt()
    else:
        print(f'No valid mode input (sys.argv = {mode})...')

    if plot:
        plt.figure()
        for j in range(lp):
            plt.plot(i_arr, res[j,-1,:])
        plt.legend()
        plt.show()

    return


def plr_sngl(sim, T_in, i_spcs, p_spcs):
    '''
    process plr-function for single parameter set
    '''
    [i_s, i_e, inum] = i_spcs
    i = np.linspace(i_s, i_e, inum)
    li = len(i)
    lo = 2 # lenght of plr output array

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

    res = np.zeros((lp,3,li))

    for j in range(lp):

        if p_ca_sngl:
            p0_in = p_spcs[0][0],p_spcs[1][j]
        elif p_an_sngl:
            p0_in = p_spcs[0][j],p_spcs[1][0]
        else:
            p0_in = p_spcs[0][j],p_spcs[1][j]

        for k in range(li):
            res[j,:,k] = sim.clc_m.plr_clc.voltage_cell(sim, sim.pec, T_in, i[k], p0_in, pp=None)

    return res

def plr_mlt():
    '''
    process plr-function for multiple parameter sets
    '''

    #for ...:
    #    run_plr_sng()

    return

if __name__ == '__main__':
    scen_path = 'data/scen/dftest/20201005/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'# Scenario name

    T_spcs = [333] # fixed temperature(range)
    p0_ca = [101325] # fixed pressure at cathode
    p0_an = [101325] # fixed pressure at anode
    p_spcs = [p0_ca, p0_an]

    i_stt = 0 # start Current density
    i_stp = 3.0*1e4 # Stop Current density // in A/m²
    i_num = 100
    i_spcs = [i_stt, i_stp, i_num]
    if len(sys.argv)>1:
        mod = sys.argv[1]
    else:
        mod=None

    run_plr(scen_path, T_spcs, i_spcs, p_spcs, mode=mod)
