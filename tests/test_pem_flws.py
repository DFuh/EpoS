'''
test pem_flws

'''
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from epos.main.simulation import ElSim
import epos.aux.readingfiles as rf
import epos.clc.pem_plr_v20 as plr
import epos.clc.pem_flws_v22 as flws

scn_pth = 'data/scen/test/20210325/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'

sim = ElSim(scn_pth, full_simu=False)
sim.setup_sim(test=True, testmode='flws')
sim.logger.info(f'Scen_name: {sim.name}')
sim.logger.info('Run plr ... ')

l = 300

T = 353 # // in K
i_arr = np.linspace(0,30000, l) # // in A/m²
m_H2O_in_an = 0.0033# in this case: water flow in 1 cell // in kg/s

c = np.zeros((l,4))
n = np.zeros((l,4))
c_in = (0,0,0,0)
n_in = (0,0,0,0)
for j,i_val in enumerate(i_arr):

    fout = flws.materialbalance(sim, T, i_val, m_H2O_in_an, sim.p, c_in, n_in,
                        stf=1, sns=True, ntd=None, m=None)
    c_res = fout.c_H2_out_an, fout.c_H2_out_ca, fout.c_O2_out_an, fout.c_O2_out_ca
    n_res = fout.n_H2_out_an, fout.n_H2_out_ca, fout.n_O2_out_an, fout.n_O2_out_ca
    #Y[j] = fout.x_H2_out_an

    c_in = c[j,:] = c_res
    n_in = n[j,:] = n_res

dat_lst = [c,n]
nms = [['c_H2_out_an', 'c_H2_out_ca', 'c_O2_out_an', 'c_O2_out_ca'],
        ['n_H2_out_an', 'n_H2_out_ca', 'n_O2_out_an', 'n_O2_out_ca']]
fig,axes = plt.subplots(len(dat_lst),1,sharex=True)
#ax0, ax1 = axes
for k,ax in enumerate(axes):
    for m in range(4):#enumerate(dat_lst[k]):
        ax.plot(i_arr, dat_lst[k][:,m], label=nms[k][m])
    ax.legend()
    ax.grid()
plt.show()
