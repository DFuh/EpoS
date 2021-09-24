'''
test pem_flws

'''
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from epos.main.simulation import ElSim
import epos.auxf.readingfiles as rf
import epos.clc.pem_plr_v20 as plr
import epos.clc.pem_flws_v22 as flws

#scn_pth = 'data/scen/test/20210325/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'
scn_pth = 'data/scen/test/20210526/Scen_PEM_0.6_1_sig_05_WEAoff_2000_01_.json'
sim = ElSim(scn_pth, full_simu=False)
sim.setup_sim(test=True, testmode='flws') # testmode only for refvals
sim.logger.info(f'Scen_name: {sim.name}')
sim.logger.info('Run plr ... ')

l = 1000

T = 353 # // in K
t_arr = np.arange(l)
#i_arr = np.ones(l) * 25000#
i_arr = np.linspace(1,30000, l) # // in A/m²
#i_arr[0] = 0
m_H2O_in_an = 2*0.0033# in this case: water flow in 1 cell // in kg/s

c = np.zeros((l,4))
n = np.zeros((l,4))
xH2O = np.zeros((l,2))
perm = np.zeros((l,2))
c_H = np.zeros((l,1)) # Henry

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
    perm[j,:] = sim.av.n_H2_prm, sim.av.n_O2_prm
    xH2O[j,:] = sim.av.x_H2O_ch_ca, sim.av.x_H2O_ch_an
    c_H[j,:] = sim.av.c_H2_henry
    #print('sim.av.c_H2_henry: ',sim.av.c_H2_henry)

dat_lst = [c,n, xH2O,perm,c_H]
dat_len = [4,4,2,2,1]
nms = [['c_H2_out_an', 'c_H2_out_ca', 'c_O2_out_an', 'c_O2_out_ca'],
        ['n_H2_out_an', 'n_H2_out_ca', 'n_O2_out_an', 'n_O2_out_ca'],
        ['x_H2O_ca', 'x_H2O_an'],['n_H2_prm','n_O2_prm'],['c_H2_henry',]]


fig,axes = plt.subplots(len(dat_lst)+1,1,sharex=True)
#ax0, ax1 = axes
for k,ax in enumerate(axes):
    if k < len(dat_lst):
        for m in range(dat_len[k]):#enumerate(dat_lst[k]):
            ax.plot(t_arr, dat_lst[k][:,m], label=nms[k][m])
    else:
        ax.plot(t_arr,i_arr, label='i')
    ax.legend()
    ax.grid()
plt.show()

fig, ax = plt.subplots(1,1)

ax.plot(i_arr, (n[:,0]/ (n[:,0]+n[:,2])), label='X_H2inO2')
ax.legend()
ax.grid()
plt.show()
