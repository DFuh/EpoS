from collections import namedtuple
import time
from epos.main.simulation import ElSim
import epos.aux.readingfiles as rf
import epos.clc.gnrl_pwr_v20 as pwr

import epos.clc.ael_plr_v20 as plr
import epos.clc.ael_flws_v20 as flws

### Setup test sim-instance
#scn_pth = "data/scen/newtest/20210210/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json"
scn_pth = 'data/scen/newtest/20210211/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'

sim = ElSim(scn_pth, full_simu=False)
sim.setup_sim(test=True, testmode='plr')
sim.logger.info(f'Scen_name: {sim.name}')
sim.logger.info('Run plr ... ')

T = 353 # // in K
i = 1500 # // in A/m²
pp = flws.partial_pressure(sim, sim.pec, T, sim.p)


i_max = sim.pcll.current_density_max

P0 = 200#250 #200
t1 = time.time()
#pout = pwr.op_opt(sim, sim.pec, T, i, i_max, sim.p, pp, P_in=P0, u_mx=None, ifun=None, ini=False)
P_prv = 270#100#270
u_prv= 1.9
dt = 10
pout = pwr.cntrl_pow_clc(sim, sim.pec, T, i, sim.p, pp, P0, P_prv, u_prv, dt)

t2 = time.time()

print('Output: ',pout)
print(f'Elapsed time:  {t2-t1} seconds')
