#from collections import namedtuple
#import time
from epos.main.simulation import ElSim
#import epos.aux.readingfiles as rf
import epos.clc.gnrl_pwr_v20 as pwr

import epos.clc.ael_plr_v20 as plr
import epos.clc.ael_flws_v20 as flws

from epos.clc import bsc
### Setup test sim-instance
#scn_pth = "data/scen/newtest/20210210/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json"
#scn_pth = 'data/scen/test/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'
scn_pth = 'data/scen/test/20210225/Scen_AEL_0.6_1_sig_05_WEAoff_2000_00_.json'

sim = ElSim(scn_pth, full_simu=False)
sim.setup_sim(test=True, testmode='plr')
#sim.logger.info(f'Scen_name: {sim.name}')
#sim.logger.info('Run plr ... ')

bsc_par = sim.prms['bsc_par']
par_dct = sim.prms['parameters_tec_el']
ret = bsc.clc_pwr_vls(sim, bsc_par, par_dct)
