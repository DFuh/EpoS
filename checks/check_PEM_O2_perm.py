'''
check implementation of oxygen permeation in PEM flws
'''


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from epos.main.simulation import ElSim
import epos.auxf.readingfiles as rf
import epos.clc.pem_flws_v22 as flws


#scn_pth = 'data/scen/test/20210325/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'
scn_pth = 'data/scen/test/20210522/Scen_PEM_0.6_1_sig_05_WEAoff_2000_00_.json'
sim = ElSim(scn_pth, full_simu=False)
sim.setup_sim()

T=353
i = 2e4     # // A/m²
A_cell = 300/1e4 # // m²
flws.xflws.clc_flws_auxpars(sim, T )
n_H2O_prm =flws.xflws.clc_crssflw_membrane_H2O_chandesris(sim,sim.pec, i, T, A_cell)
out = flws.xflws.clc_oxygen_permeation(sim, T, i)
print('out: ', out)
