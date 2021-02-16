'''
testing all plr modules on test scenario
-> check, if output in given range

test the following functions:
----
voltage_cell(obj, pec, T, i, p, pp=None, ini=False,
                appl_cv_rev=True, appl_ov_act=True, appl_ov_cnc=False,
                appl_ov_ohm=True)
-> ret: (U_an, U_ca, U_cell)

pwr_cell(u, i, A_cell=None, N_cells=None)
-> ret

cv_rev(obj, pec, T, pp)
-> ret: ((dE_rev_ca, dE_rev_an),dE_rev,U_tn)

ov_act(obj, pec, T, i, p, pp)
-> ret: (dU_act_an, dU_act_ca)

ov_cnc(apply_funct=True)
-> ret: (dU_cnc_ca, dU_cnc_an)

ov_ohm(obj, pec, T, i, pp)
-> ret: (dU_ohm_ca, dU_ohm_an, dU_ohm_sep)

### AUX
clc_conductivity_KOH(obj, pec, T, w_KOH)

clc_bubble_cvrg(obj, pec, T, i, p, pp)

clc_gibbs_free_energy(obj, pec, T)

'''
from collections import namedtuple
from epos.main.simulation import ElSim
import epos.aux.readingfiles as rf
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
i = 1000 # // in A/m²
pp = flws.partial_pressure(sim, sim.pec, T, sim.p)


#pwr = plr.pwr_cell(u, i, A_cell=None, N_cells=None)
#def test_voltage_cell():
cv_full = plr.voltage_cell(sim, sim.pec, T, i, sim.p, pp, ini=True)
cvrev = plr.cv_rev(sim, sim.pec, T, pp)
ovact = plr.ov_act(sim, sim.pec, T, i, sim.p, pp)
ovcnc = plr.ov_cnc(apply_funct=True)
ovohm = plr.ov_ohm(sim, sim.pec, T, i, pp)

def test_voltage_cell():
    errors = []
    if not (sim.refvals.U_an_lo <= cv_full[0] <= sim.refvals.U_an_hi):
        errors.append('U_an')
    if not (sim.refvals.U_ca_lo <= cv_full[0] <= sim.refvals.U_ca_hi):
        errors.append('U_ca')
    if not (sim.refvals.U_cell_lo <= cv_full[-1] <= sim.refvals.U_cell_hi):
        errors.append('U_cell')
    assert not errors, "errors occured in voltage_cell():\n{}".format("\n".join(errors))

def test_cv_rev():
    #Return of cv_rev(): ((dE_rev_ca, dE_rev_an),dE_rev,U_tn)
    errors = []

    #dE_rev_ca
    if not (sim.refvals.dE_rev_ca_lo <= cvrev[0][0] <= sim.refvals.dE_rev_ca_hi):
        errors.append('dE_rev_ca: '+str(cvrev[0][0]))
    # dE_rev_an
    if not (sim.refvals.dE_rev_an_lo <= cvrev[0][1] <= sim.refvals.dE_rev_an_hi):
        errors.append('dE_rev_an: '+str(cvrev[0][1]))
    # dE_rev
    if not (sim.refvals.dE_rev_lo <= cvrev[1] <= sim.refvals.dE_rev_hi):
        errors.append('dE_rev: '+str(cvrev[1]))
    # U_tn
    if not (sim.refvals.U_tn_lo <= cvrev[-1] <= sim.refvals.U_tn_hi):
        errors.append('U_tn: '+str(cvrev[-1]))

    assert not errors, "errors occured in cv_rev():\n{}".format("\n".join(errors))

der test_ov_act():
    # (dU_act_an, dU_act_ca)
    errors = []
    # dU_act_ca
    if not (sim.refvals.dU_act_ca_lo <= cvact[0] <= sim.refvals.dU_act_ca_hi):
        errors.append('dU_act_ca: '+str(cvact[0]))
    # dU_act_an
    if not (sim.refvals.dU_act_an_lo <= cvact[1] <= sim.refvals.dU_act_an_hi):
        errors.append('dU_act_an: '+str(cvact[1]))

    assert not errors, "errors occured in ov_act():\n{}".format("\n".join(errors))


def test_ov_cnc():
    # (dU_cnc_ca, dU_cnc_an)
    errors = []
    # dU_cnc_ca
    if not (sim.refvals.dU_cnc_ca_lo <= ovcnc[0] <= sim.refvals.dU_cnc_ca_hi):
        errors.append('dU_cnc_ca: '+str(ovcnc[0]))
    # dU_cnc_an
    if not (sim.refvals.dU_cnc_an_lo <= ovcnc[1] <= sim.refvals.dU_cnc_an_hi):
        errors.append('dU_cnc_an: '+str(ovcnc[1]))

    assert not errors, "errors occured in ov_cnc():\n{}".format("\n".join(errors))


def test_ov_ohm():
    # -> ret: (dU_ohm_ca, dU_ohm_an, dU_ohm_sep)
    errors = []
    # dU_ohm_ca
    if not (sim.refvals.dU_ohm_ca_lo <= ovohm[0] <= sim.refvals.dU_ohm_ca_hi):
        errors.append('dU_ohm_ca: '+str(ovohm[0]))
    # dU_ohm_an
    if not (sim.refvals.dU_ohm_an_lo <= ovohm[1] <= sim.refvals.dU_ohm_an_hi):
        errors.append('dU_ohm_an: '+str(ovohm[1]))
    # dU_ohm_sep
    if not (sim.refvals.dU_ohm_sep_lo <= ovohm[-1] <= sim.refvals.dU_ohm_sep_hi):
        errors.append('dU_ohm_sep: '+str(ovohm[-1]))

    assert not errors, "errors occured in ov_ohm():\n{}".format("\n".join(errors))

#def test_ov_act():

print('Output > voltage_cell < ', cv_full)
print('Output > cv_rev < ', cvrev)
print('Output > ov_act < ', ovact)
print(' -> Theta (an, ca): ', sim.av.theta_an, sim.av.theta_ca)
print('Output > ov_cnc < ', ovcnc)
print('Output > ov_ohm < ', ovohm)
