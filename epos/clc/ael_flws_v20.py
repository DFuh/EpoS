'''
AEL
calculation: material balance in cells
'''
import os
import numpy as np
import importlib as impl
print(__name__ + ' imported...')

import epos.aux.faux as fx

xflws = fx.dyn_aux_import(__file__, __name__)

def testfun():
    res = xflws.testcalc(2,3)
    print('testfun --> result from testcalc:', res)
    return

def clc_matbal_params_Tbased(obj,):
    ### clc T-params
    obj.D_ik       = xflws.clc_D_ik(obj.T, obj.w_KOH) #* (1 + (sns * (fctr -1)))

    obj.Vhcell     = xflws.clc_Vhcell(obj.T, obj.w_KOH)

    obj.gamma      = xflws.clc_gamma(obj.T, obj.w_KOH)

    obj.rho_H2     = xflws.clc_rho_H2(obj.T, obj.w_KOH)

    obj.rho_KOH    = xflws.clc_rho_KOH(obj.T, obj.w_KOH)

    obj.beta       = xflws.clc_beta(fctr=obj.betafctr)

    obj.S_ik       = xflws.clc_S_ik(obj.T, obj.w_KOH)

    obj.pp_H2O     = xflws.clc_pp_H2O(obj.T, obj.w_KOH)

    return

def clc_matbal_params_ibased(obj,i):
    ### clc i-params

        ### clc flow of electrolyte (VL)
    clc_VL(obj, i=i,dynVely=obj.dyn_ely, balVely=obj.bal_ely)

    obj.epsilon    = xflws.clc_epsilon(obj, i)

    obj.d_b        = xflws.clc_d_b(obj, i, obj.i_crit_db, obj.gamma, obj.beta,
                                    obj.rho_H2, obj.rho_KOH)

    obj.k_Li       = xflws.clc_kLi(obj, obj.T, i, d_b=obj.d_b,
                                        D_ik=obj.D_ik, epsilon=obj.epsilon,
                                        fctr=obj.klifctr)

    obj.A_GL       = xflws.clc_A_GL(obj.Vhcell, obj.epsilon, obj.gamma, obj.d_b, obj.p_an, obj.p_ca)

    obj.f_G        = xflws.clc_f_G(obj, i, fctr=obj.fGfctr)
    return


def clc_VL(self, i=None, dynVely=False, balVely=True):
    '''
    calc electrolytre volumetric flowrate
    orig: parameter_sensibility_bsc_matbal_v1002.py
    '''
    def clc_dynVely(V0,i):
        #i_max   = 0.5#*1e4         #m_par.i_max

        #Vely_min = 0.2 * V0 #m_par.Vely_min
        #print('i / self.i_max: ', (i / self.i_max))
        #print(f'i= {i}, i_max={self.i_max}')
        V_is    = V0 * i / self.i_max
        if V_is < self.Vely_min:
            V_is = self.Vely_min
        return V_is

    #print('dynVely: ', dynVely)

    if not balVely:
        self.VL_an = self.V0_ely * self.fac_Vely_an
        self.VL_ca = self.V0_ely * self.fac_Vely_ca
    else:
        self.VL_an           = self.V0_ely # liquid (electrolyte) flowrate; anode // in m³/s | Haug 2017, Tab2: 330mL/min
        self.VL_ca           = self.V0_ely # liquid (electrolyte) flowrate; cathode // in m³/s

    if i and dynVely:
        self.VL_an = clc_dynVely(self.VL_an, i)
        self.VL_ca = clc_dynVely(self.VL_ca, i)
    self.VL = (self.VL_an + self.VL_ca) / 2 # for clc c_mix
    #print(f'VL = {self.VL}, VL_an= {self.VL_an}, VL_ca= {self.VL_ca}') #'??? line 155 in clcVL')
    return


def materialbalance(self, T, i, c_in, n_in): #sns=False):

    t = 0 # dummy variable (for matbal; inc ase of dyn. calc)
    ### pre-clc for matbal
    xflws.matbal_preclc(self, T,i, )

    ### clc materialbalance
    # Returns: concentrations
    c_out = xflws.matbal(self, c_in, t, T, i, n_in, prm_dff=False, prm_drc=False)

    ### clc molar flows (n) from concentrations
    n_out = xflws.clc_molarflows(self, c_out, n_in)

    return c_out, n_out
