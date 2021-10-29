'''
AEL
calculation: material balance in cells
'''
import os
import numpy as np
import importlib as impl
from collections import namedtuple
print(__name__ + ' imported...')

import epos.auxf.faux as fx

xflws = fx.dyn_aux_import(__file__, __name__)

def testfun():
    res = xflws.testcalc(2,3)
    print('testfun --> result from testcalc:', res)
    return



# TODO: ini auxvals ???
# TODO: Check >time< of calculating partial pressure

def materialbalance(obj,T, i, m_H2O_in_an, p, c_in, n_in, stf=1,
                                ntd=None, sns=False, m=None): #sns=False):

    t = 0 # dummy variable (for matbal; inc ase of dyn. calc)
    ### pre-clc for matbal
    if sns==False:
        '''
        CHECK -> wKOH
        '''

        clc_matbal_params_Tbased(obj, T, obj.pec.w_KOH)
        clc_matbal_params_ibased(obj, T, i)
    #print('self.d_b: ', self.d_b)
    xflws.matbal_preclc(obj, T,i, )

    ### clc materialbalance
    # Returns: concentrations
    c_out = xflws.matbal(obj, c_in, t, T, i, n_in, prm_dff=False, prm_drc=False)
    #print('c_out: ', c_out)
    ### clc molar flows (n) from concentrations
    n_out, pp_out = xflws.clc_molarflows(obj, c_out, n_in)

    #print('---n_out: ', n_out)
    x_H2inO2 = n_out[0] / (n_out[0] + n_out[2])
    #print('calc materialbalnce for ', self.name, 'i= ', i, 'x_= ', x_H2inO2)


    if (not sns) and (ntd!=None):


        ntd.n_H2_an[m], ntd.n_H2_ca[m] = n_out[0]*stf, n_out[1]*stf
        ntd.n_O2_an[m], ntd.n_O2_ca[m] = n_out[2]*stf, n_out[3]*stf
        ntd.c_H2_an[m], ntd.c_H2_ca[m], ntd.c_O2_an[m], ntd.c_O2_ca[m] = c_out
        ntd.x_H2_ca[m], ntd.x_H2_an[m], ntd.x_O2_ca[m], ntd.x_O2_an[m] = 0, x_H2inO2, 0,0
        ntd.pp_H2_an[m], ntd.pp_H2_ca[m], ntd.pp_O2_an[m], ntd.pp_O2_ca[m] = pp_out
        return
    else:
        FLWS = namedtuple('FLWS', '''n_H2_out_ca n_H2_out_an n_O2_out_ca n_O2_out_an
                                    c_H2_out_ca c_H2_out_an c_O2_out_ca c_O2_out_an
                                    x_H2_out_ca x_H2_out_an x_O2_out_ca x_O2_out_an
                                    pp_H2_mem_ca pp_H2_mem_an pp_O2_mem_ca pp_O2_mem_an
                                    pp_H2O_mem_an''')
        n_H2_ch_an, n_H2_ch_ca, n_O2_ch_an, n_O2_ch_ca = n_out
        c_H2_ch_an, c_H2_ch_ca, c_O2_ch_an, c_O2_ch_ca = c_out
        x_H2_ch_ca, x_H2_ch_an, x_O2_ch_ca, x_O2_ch_an = 0, x_H2inO2, 0,0
        pp_H2_mem_an, pp_H2_mem_ca, pp_O2_mem_an, pp_O2_mem_ca = pp_out

        flws_out = FLWS(n_H2_ch_ca, n_H2_ch_an, n_O2_ch_ca, n_O2_ch_an,
                        c_H2_ch_ca, c_H2_ch_an, c_O2_ch_ca, c_O2_ch_an,
                        x_H2_ch_ca, x_H2_ch_an, x_O2_ch_ca, x_O2_ch_an,
                        pp_H2_mem_ca, pp_H2_mem_an, pp_O2_mem_ca, pp_O2_mem_an,
                        obj.av.pp_H2O ) # pp_H2O -> pp_sat_KOH


        return flws_out #c_out, n_out, x_H2inO2


def clc_matbal_params_Tbased(obj,T, w_KOH):
    ### clc T-params
    obj.av.D_ik       = xflws.clc_D_ik(obj,T, w_KOH)

    obj.av.Vhcell     = xflws.clc_Vhcell(obj,T, w_KOH,
                        fctr=getattr(obj,'vhcfctr',1))

    obj.av.gamma      = xflws.clc_gamma(obj,T, w_KOH)

    obj.av.rho_H2     = xflws.clc_rho_H2(obj,T, w_KOH)

    obj.av.rho_KOH    = xflws.clc_rho_KOH(obj,T, w_KOH)

    obj.av.beta       = xflws.clc_beta(obj, T, fctr=getattr(obj, 'betafctr',1))

    obj.av.S_ik       = xflws.clc_S_ik(obj, T, w_KOH)

    obj.av.pp_H2O     = xflws.clc_pp_H2O(obj, obj.pec, T)#, w_KOH)

    return

def clc_matbal_params_ibased(obj,T,i):
    ### clc i-params

        ### clc flow of electrolyte (VL)
    clc_VL(obj, i=i,dynVely=getattr(obj,'dyn_ely', False),
                balVely=getattr(obj, 'bal_ely', True))

    obj.av.epsilon    = xflws.clc_epsilon(obj, T, i, fctr=getattr(obj, 'epsfctr',1 ))

    obj.av.d_b        = xflws.clc_d_b(obj, T, i, )
    #print('obj.d_b: ', obj.d_b)
    obj.av.k_Li       = xflws.clc_k_Li(obj, T, i, d_b=obj.av.d_b,
                                        D_ik=obj.av.D_ik, epsilon=obj.av.epsilon,
                                        fctr=getattr(obj, 'klifctr',1))

    obj.av.A_GL       = xflws.clc_A_GL(obj, T, i)

    obj.av.f_G        = xflws.clc_f_G(obj, T, i, fctr=getattr(obj, 'fGfctr',1 ))
    return


def clc_VL(self, i=None, dynVely=False, balVely=True):
    '''
    calc electrolytre volumetric flowrate
    orig: parameter_sensibility_bsc_matbal_v1002.py
    '''
    def clc_dynVely(i, dV0, dV_min):
        #i_max   = 0.5#*1e4         #m_par.i_max

        #Vely_min = 0.2 * V0 #m_par.Vely_min
        #print('i / self.i_max: ', (i / self.i_max))
        #print(f'i= {i}, i_max={self.i_max}')
        #V_is    = V0 * i / self.i_max
        V_is    = dV0 * i / self.pcll.current_density_max

        # if V_is < self.Vely_min:
        if V_is < dV_min:
            # V_is = self.Vely_min
            V_is = dV_min
        return V_is

    #print('dynVely: ', dynVely)
    # dV0_ely = (self.bop.volumetricflow_ely_nominal/
    #             self.pplnt.number_of_cells_in_stack_act)
    dV0_ely = self.bop.volumetricflow_ely_nominal_matbal


    if not balVely:
        self.av.VL_an = dV0_ely * self.fac_Vely_an
        self.av.VL_ca = dV0_ely * self.fac_Vely_ca
    else:
        self.av.VL_an           = dV0_ely # liquid (electrolyte) flowrate; anode // in m³/s | Haug 2017, Tab2: 330mL/min
        self.av.VL_ca           = dV0_ely # liquid (electrolyte) flowrate; cathode // in m³/s

    if i and dynVely:
        print(' flws_aux l. 158 ---> dynEly active !')
        dV_ely_min = dV0_ely * self.bop.fctr_volumetricflow_ely_min
        self.av.VL_an = clc_dynVely(i, self.av.VL_an, dV_ely_min)
        self.av.VL_ca = clc_dynVely(i, self.av.VL_ca, dV_ely_min)
    self.av.VL = (self.av.VL_an + self.av.VL_ca) / 2 # for clc c_mix
    #print(f'VL = {self.VL}, VL_an= {self.VL_an}, VL_ca= {self.VL_ca}') #'??? line 155 in clcVL')
    return

def partial_pressure(obj, pec, T, p):
    '''
    simple calculation (!)
    partial pressure of product gases dependend on water-vapor-pressure
    '''
    #T in K
    #p_ca, p_an = p_in
    pp_H2O = xflws.clc_pp_H2O(obj, pec, T, )
    pp_H2_ca = p.cathode - pp_H2O

    pp_O2_an = p.anode - pp_H2O

    return pp_H2_ca, pp_O2_an, pp_H2O #av.plr_ppr[0][1:3]
