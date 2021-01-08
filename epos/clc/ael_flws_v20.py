'''
AEL
calculation: material balance in cells
'''
import os
import numpy as np
import importlib as impl
print(__name__ + ' imported...')

'''
# -> import aux-funct.:
ref = 'Projects.EpoS.epos.clc.aux' #'.EpoS.epos.clc'
prfx, sffx = os.path.basename(__file__).split('_v')
s = '.'+prfx+'_aux_v' + sffx.split('.py')[0]
print('s: ', s)
print('cwd:', os.getcwd())
xflws = impl.import_module(s, ref)
'''
def testfun():
    res = xflws.testcalc(2,3)
    print('testfun --> result from testcalc:', res)
    return

def materialbalance(obj, ): #sns=False):

    ### clc T-params
    obj.D_ik       = xflws.clc_D_ik(obj.T, obj.w_KOH) #* (1 + (sns * fctr -1))

    obj.Vhcell     = xflws.clc_Vhcell(obj.T, obj.w_KOH)

    obj.gamma      = xflws.clc_gamma(obj.T, obj.w_KOH)

    obj.rho_H2     = xflws.clc_rho_H2(obj.T, obj.w_KOH)

    obj.rho_KOH    = xflws.clc_rho_KOH(obj.T, obj.w_KOH)

    obj.beta       = xflws.clc_beta(fctr=obj.betafctr)

    obj.S_ik       = xflws.clc_S_ik(obj.T, obj.w_KOH)

    obj.pp_H2O     = xflws.clc_pp_H2O(obj.T, obj.w_KOH)


    ### clc i-params

        ### clc flow of electrolyte (VL)
    obj.clc_VL(i=i,dynVely=obj.dyn_ely, balVely=obj.bal_ely)

    obj.epsilon    = par.clc_epsilon(obj, i)

    obj.d_b        = par.clc_d_b(obj, i, obj.i_crit_db, obj.gamma, obj.beta,
                                    obj.rho_H2, obj.rho_KOH)

    obj.k_Li       = par.clc_kLi(obj, obj.T, i, d_b=obj.d_b,
                                        D_ik=obj.D_ik, epsilon=obj.epsilon,
                                        fctr=self.klifctr)

    obj.A_GL       = par.clc_A_GL(obj.Vhcell, obj.epsilon, obj.gamma, obj.d_b)

    obj.f_G        = par.clc_f_G(obj, i, fctr=obj.fGfctr)



    ### pre-clc for matbal
    matbal_preclc(self, T,i, )

    ### clc materialbalance
    # Returns: concentrations
    matbal(self,y_in, t, T, i, n_in, prm_dff=False, prm_drc=False)

    ### clc molar flows (n) from concentrations
    clc_molarflows(self, c_in, n_in)

    return
