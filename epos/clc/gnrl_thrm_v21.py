'''
calculation: thermal behaviour
'''
print(__name__ + ' imported...')

import numpy as np
from scipy.integrate import odeint

from epos.clc import ctrl

def heatbalance(obj, T_st_in, m_ely_in, m_c_in, u_cell, i_cell,
                n_H2_ca, n_O2_an, n_H2O_cns,
                t_arr, ntd=None, Tconst=False,
                par_out=False):
    '''
    mainfunction for thermal calc.
    -> clc. Stack Temperature
    --> stack water inflow conditioning:
        -> clc. water reservoir Temp. (?)
        -> clc. coolant massflow
        -> power of preheater

    '''
    if not Tconst:
        #pass

        ### ctrl values

        ctrl.plnt_ctrl(obj, i_cell)
        # print('stndby (heatbal): ', obj.av.stndby)
        ### PID control
        # print(f'---> thrm || u_cell= {u_cell} | i_cell= {i_cell}')
        # print('--> PID: ', obj.pid_ctrl.get_tuning_par(), 'flw_mx= ', obj.mssflw_max)
        if obj.av.swtch_to_stndby:
            obj.pid_ctrl.reset_Icmpnt()
        if obj.av.stndby:
            obj.pid_ctrl.set_SP(obj.bop.temperature_standby)
        else:
            obj.pid_ctrl.set_SP(obj.pcll.temperature_nominal)
        obj.pid_ctrl.clc_components(T_st_in, t_arr[1], t_arr[0])
        u_pid = obj.pid_ctrl.clc_output()
        # print('u_pid= ', u_pid)
        m_c_out, P_heat = ctrl.plnt_thrm_ctrl(obj, T_st_in, u_pid, obj.av.stndby)
        #td = np.array(t_arr)
        # print('m_c_out=', m_c_out)
        # P_heat = 0

        m_ely_out = obj.bop.volumetricflow_ely_nominal * obj.pplnt.number_of_cells_in_stack_act
        # m_ely_out = 10*2*0.18*0.01801528*obj.pplnt.number_of_cells_in_stack_act
        # m_ely_out = m_ely_in # 0.30*obj.av.stckfctr #?

        if np.isnan(n_H2_ca) or (n_H2_ca<0):
            n_H2_ca = 0
        if np.isnan(n_O2_an) or (n_O2_an < 0):
            n_O2_an = 0


        T_out, params_thrm = clc_temperature_stack(obj, T_st_in, m_c_out, P_heat,
                                        u_cell, i_cell,
                                        n_H2_ca, n_O2_an, n_H2O_cns,
                                        t_arr, ntd=ntd,
                                        par_out=par_out)
        # print('T_Stack_in = {0} | T_Stack_out = {1}'.format(T_st_in, T_out))
        # print('P_heat: ', P_heat/1e3)
    else:
        T_out = obj.pcll.temperature_nominal
        obj.clc_m.flws.xflws.clc_flws_auxpars(obj, T_out)#ntd.T_st[m]) #???
        m_ely_out = (obj.bop.volumetricflow_ely_nominal * obj.av.rho_ely
                        * obj.pplnt.number_of_stacks_act) #V0: on Stack level
        m_c_out = m_ely_out # Check level!
        P_heat = 0 # Check level!
    #obj.
    return (T_out,
            m_ely_out,
             m_c_out,
              P_heat/1e3,
              params_thrm) # Output on Stack level

# ----------------------- Temperature Stack ---------------------------------- #
def clc_temperature_stack(obj, T_act, dm_cw, dQ_heat,
                            u_cell, i_cell, n_H2_ca, n_O2_an, n_H2O_cns,
                            t_arr,
                            ntd=None, dyn=False,
                            par_out=False):

    # dyn=True
    # clc T_St
    n_H2O_resid_out = 0


    eta_e = obj.pec.u_tn/u_cell if u_cell >0 else 0
    cpm_O2 = f_cpm(T_act, *obj.bop.args_cpm_O2) # Espinosa-Lopez Tab 2 // in J/molK
    cpm_H2 = f_cpm(T_act, *obj.bop.args_cpm_H2)
    if obj.scn_setup:
        kA_fctr = 1 #kA_fun(dm_cw, m0=obj.bop.massflow_coolant_max)
    else:
        kA_fctr = kA_fun(dm_cw, m0=obj.bop.massflow_coolant_max)
        kA_fctr = kA_fctr if kA_fctr<1 else 1
    U_HAx = obj.pplnt.UA_hx0_st * kA_fctr #kA_fun(dm_cw, m0=obj.bop.massflow_coolant_max)

    # print(f'kA = {kA_fctr} ; UHAx={U_HAx}')

    n_c     = obj.pplnt.number_of_cells_in_stack_act
    n_st    = obj.pplnt.number_of_stacks_act
    A       = obj.pcll.active_cell_area
    C_cw    = dm_cw * obj.bop.cp_coolant #// in ?
    # exp_f = C_cw * (1 - np.exp(-U_HAx / C_cw)) if C_cw >0 else 0
    exp_f   = (1 - (np.exp(-U_HAx / C_cw) if C_cw >0 else 0))
    # print('exp-f: ', exp_f)
    Tcorr = 0. #273.15
    # print('n_i in heatbal: n_H2_ca, n_O2_an, n_H2O_cns, n_H2O_resid_out', n_H2_ca, n_O2_an, n_H2O_cns, n_H2O_resid_out)
    par_b = (obj.bop.temperature_ambient-Tcorr, obj.bop.temperature_coolant-Tcorr, C_cw,
                n_H2_ca, n_O2_an, n_H2O_cns,
                n_H2O_resid_out,
                u_cell, i_cell, eta_e, n_c, n_st, A,
                obj.pplnt.heat_capacity_st,
                obj.pplnt.thermal_resistance_st, obj.pplnt.UA_hx0_st, dQ_heat,
                cpm_H2 ,cpm_O2, obj.bop.cpm_H2O, exp_f)

    par_a = (C_cw, n_H2_ca, n_O2_an, n_H2O_cns,
            n_H2O_resid_out,  n_c, n_st, A,
            obj.pplnt.heat_capacity_st, obj.pplnt.thermal_resistance_st,
            obj.pplnt.UA_hx0_st, cpm_H2 ,cpm_O2, obj.bop.cpm_H2O, exp_f)

    # print('Q_cool=', C_cw*(T_act-obj.bop.temperature_coolant))
    #t_arr = [0, (t_act-t_prv)]

    ### clc Stack Temperature
    if dyn:
        T_ = odeint(dydt, T_act-Tcorr, t_arr, args=(par_a, par_b))[-1] +Tcorr
    else:
        T_ = dydt_slvd(T_act-Tcorr, np.diff(t_arr), par_a, par_b) +Tcorr

    if par_out:
        return T_, (par_a, par_b)
    else:
        return T_, None


def f_cpm(T, a,b,c,d):
    return a + b*T + c*T**2 + d*T**3

def dydt(T, t, par_a, par_b):
    '''
    ODE based on Ulleberg [] (and Espinosa-Lopez [])
    '''
    return eq_b(*par_b) - eq_a(*par_a)*T


def dydt_slvd(T_ini, t, par_a, par_b):
    '''
    analytically solved ODE based on Ulleberg [] (and Espinosa-Lopez [])
    '''
    res_a = eq_a(*par_a)
    res_b = eq_b(*par_b)
    return (T_ini - (res_b/res_a)) * np.exp(-res_a * t) + (res_b/ res_a)


def eq_b(T_a, T_cwi, C_cw, n_H2, n_O2, n_H2O_cns_in, n_H2O_resid_out,
         U, i_cell, eta_e, n_c, n_st, A_c, C_t, R_t, U_HAx, dQ_heat,
         cp_mH2 ,cp_mO2, cp_mH2O, exp_f):
    # print('Q_gen = ', (n_c*A_c * U * i_cell * (1-eta_e)))

    # print(locals())
    # print('heat gain flws: ', (n_O2 * cp_mO2 + n_H2 * cp_mH2-n_H2O_cns_in * cp_mH2O)*T_a)
    # print('heat gain Ui: ', n_c*A_c * U * i_cell * (1-eta_e))
    # print( 'heat gain Ccw: ', C_cw * T_cwi*exp_f)
    return 1/C_t * ( (n_c*A_c * U * i_cell * (1-eta_e)) + (C_cw * T_cwi*exp_f )
                    + (T_a/R_t) #+
                     + (n_O2 * cp_mO2 + n_H2 * cp_mH2)                 # CHECK sign of molar flows !!
                    #- (n_H2O_cns_in * cp_mH2O)
                    *T_a
                    + dQ_heat)



def eq_a(C_cw, n_H2, n_O2, n_H2O_cns_in, n_H2O_resid_out,
         n_c, n_st, A_c, C_t, R_t, U_HAx, cp_mH2, cp_mO2, cp_mH2O, exp_f):
    # print(locals())
    # print('heat loss exp-f: ',exp_f*C_cw)
    # print('heat loss flws: ', (n_O2 * cp_mO2 + n_H2 * cp_mH2-n_H2O_cns_in * cp_mH2O))
    return 1/C_t * ((1/R_t + exp_f*C_cw )
                    + n_O2 * cp_mO2 + n_H2 * cp_mH2)          # CHECK sign of molar flows !!
                    # -n_H2O_cns_in * cp_mH2O)

def kA_fun(m_act, m0=1, ):
    val_x0 = 0.01
    pars = (1.43144695, 0.63439993, 2.47031518, 0.12110944)

    if m_act <= 1e-9:
        kA_fctr = val_x0
    else:

        kA_fctr = log_growth(m_act/m0, *pars)
    return kA_fctr

def log_growth(x,a,b,c,d):
    return a +b*np.log(d+x/c)



def scl_thrm():
    ''' Returns scaling factor for thrml clc.'''
    # -> either stack level
    # -> or plant level
    return

# ----------------------- Water inflow conditioning -------------------------- #

def water_inflow_conditioning():

    # clc temp of water reservoir

    # clc heat exchanger

    # clc coolant massflow // ventilator power

    return

# ----------------------- Coolant massflow ----------------------------------- #

def clc_coolant_massflow():
    return


# ----------------------- PID control ---------------------------------------- #
