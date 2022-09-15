'''
control of plant (Stack?)
'''

import numpy as np
'''
### check temp
if T_min < T_st < T_max:
    temp_ok = True
else:
    temp_ok = False

### cooler/ heater
if plnt_on:
    dT = T_tar - T_st
    if dT >0:
        cooler = True
        heater = False
    elif dT <0:
        heater = True
        cooler = False

###  count idling time
if plnt_on & not plnt_running:
    idling = True
else:
    idling = False
'''

def plnt_ctrl(obj, i_cell):
    '''
    Select operational State of plant

    Parameters
    ----------
    i_cell : float
        Current density of cell     | in A/m²


    Retruns
    -------
    None.
    '''

    obj.av.plnt_on = True # True/ False
    if obj.av.plnt_on & (i_cell>0):
        obj.av.stndby = False # True / False
    elif obj.av.plnt_on & (i_cell<=0):
        obj.av.stndby=True
    else:
        obj.av.stndby=False


    # plnt_on * stndby        -> standby
    # plnt_on & stndby==False -> run
    # plnt_on==False          -> off -> stndby=False
    return

def plnt_swtch_plnt_state(obj, P_in_act, P_in_pre):
    '''
    Detects changes in plant-(operation)state

    P_in_act : Float
        Actual Power input | in kW
    P_in_pre :
        Previous Power input | in kW

    Returns
    -------
    None
    '''


    ### switch to standby
    if (P_in_act <= 0.) & (P_in_pre>0):
        obj.av.swtch_to_stndby = True
    else:
        obj.av.swtch_to_stndby = False

    if True:
        if (P_in_act <= 0.) & (P_in_pre>0):
            obj.av.swtch_to_on = True
        else:
            obj.av.swtch_to_on = False
    return


def plnt_crrnt_ctrl(obj, T):
    '''
    Control (limit) current density according to cell-temperature

    Parameters
    ----------
    T : Float
        Cell temperature | in K


    Returns
    -------
    None
    '''

    if True:
        ### warm up
        if not obj.av.stndby and (T> obj.pcll.temperature_nominal-10):
            obj.av.warmup = False
        elif not obj.av.stndby and (T< obj.pcll.temperature_nominal-10):
            obj.av.warmup = True
        else:
            obj.av.warmup = True
    # obj.av.warmup = False

    if not obj.av.stndby and obj.av.warmup:
        obj.av.i_max_act = obj.pcll.current_density_warmup
        obj.pid_ctrl.reset_I()
    # elif not stndby and not obj.av.warmup:
    else:
        obj.av.i_max_act = obj.pcll.current_density_nominal
    #else:

    #print('Set i_max: ', obj.av.i_max_act)
    return

def plnt_thrm_ctrl(obj, T, u_pid, stndby):
    '''
    Thermal control of plant
    - Heat and coolingwater limits

    Parameters
    ----------
    T : Float
        Cell temperature | in K
    u_pid : Float
        Outputvalue of PID
    stndby : bool
        State of plant (Standby)

    Returns
    -------
    dm_cw : Float
        Massflow of coolant | kg/s
    dQ_heat : Float
        Actual heating power | in kW
    '''
    dm_max = obj.bop.massflow_coolant_max
    dm_min = obj.bop.massflow_coolant_min

    if (u_pid < 0) & (stndby or (T<obj.bop.temperature_standby)):
        dQ_heat = -u_pid*obj.bop.fctr_pid_Pheat
        if dQ_heat > obj.pplnt.power_of_stack_act*1e3*0.002:
            dQ_heat = obj.pplnt.power_of_stack_act*1e3*0.002
        dm_cw = 0
    elif (u_pid > 0) & (stndby):
        dQ_heat = 0 #-u_pid*obj.bop.fctr_pid_Pheat
        dm_cw = 0
    else:
        dQ_heat = 0
        dm_cw = u_pid
    if dm_cw > dm_max:
        dm_cw = dm_max
    elif dm_cw <= dm_min:
        dm_cw = dm_min
    return dm_cw, dQ_heat

# ==============================================================================

### PID
class PID_controller():
    '''
    Class of PID-controller for thermal control of Stack
    '''
    def __init__(self, **kwargs):
        self.enable_P = kwargs.get('enable_P', True)
        self.enable_I = kwargs.get('enable_I', True)
        self.enable_D = kwargs.get('enable_P', True)

        self.Kp = 0
        self.KI = 0
        self.Kd = 0

        self.I_temp_tol = kwargs.get('I_tol', 22)
        self.reset_I_ext = False

        self.SP = kwargs.get('setpoint', 353) # Setpoint

        self.u_min = 0

    def reset(self, ):
        '''
        Reset Controller-Paramters
        (all to 0)
        '''
        self.u = 0
        self.e = 0
        self.e_prv = 0
        self.I_prv = 0
        self.delta_t = 0
        self.P=0
        self.I = 0
        self.D = 0
        return

    def reset_I(self):
        '''
        Reset I-value of PID-control
        (set to 0)
        '''
        self.I_prv = 0
        return

    def set_SP(self,sp):
        '''
        Set setpoint of controller

        Parameters
        ----------
        sp : Float
            Setpoint (Target Temp)
        '''
        self.SP = sp
        return

    def get_tuning_par(self,):
        '''
        Return actual PID-ctrl-Tuning

        Returns
        -------
        Kp : Float
            Linear Gain
        KI : Float
            Integral Gain
        Kd : Float
            Differential Gain
        '''
        return self.Kp, self.KI, self.Kd

    def tune_man(self, kp, ki, kd):
        '''
        Manual tuning of Gain values

        Parameters
        ----------
        kp : Float
            Linear Gain
        ki : Float
            Integral Gain
        kd : Float
            Differential Gain
        '''
        self.Kp = kp
        self.KI = ki
        self.Kd = kd
        if kp <=0:
             self.enable_P= False
        if ki <=0:
             self.enable_I= False
        if kd <=0:
             self.enable_D= False

        return

    def reset_Icmpnt(self,):
        '''
        +++ DUPLICATE ??? +++
        '''
        self.I = 0
        return

    def tune_aut(self):
        ''' +++ ? '''
        return

    def aw_I_6(self, ):
        '''
        anti-windup
        '''
        return self.u

    def clc_components(self, MV, t, t_prev):
        '''
        Calculate actual PID-Components

        Parameters
        ----------
        MV : Float
            Measured Value (of Temperature)
        t : float
            Current time (of calculation step)
        t_prev : Float
            Time-value of previous calculation step
        '''

        self.e_prv = self.e
        self.e       = -(self.SP - MV)                               #SP : Setpoint/ MV: measured value
        self.P         = self.Kp * self.e                                      # Propotional control
        #print('e = ', self.e)
        self.delta_t   = (t - t_prev)                                 #sampling

        self.I_prv = self.I
        if (self.SP-self.I_temp_tol) < MV < (self.SP+self.I_temp_tol) and not self.reset_I_ext:
            self.I         = self.I_prv + (self.KI * self.e * self.delta_t) -0.1*self.aw_I_6()       #Integral control
        else:
            self.I         = 0
        # if self.I <-50:
        #    self.I = 0
        #print('PID-I: ', I)
        self.D        =  self.Kd * (self.e - self.e_prv) / self.delta_t           #Derivative control
        return


    def clc_output(self):
        '''
        Calculate final PID-Output (u-value)

        Returns
        -------
        u : Float
            -> Massflow of coolant
        '''
        self.u = self.P*self.enable_P + self.I * self.enable_I + self.D*self.enable_D
        if np.isnan(self.u):
            self.u = 0
        return self.u
