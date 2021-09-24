'''
control of plant (Stack?)
'''

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

def plnt_thrm_ctrl(obj, T, u_pid, stndby):
    dm_max = obj.bop.massflow_coolant_max
    dm_min = obj.bop.massflow_coolant_min

    if (u_pid < 0) & (stndby or (T<333)):
        dQ_heat = -u_pid*obj.bop.fctr_pid_Pheat
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
    def __init__(self, **kwargs):
        self.enable_P = kwargs.get('enable_P', True)
        self.enable_I = kwargs.get('enable_I', True)
        self.enable_D = kwargs.get('enable_P', True)

        self.Kp = 0
        self.KI = 0
        self.Kd = 0

        self.SP = kwargs.get('setpoint', 353) # Setpoint

        self.u_min = 0

    def reset(self, ):
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
        self.I_prv = 0
        return

    def set_SP(self,sp):
        self.SP = sp
        return

    def get_tuning_par(self,):
        return self.P, self.I, self.D

    def tune_man(self, kp, ki, kd):
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


    def tune_aut(self):

        return

    def aw_I_6(self, ):
        '''
        anti-windup
        '''
        return self.u

    def clc_components(self, MV, t, t_prev):

        self.e_prv = self.e
        self.e       = -(self.SP - MV)                               #SP : Setpoint/ MV: measured value
        self.P         = self.Kp * self.e                                      # Propotional control
        #print('e = ', self.e)
        self.delta_t   = (t - t_prev)                                 #sampling

        self.I_prv = self.I
        self.I         = self.I_prv + (self.KI * self.e * self.delta_t) -0.1*self.aw_I_6()       #Integral control

        #print('PID-I: ', I)
        self.D        =  self.Kd * (self.e - self.e_prv) / self.delta_t           #Derivative control
        return


    def clc_output(self):
        self.u = self.P*self.enable_P + self.I * self.enable_I + self.D*self.enable_D
        return self.u
