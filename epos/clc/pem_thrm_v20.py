'''
calculation: thermal behaviour
'''
print(__name__ + ' imported...')

def heatbalance(obj, T_st_in, m_c_in, m_ely_in, ntd=None, Tconst=False):
    '''
    mainfunction for thermal calc.
    -> clc. Stack Temperature
    --> stack water inflow conditioning:
        -> clc. water reservoir Temp. (?)
        -> clc. coolant massflow
        -> power of preheater

    '''
    if not Tconst:
        pass
    else:
        T_out = obj.pcll.temperature_nominal
        obj.clc_m.flws.xflws.clc_flws_auxpars(obj, T_out)#ntd.T_st[m]) #???
        m_ely_out = (obj.bop.volumetricflow_ely_nominal * obj.av.rho_ely
                        * obj.pplnt.number_of_stacks_act) #V0: on Stack level
        m_c_out = m_ely_out # Check level!
        P_heat = 0 # Check level!
    #obj.
    return T_out, m_ely_out, m_c_out, P_heat # Output on plant level

# ----------------------- Temperature Stack ---------------------------------- #
def clc_temperature_stack():
    # clc T_St

    # solve ode

    return


def dydt():
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
