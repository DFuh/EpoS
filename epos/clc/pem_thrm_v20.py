'''
calculation: thermal behaviour
'''
print(__name__ + ' imported...')

def heatbalance(obj, Tconst=False):
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
        T = obj.pcll.temperature.nominal
        m_c = 0
        m_ely = 0
    #obj.
    return T, m_c, m_ely

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

def clc_coolant_maffflow():
    return


# ----------------------- PID control ---------------------------------------- #
