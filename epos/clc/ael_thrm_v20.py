'''
calculation: thermal behaviour
'''
print(__name__ + ' imported...')

def heatbalance(obj, T_in, m_ely_in, m_c_in, Tconst=False):
    '''
    calc. heat balance of Stack(s) and water cycle

    Parameters
    ----------
    obj: simu object
    T_in: float
        Stack temperature in K
    m_c_in: float
        massflowrate of coolant (in hex) in kg/s
    m_ely_in: float
        massflowrate of water/electrolyte entering stack


    returns
    -------
    T_out: float
        Actual temperature of stack in K
    m_c_out: float
        Adapted coolant massflowrate in kg/S


    '''
    if not Tconst:
        pass
    else:
        T_out = 353#T_in
        obj.clc_m.flws.xflws.clc_flws_auxpars(obj, T_out)#ntd.T_st[m]) #???
        m_ely_out = (obj.bop.volumetricflow_ely_nominal * obj.av.rho_ely
                        * obj.pplnt.number_of_stacks_act)
        m_c_out = m_ely_out
        P_heat = 0

    return T_out, m_ely_out, m_c_out, P_heat
