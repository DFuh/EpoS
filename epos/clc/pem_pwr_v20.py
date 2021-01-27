'''
calculation: power consumption and control
    -> optimal power/ operation point
        -> suitable power gradient
        ->
    -> power consumption BoP
        -> power consumption pumps (feedwater, (coolant?, compressor?)
        -> power characteristics (efficiencies): rectifier, el-motor, fan?
        -> power consumption gas treatment (dryer)

'''
import scipy.optimize as scpt

print(__name__ + ' imported...')





### Balance of plant
def BoP_power(n_H2, m_c_in, Vely, p_in, ):
    ''' power consumption of plant/ periphery'''

    ### gas treatment (drying)
    P_gt    = pow_gas_trt(n_H2, ) # power input gas treatment         // in W

    m_c0    = pv.m_cmax
    P_pmp   = pow_pump_m(m_c_in, m_c0, p_A, )    # power input feed/cooling pump     // in W
    #eff_pe  = ip_rect()     # rectifier efficiency
    #P_pe    = P_st * (1-eff_pe) # power input/loss power electronics

    P_bop = P_gt + P_pmp + P_ely # // in W

    return P_bop
