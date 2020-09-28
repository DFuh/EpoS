'''
tools for unit conversion
'''

def convert_unit(name = '', num=None):
    '''
    return factor in order to convert to basic unit

    basic units:
    power: kW   | kilowatts
    energy: kWh | kilowatthours
    time: s     | seconds
    length: m   | meters

    '''
    def u_sw(name):
        u_switch = {
                    "MW":1e3
                    "kW": 1,
                    "W": 1e-3,

                    "kWh":1,
                    "Wh": 1e-3,
                    "J": ,
                    "kJ" ,

                    }

        return u_switch.getattr(name, str(convert_unit.__name__)+': unit not found ...')

    fac = float(u_sw(name))
    if not num:
        return fac
    else:
        ret = num * fac
        return ret
