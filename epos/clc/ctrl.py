'''
control of plant
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
