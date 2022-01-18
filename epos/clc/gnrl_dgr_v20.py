#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:45:01 2021

@author: dafu_res
"""

def voltage_increase(obj, ):

    ### absolute voltage increase
    dU_dgr_abs = obj.av.t_op/3600 * obj.pec.fctr_vlr # h * V/h

    ### incremental voltage increase
    dU_dgr_incr = obj.av.t_diff/3600 * obj.pec.fctr_vlr # h * V/h

    return dU_dgr_incr, dU_dgr_abs
