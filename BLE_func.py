# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 20:15:01 2020

@author: dalgl
"""

from BLE_objects import *


def find_ideal_shield(func, rho_A_max, success_rate):
    # define initial values
    S1 = BLE('HC SW', 'Al 6061-T6').S1
    S2 = BLE('HC SW', 'Al 6061-T6').S2
    t_ob = BLE('HC SW', 'Al 6061-T6').t_ob
    t_b  =BLE('HC SW', 'Al 6061-T6').t_b
    t_w = BLE('HC SW', 'Al 6061-T6').t_wall
    Sigma = BLE('HC SW', 'Al 6061-T6').Sigma
    rho_b = BLE('HC SW', 'Al 6061-T6').rho_b
    Arho_b = BLE('HC SW', 'Al 6061-T6').Arho_b
    # define average inner structure density
    rho_f = rho_b * 0.25
    impact_angle = 0
    
    