# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 17:48:48 2021

@author: dalgl
"""
from BLE_objects import *
from BLE_func import *
import time
import matplotlib.pyplot as plt
import numpy as np


def op1(master_file, ballistic_model, shieldsheet, projsheet, success, succ_meth=1):
    start_time = time.time()
    base_shield = ballistic_model(shieldsheet, projsheet)
    print(f'Base Shield created: {time.time()-start_time:.2f}s\nFinding Ideal...\n')
    if succ_meth is 1:
        id_shield = base_shield.ideal(success, succ_it = succ_meth, file = master_file)
    else:
        id_shield = base_shield.ideal(success, succ_it = succ_meth)
    print(f'Ideal shield found: {time.time()-start_time:.2f}s\n')
    save_shield_data(base_shield, id_shield)
    print(f'Ideal shield saved: {time.time()-start_time:.2f}s\n')
    new_shield = ballistic_model(base_shield.lbl+str(round(id_shield[3],2)), projsheet)
    print(f'Ideal shield created: {time.time()-start_time:.2f}s\n')
    new_shield.graph()
    base_shield.graph()
    print(f'Graphs plotted: {time.time()-start_time:.2f}s\n Operation Complete\n')
    return base_shield, new_shield
    


def op2(ballistic_model, shieldsheet, projsheet, shield_attr, val_rng, num_steps = 10):
    base_shield = ballistic_model(shieldsheet, projsheet)
    stp = (val_rng[1]-val_rng[0])/(num_steps-1)
    vals = np.arange(start = val_rng[0], stop = val_rng[1]+stp, step=stp)
    if shield_attr in base_shield.__dict__:
        for i in vals:
            setattr(base_shield, shield_attr, i)
            setattr(base_shield, 'lbl', str(round(i,2)))
            setattr(base_shield, 'col', (i/val_rng[1]))
            base_shield.graph()
