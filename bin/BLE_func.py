# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 20:15:01 2020

@author: dalgl
"""

import BLE_objects as bo
import numpy as np
import matplotlib.pyplot as plt

# directory to model inputs
ipt_dir = 'C:/Users/dalgl/OneDrive/Documents/4th year/GDP/Space-Shield-GDP/inputs'


def read_master_data(text_file):
    f = open(text_file, 'r')
    content = f.readlines()
    f.close()
    data = {}
    n = 0
    for i in content:
        i.strip('\n').split('\t')
        data[n] = float(i[0]), float(i[1]), float(i[19])
        n+=1
    return np.array(list(data.values()))


def plot_heatmap(data_file):
    data = read_master_data(data_file)
    plt.plot(data[:, 0], data[:,1], 'c')
    plt.axis([0, BLE.maxvel+1, 0, data[0, 1]+0.1])
    plt.xlabel('Velocity (km.s-1)')
    plt.ylabel('Crit Diameter (cm)')

