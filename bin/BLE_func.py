# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 20:15:01 2020

@author: dalgl
"""

# import BLE_objects as bo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import master_output_reader as mor
import directories as dirs
import math, cProfile, pstats


def read_master_data(text_file):
    f = open(dirs.path+text_file, 'r')
    content = f.readlines()
    f.close()
    # Define cross-sectional area of satellite
    A = 5       # m^2
    # Define the time on orbit for the satellite
    t = 9       # years
    # determine sum of probabilities for this sample in order to normalise probability values
    sump = 0
    # create dictionaries for flux and probability data
    flux_dat = {}
    prob_dat = {}
    n = 0
    for i in content:
        try:
            lin = i.split()
            F = float(lin[19])#/float(lin[0])             # Takes the flux value 
            p = 1 - math.exp((-1)*F*A*t)    # Calculates the probability of this impact
            prob_dat[n] = float(lin[0]), float(lin[1])*1E2, p # append values to sheets to arrays and consert diameter to cm
            flux_dat[n] = float(lin[0]), float(lin[1])*1E2, F
            n+=1
            sump += p
        except:
            pass
    # normalise probabilities
    prob_dat = [(*i, i[2]/sump) for i in list(prob_dat.values())]

    return np.array(prob_dat), np.array(list(flux_dat.values()))


def plot_heatmap(data_file):
    data = read_master_data(data_file)[0]
    # initialise 1D lists for vel(x), di(y) and probability()
    N = data.shape[0]
    x_start = data[0, 0]
    y_start = data[0,1]
    z_start = data[0,3]
    X = [x_start]
    Y = [y_start]
    Z = [[] for i in range(N)]
    Z[0].append(z_start)
    a = 0
    b = 0
    # c = 0
    count = 1
    
    # create 1D lists of X and Y and a square matrix of probability values 
    while count < N:
        x = data[count,0]
        y = data[count,1]
        z = data[count,3]
        if x != data[count-1, 0]:
            X.append(x)
            a+=1
        elif y not in Y:
            Y.append(y)
            b+=1
        Z[a].append(z)
        count+=1
    # remove empty lists in Z
    for j in range(N):
        i = j+1
        if not Z[N-i]:
            Z.pop(N-i)
    # create empty square matrices of X and Y for contourf plot
    sq_vel = [[] for i in range(a+1)]
    sq_di = [[] for i in range(b+1)]

    #  fill matrices with values
    for i in range(a+1):
        u = 0
        while u<a+1:
            sq_vel[i].append(X[i])
            sq_di[i].append(Y[i])
            u+=1

    sq_vel = np.array(sq_vel).transpose()
    sq_di = np.array(sq_di)
    Z=np.array(Z).transpose()

    # plot heatmap
    plt.contourf(sq_vel, sq_di, Z, 100, cmap = 'RdPu')
    plt.colorbar()
    # fig, ax = plt.subplots()
    # im = ax.imshow(data)
    # return Z, sq_vel, sq_di


def plot_heat_and_limit(data_file, BLE_conditions):
    if isinstance(BLE_conditions, tuple):
        for i in BLE_conditions:
            i.graph()
    else:
        BLE_conditions.graph()
    plot_heatmap(data_file)


def save_shield_data(BLE_object, ideal_output_arr):
    S1 = ideal_output_arr[0]
    tob = ideal_output_arr[1]
    tb = ideal_output_arr[2]
    succ = ideal_output_arr[3]
    Arho_b = ideal_output_arr[4]
    rho_b = BLE_object.rho_b
    Sigma = BLE_object.Sigma
    S2 = BLE_object.S2
    t_wall = BLE_object.t_wall
    # create sheetname
    #BLE_object.lbl
    sht_me = BLE_object.lbl+str(round(succ,2))
    # create index and columns list for pandas dataframe
    idx = ['Distance between first and second bumper/back wall',
           'Thickness of front bumper',
           'Volumetric density of front bumper',
           'Areal density of shield (not inc. wall)',
           'rear wall yield stress',
           'Thickness of second bumper',
           'Distance between second bumper and back plate',
           'Thickness of back plate/wall']
    cols = ['Denotion', 'Value', 'Units']
    # create pandas dataframe
    df = pd.DataFrame(np.array([['S1', S1, 'cm'],
                                ['t_ob', tob, 'cm'],
                                ['rho_b', rho_b, 'g.cm^-3'],
                                ['Arho_b', Arho_b, 'g.cm^-2'],
                                ['sigma', Sigma, 'ksi'],
                                ['t_b', tb, 'cm'],
                                ['S2', S2, 'cm'],
                                ['t_wall', t_wall, 'cm']]),
                      index = idx, columns=cols)
    with pd.ExcelWriter(dirs.in_dir+'ShieldProperties.xlsx', engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name = sht_me)


def time_stats(func,*args):
    pr = cProfile.Profile()
    pr.runcall(func, *args)
    ps = pstats.Stats(pr).sort_stats('cumulative')
    ps.print_stats(10)
    