# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 11:11:15 2020

@author: dalgl
"""
import numpy as np
import pandas as pd
import pylab
from itertools import product
import time


class BLE:
    def __init__(self, shieldsheet, projsheet):
        """
        This function defines the variables from a given shield sheet and
        projectile sheet
        """
        varis = pd.DataFrame(pd.read_excel(
            r'C:\Users\dalgl\OneDrive\Documents\4th year\GDP\Space-Shield-GDP\ShieldProperties.xlsx',
            sheet_name=shieldsheet, index_col=0,
            header=0, dtype={'denotion': str,
                             'Value': float}))
        # import data from excel sheet and convert to numpy array
        varip = pd.DataFrame(pd.read_excel(
            r'C:\Users\dalgl\OneDrive\Documents\4th year\GDP\Space-Shield-GDP\ProjectileProperties.xlsx',
            sheet_name=projsheet, index_col=0, header=0,
            dtype={'Value': float}))

        # convert pandas dataframe to numpy array
        varo = varis.values  
        vara = varip.values

        # Assign variables for shield
        self.S1 = varo[0, 1]
        self.t_ob = varo[1, 1]
        self.rho_b = varo[2, 1]
        self.Arho_b = varo[3, 1]
        self.Sigma = varo[4, 1]
        self.t_b = varo[5, 1]
        self.S2 = varo[6, 1]
        self.t_wall = varo[7, 1]
        # assign variablesfor projectile
        self.di = vara[0, 0]
        self.rho_p = vara[1, 0]
        # assign input variables 
        self.impact_angle = 0
        self.vel = 1
        # create dcrit arrays
        self.rngdcrit = []
        # define iteration inputs
        self.Arho_max = 3 # g
        self.succ_rate = 95 # %
        # define label inputs
        self.lbl = 0
        self.maxvel = 12
        self.minvel = 0.5
        # set shatter region limits
        # ideally shock analysis done to work this out
        self.v_shat = 4.2 # kms^-1
        self.v_vap = 8.4 # kms^-1


    def drange(self, vel):
        self.vel = vel
        # for ballistic region - Vn <= 4.2kms^-1
        if self.vel <= self.v_shat:
            dc = self.crit_d_bal(self.vel)
        # for transition region
        elif self.v_shat < vel < self.v_vap:
            v1 = self.v_shat
            v2 = self.v_vap
            dbal = self.crit_d_bal(v1)
            dvap = self.crit_d_vap(v2)

            dc = dbal + ((dvap - dbal)/(v2 - v1)) * (vel - v1)
        # for vapourised region
        else:
            dc = self.crit_d_vap(self.vel)
        return dc


    def graph(self):
        vd = np.array(self.arr())
        pylab.plot(vd[:, 0], vd[:,1], 'c')
        pylab.axis([0, self.maxvel+1, 0, vd[0, 1]+0.1])
        pylab.xlabel('Velocity (km.s-1)')
        pylab.ylabel('Crit Diameter (cm)')


    def arr(self):
        """
        This function returns an array of normal projectile velocities and
        their corresponding dcrit for the SRL BLE
        """
        # determine critical diameter for each velocity with given shield
        # and projectile properties
        d = {}
        for i in np.linspace(self.minvel, self.maxvel, 100):
            d[i] = float(self.drange(i))
        self.rngdcrit = list(d.items())
        return self.rngdcrit


    def ideal(self, Arho_max, succ_rate):
        start_time = time.time()
        # define % of non-perforation over velocity range
        self.succ_rate = succ_rate
        self.Arho_max = Arho_max
        # define output list containing successful variables
        a = {}
        # define range of shield values to vary and the ranges of variation
        S1rng = np.arange(start = 0.1, stop = 5, step = 0.1)
        t_obrng = np.arange(start = 0.1, stop = 1, step = 0.05)
        t_brng = np.arange(start = 0.1, stop = 1, step = 0.05)
        # define ratio between bumper density and inner structure density
        rho_ratio = 0.5
        # iterate for all combinations
        # prog_count = 0
        itr_count = 0
        for i, j, k in product(S1rng, t_obrng, t_brng):
            self.S1 = i
            self.t_ob = j
            self.t_b = k
            # filter out area ratios greater than Arho_max
            A_rho = self.rho_b * (j + k + rho_ratio * i)
            if A_rho <= self.Arho_max:
                g = self.arr()
                succ = 0
                for p in g:
                    if p[1] > self.di:
                        succ += 1
                sucr = (succ/len(g)) * 100
                if sucr >= self.succ_rate:
                    a[itr_count] = i, j, k, sucr, A_rho
                itr_count += 1
            # prog_count += 1
            # prog_perc = (prog_count/len(S1rng))*100
            # print(prog_perc)   
            if not itr_count%100:
                print(f'itr {itr_count} in {(time.time() - start_time):.3f} s')
        b = np.array(list(a.values())) 
        c = b[np.argmin(b[:][:,-1])]

        print(f"Bumper spacing = {c[0]:.2f}cm\nOuter bumper thickness ="
              f" {c[1]:.2f}cm\nInner bumper thickness = {c[2]:.2f}cm\n"
              f"These parameters give a success rate of {c[3]}% "
              f"and the area density is {c[4]}g.cm^-2")

        return c
                    

class SRL(BLE):
    """
    The SRL BLE is based on a 2-plate shield with a honeycombe between each
    plate.
    This BLE has a 80% success rate in prediction based on 55 experiments
    """

    def crit_d_bal(self, vel, impact_angle=0):
        self.vel = vel
        self.impact_angle = impact_angle
        # define SRL BLE shield constant
        K3S = 1.1
        # determine normal velocity
        theta = self.impact_angle * (np.pi/180)

        dc = (((((self.t_wall ** 0.5 + self.t_b) / K3S) *
                 ((self.Sigma/40) ** 0.5) + self.t_ob) /
            (0.6 * (np.cos(theta) ** (4/3)) * (self.rho_p ** 0.5)
            * self.vel ** (2/3))) ** (18/19))
        return dc


    def crit_d_vap(self, vel, impact_angle=0):
        self.vel = vel
        self.impact_angle = impact_angle
        # define SRL BLE shield constant       
        K3D = 0.4
        # determine normal velocity
        theta = self.impact_angle * (np.pi/180)

        dc = ((1.155 * ((self.S1 ** (1/3)) * ((self.t_b + self.t_wall)
            ** (2/3)) + (self.S2 ** (1/3))
            * (self.t_wall ** (2/3))) * (self.Sigma/70) ** (1/3)) /
            ((K3D ** (2/3)) * (self.rho_p ** (1/3)) * (self.rho_b ** (1/9))
            * (self.vel ** (2/3)) * (np.cos(theta) ** (4/3))))
        return dc


class Christiansen(BLE):
    lbl = 'Christiansen'

    def crit_d_bal(self, vel, impact_angle=0):
        """
        Inputs: proj_vel is the particle velocity in the vapourised region
                the impact angle is a positive angle between 0 and 90 degrees
        Outputs:critcal projectile diameter for given velocity in the
                ballistic region
        """
        self.vel = vel
        self.impact_angle = impact_angle
        # determine angle in radians
        theta = self.impact_angle * (np.pi/180)
        # determine crit diameter
        dc = (((self.t_wall*((self.Sigma/40) ** 0.5) + self.t_ob) /
           (0.6*(np.cos(theta) ** (5/3)) * (self.rho_p ** 0.5)
            * ((self.vel) ** (2/3)))) ** (18/19))
        return dc


    def crit_d_vap(self, vel, impact_angle=0):
        """
        Inputs: proj_vel is the particle velocity in the vapourised region
                the impact angle is a positive angle between 0 and 90 degrees
        Outputs:critcal projectile diameter for given velocity in the
                ballistic region
        """
        self.vel = vel
        self.impact_angle = impact_angle
        # determine angle in radians
        theta = self.impact_angle * (np.pi/180)
        # determine crit diameter
        dc = (3.918 * self.t_wall ** (2/3) * self.S1 ** (1/3)
        * (self.Sigma/70) ** (1/3) * self.rho_p ** (-1/3)
        * self.rho_b ** (-1/9) * (self.vel * np.cos(theta)) ** (-2/3))
        return dc



