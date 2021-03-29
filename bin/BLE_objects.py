# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 11:11:15 2020

@author: dalgl
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from itertools import product
import time, random
import directories as dirs
from BLE_func import read_master_data


class BLE:
    def __init__(self, shieldsheet, projsheet):
        """
        This function defines the variables from a given shield sheet and
        projectile sheet
        """
        varis = pd.DataFrame(pd.read_excel(
            dirs.in_dir+'ShieldProperties.xlsx',
            sheet_name=shieldsheet, index_col=0,
            header=0, dtype={'denotion': str,
                             'Value': float}, engine='openpyxl'))
        # import data from excel sheet and convert to numpy array
        varip = pd.DataFrame(pd.read_excel(
            dirs.in_dir+'ProjectileProperties.xlsx',
            sheet_name=projsheet, index_col=0, header=0,
            dtype={'Value': float}, engine='openpyxl'))

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
        # define label and plot colour no.
        self.col = random.randint(0,255)/255
        self.lbl = shieldsheet
        # define minimum and maximum velocities to scan over
        self.maxvel = 16
        self.minvel = 0.001
        # set shatter region limits
        # ideally shock analysis done to work this out
        self.v_shat = 3 # kms^-1
        self.v_vap = 7 # kms^-1


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
        #start_time = time.time()
        vd = np.array(self.arr())
        cmap = plt.cm.rainbow
        plt.plot(vd[:, 0], vd[:,1], color=cmap(self.col), label=self.lbl)
        plt.axis([0, self.maxvel, 0.1, 0.8])
        plt.xlabel('Velocity (km.s-1)')
        plt.ylabel('Crit Diameter (cm)')
        plt.legend()
        #print(f'Time: {(time.time() - start_time):.3f} s')


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
    
    
    def u_succ(self, MASTER_arr):
        """
        Parameters
        ----------
        MASTER_arr : Array of floats
            An array output from MASTER analysis. In first column is impactor 
            velocity,in second column is impactor diameter and in the third
            column is collision probability. Length of this array is
            independant of shield_arr.
        shield_arr : List
            list array of two columns with impactor velocity and crtical
            diameter in the first and second columns respectively.

        Returns
        -------
        None.

        """
        cum_prob = 0
        g = np.array(self.arr())
        for i in MASTER_arr:
            # find the point in shield_arr with the highest velocity below that of the point in MASTER_arr
            n = 0
            v = i[0]
            idx = (np.abs(g[:,0]-v)).argmin()
            if g[idx][0]<= v:
                p1 = g[idx]
                p2 = g[idx+1]
            else:
                p1 = g[idx-1]
                p2 = g[idx]
            m = (p2[1] - p1[1])/(p2[0] - p1[0])
            if i[1] <= (p1[1] + m * (v-p1[0])):
                cum_prob += i[3]
        return cum_prob


    def ideal(self, succ_rate, Arho_max = 3, bumper_dependancy = 0, inc_wall = True, succ_it = 1, file = 'master_t.__1'):
        """
        The definition of 'ideal' used in this function is:
            An ideal shield should have the lowest possible area density while
            preventing perforation caused by a given diameter particle that
            could be travelling over a range of velocities.
        Parameters
        ----------
        succ_rate : int
            The percentage of velocities over the velocity range for which 
            perforation has been successfully prevented. Given in %.
        Arho_max : int, optional
            A filter. The maximum area density that function should scan over.
            The default is 3. Given in g.cm^-2.
        bumper_thickness_dependancy : int, optional
            The ratio of thickness between the outer and inner bumper. If both
            bumper thicknesses should be independant set to 0. 
            If multiple bumpers: enter list of floats where the first number
            corresponds to thickness ratio between fist and second bumper,
            the second number corresponds to second and third e.t.c.
            Independant ratios in this list should be set to 0.
            The default is 0.

        Returns
        -------
        Array of floats64
            DESCRIPTION.

        """
        #start_time = time.time()
        # save starting shield properties
        S1 = self.S1
        t_ob = self.t_ob
        t_b = self.t_b
        # define % of non-perforation over velocity range
        self.succ_rate = succ_rate
        self.Arho_max = Arho_max
        # define output list containing successful variables
        a = {}
        # create probabilities array from MASTER data
        if succ_it == 1:
            MASTER_arr = read_master_data(file)[0]
        # define shield values to vary and the ranges of variation
        S1rng = np.arange(start = 0.1, stop = 12, step = 0.1)
        t_obrng = np.arange(start = 0.01, stop = 5, step = 0.05)
        if self.t_b == 0: # allows analysis of whipple shields
            t_brng = np.zeros(1)
            rho_ratio = 0
        else: # used for shields with internal structure like a stuffed Whipple
            t_brng = np.arange(start = 0.1, stop = 5, step = 0.05)
            # approximate ratio between average inner structure density and 
            # bumper density
            rho_ratio = 0.1
        # is wall considered in area density calculation?
        if inc_wall == False:
            wall = 0
        elif inc_wall == True:
            wall = self.t_wall
        else:
            print('ERROR: inc_wall poorly defined. Define as either True or False')
            return 0
        # iterate for all combinations
        itr_count = 0
        # define error counters
        Arho_ct = 0
        succ_ct = 0
        for i, j, k in product(S1rng, t_obrng, t_brng):
            self.S1 = i
            self.t_ob = j
            if bumper_dependancy == 0:
                self.t_b = k
            else:
                self.t_b = bumper_dependancy * j
            A_rho = self.rho_b * (j + k + rho_ratio * i + wall)
            # filter out area ratios greater than Arho_max
            if A_rho <= self.Arho_max:
                Arho_ct += 1
                # use simple success criteria
                if succ_it == 0:
                    g = self.arr()
                    succ = 0
                    for p in g:
                        if p[1] > self.di:
                            succ += 1
                    sucr = (succ/len(g)) * 100
                # use probability based criteria using MASTER data
                elif succ_it == 1:
                    sucr = self.u_succ(MASTER_arr)*100
                    #print(f'{sucr}')
                #else:
                    #print('ERROR: Invalid succ_it value')
                    #return 0
                # determine whether condition satisfies the minimum success rate
                if sucr >= self.succ_rate:
                    a[itr_count] = self.S1, self.t_ob, self.t_b, sucr, A_rho
                    succ_ct += 1
                itr_count += 1
                #print(f'Performing Iterations: {itr_count}', end="\r")
                
        # reassign original shield properties
        self.S1 = S1
        self.t_ob = t_ob
        self.t_b = t_b
        if Arho_ct == 0:
            print('ERROR: No combination with Arho <= Arho_max\n'
                  'Try increasing Arho_max')
            return 0
        elif succ_ct == 0:
            print('ERROR: No condition with this Arho_max was successful\n'
                  'Try decreasing succ_rate or increasing Arho_max')
            return 0
        else:
            b = np.array(list(a.values()))
            c = b[np.argmin(b[:][:,-1])][4]
            # find all conditions with lowest area density then find the greates success rate among them
            d = []
            for i in b:
                if (i[4]-c)/c <= 0.05: # include all area densities within 5% of min value
                    d.append(i)
            e = np.array(d)
            ans = e[np.argmax(e[:][:,-2])]
#            print(f'itr {itr_count} in {(time.time() - start_time):.3f} s')
#            print(f"Bumper spacing = {ans[0]:.2f}cm\nOuter bumper thickness ="
#                  f" {ans[1]:.2f}cm\nInner bumper thickness = {ans[2]:.2f}cm\n"
#                  f"These parameters give a success rate of {ans[3]:.2f}% "
#                  f"and the area density is {ans[4]:.2f}g.cm^-2")
            return ans


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


class Whipple(BLE):
    
    #def graph_prop(self):
     #   self.lbl = 'Whipple'

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



