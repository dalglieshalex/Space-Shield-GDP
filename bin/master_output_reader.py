# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 12:24:53 2021

@author: Ollie Marston
"""

import numpy as np
import math

"""
INPUTS
These will be what the function takes as input once this script is converted to a function
"""

"""Define the MASTER output file to be read from"""
path = 'C:/Users/dalgl/OneDrive/Documents/default/output/'

"""Define cross-sectional area of satellite"""
A = 5       # m^2

"""Define the time on orbit for the satellite"""
t = 9       # years

"""INPUTS"""
def read_master_data(results_file_path):
    f = open(results_file_path, "r")     # Read output file

    vdfarr = [0, 0, 0]          # Set-up "Velocity, Diameter, Flux Array" (hence 'vdfarr')
    vdparr = [0, 0, 0]          # Set-up "Velocity, Diameter, Probability Array" (hence 'vdparr')

    content = f.readlines()
    f.close()
 
    for lineno in content:        # Cycle through each line of text file in turn ("lineno" is for "Line Number")
        try:
            linearr = np.array(lineno.split())                  # Splits line up into a list with each entry delineated by spaces ("linearr" is for "Line Array")
            linearr = linearr[[0,1,-1]]                         # Takes the first entry (velocity column), second entry (diameter column) and last entry...
                                                                # ...(total flux per m^2 per year column) and keeps these in an array
            linearr = np.asarray(linearr, dtype=np.float32)     # Changes the array from an array of string entries to an array of float entries
            vdfarr = np.vstack([vdfarr, linearr])               # Appends this array as a new row onto the existing VDF array
            
            F = linearr[2]                                      # Takes the flux value 
            p = 1 - math.exp((-1)*F*A*t)                        # Calculates the probability of this impact
            linearr2 = np.array([linearr[0], linearr[1], p])    # Makes an array of velocity, diameter and probability for this impact
            vdparr = np.vstack([vdparr, linearr2])              # Appends this array as a new row onto the existing VDP array
        
        except:             # Catch errors which we expect from the lines of the text file which don't contain impact information. Either these will be...
                            # ... due to the line not having spaces to delineate by, the first column being a '#' or there not being any text on the row
            pass
    
    # vdfarr = np.delete(vdfarr, 0, 0)        # Get ride of the first row of each array which is just zeros and was just used to set-up the array to append onto.
    # vdparr = np.delete(vdparr, 0, 0)        # " "
    
    # As the arrays are so long, this prevents the whole thing from being printed, so that the arrays are truncated to the first and last three entries...
    # ...for printing.
    np.set_printoptions(threshold=20)       
    
    # Print the two arrays
    return vdfarr, vdparr