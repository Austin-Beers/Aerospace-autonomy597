#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 21:49:02 2021

@author: TheDongBringer
"""


import numpy as np
import math

master_data_dictionary = {
        "config1" : {
                1: [0,0],
                2: [0, 90],
                3: [45, 185],
                4: [45, 275],
                5: [45, 45],
                6: [85, 0]
            },
        "config2": {
                1: [0, 0],
                2: [0, 90],
                3: [45, 185],
                4: [45, 275]
            },
        "config3": {
                1: [0, 0],
                2: [0, 90],
                3: [15, 185],
                4: [15, 275]
            },
        "config4": {
                1: [45, 180],
                2: [45, 0],
                3: [45, 15],
                4: [85, 150]
            }
    }

master_matrix_dictionary = {}
master_gdop_array = []

# for loop to create h (dcm) matrix for each config key in the dictionary
for key in master_data_dictionary:
    # creates temporary pointer to current config
    temp_dict = master_data_dictionary[key]
    temp_matrix = []
    # for loop to access each item in nested dictionary
    for keyTwo in temp_dict:
        # creates temporary pointer to both elevation and azimuth in current value array
        el = math.radians(temp_dict[keyTwo][0])
        az = math.radians(temp_dict[keyTwo][1])
        # creates a matrix row for each satelite
        matrix_row = [math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el), -1]
        # appends row to total matrix representing whole config
        temp_matrix.append(matrix_row)

    # assignes current total config matrix to current key in a new dictionary
    master_matrix_dictionary[key] = temp_matrix
    
# for loop to go through each configuration matrix
for mtx in master_matrix_dictionary:
    # creates a temporary pointer to the current h matrix
    temp_pointer = master_matrix_dictionary[mtx]
    # transpose of h
    transpose = np.transpose(temp_pointer)
    # multiplies h times its transpose
    product = np.dot((transpose), (temp_pointer))
    # takes the inverse of the resultant product
    inverse = np.linalg.inv(np.matrix(product))
    # print(product)
    #  takes the trace of the inverse and square roots it
    GDOP = math.sqrt(inverse.trace())
    # print(GDOP)
    print(GDOP)
    
    


