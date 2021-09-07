#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 23:52:40 2021

@author: Austin Beers
"""
import math
from math import radians, cos, sin, asin, sqrt
kunv_y = 40.8503
kunv_x = 77.8460
dest_coords = {"bellefonte": [40.8896, 77.8134], "KDCA": [38.8534, 77.0452], "KLAX": [33.9439, 118.4085], "EDDF": [50.0375, 8.5635]}
a_array = []
hav_array=[]


for key in dest_coords:
    tempArray = dest_coords[key]
    a = (math.sqrt(((tempArray[1] - kunv_x)**2) + ((tempArray[0] - kunv_y)**2))* 60)
    a_array.append(a)
    
    



def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

for newkey in dest_coords:
    newTempArray = dest_coords[newkey]
    
    hav = haversine(kunv_x, kunv_y, newTempArray[1], newTempArray[0]) * .539957 #convert to nautical miles
    # 0.539957
    hav_array.append(hav)
    
