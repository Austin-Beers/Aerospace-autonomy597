#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 17:42:33 2021

@author: Austin Beers
"""
# Question3

import math
p_pascals = [101000, 94200, 69700, 46600, 30100]
t_k = [288.85, 284.45, 268.38, 248.56, 228.8]
t_0 = 288.85
p_0 = 101325
rho_0 = 1.225
gamma = 1.4
m_array=[]
c_s_array=[]
v_t_array=[]
products=[]
for i in p_pascals:
    m = math.sqrt(5.00 * ((((584.46/i) + 1.00)**(2.00/7.00))- 1.00))
    m_array.append(m)
    
    print(m)
    
for j in t_k:
    c_s = math.sqrt(gamma * (p_0/rho_0) * (j/t_0))
    c_s_array.append(c_s)
    


for num1, num2 in zip(c_s_array, m_array):
	products.append(num1 * num2)



