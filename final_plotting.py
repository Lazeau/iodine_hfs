# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 15:21:19 2021

@author: mlazo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from constants import *
from main import i_ii_line
from main import import_data

def import_excellent():
    try:
        f = input('Enter file name:\n')
        filename = "data/{}.csv".format(f)
        data = pd.read_csv(filename)
    except:
        raise NameError
    
    print(data)
    excellent = data[data['RMSE']<0.007]
    
    return excellent

fits = import_excellent()
print("Excellent fits:\n", fits)

best_pars = np.zeros(21)
for i in range(22):
    best_pars[i] = np.average(fits.iloc[:,i])

freq, sig = import_data()
fig = plt.figure(figsize=(8,6), dpi=128)
plt.plot(freq, sig, '.', markersize=1.5)
plt.plot(freq, i_ii_line(freq, best_pars))

dNu1 = best_pars[1]
dNu2 = best_pars[2]
B_U = best_pars[3]
B_L = best_pars[4]

A_U = (14/25)*dNu1 - (8/5)*dNu2 - (31/50)*B_U + (13/25)*B_L
A_L = (8/25)*dNu1 - (6/5)*dNu2 - (11/25)*B_U + (73/200)*B_L

print("Predicted coupling coefficients:\n A_U:{}, B_U:{}, A_L:{}, B_L:{}".format(A_U, B_U, A_L, B_L))
