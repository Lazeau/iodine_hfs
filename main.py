# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 15:08:02 2021

@author: mlazo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from constants import *

# function to read data, use pd.read_csv(file, header=[rows])

def read_data():
    try:
        f = input()
        filename = "data/{}.DIODE_LIF".format(f)
        data = pd.read_csv(filename, header=34)
    except:
        raise NameError
    
    freq = data.iloc[:,0] * 10e+9 + DAQ_OFFSET
    freq -= CORR_OFFSET
    
    sig = data.iloc[:,1]
    
    
    fig = plt.figure(figsize=(8,6), dpi=128)
    plt.plot(freq, sig, ".", markersize=1.5)
    plt.ylabel('Signal (arb.)')
    plt.xlabel('Frequency (Hz)')
    plt.grid(True)
    
    return None

def fit():
    return None

read_data()