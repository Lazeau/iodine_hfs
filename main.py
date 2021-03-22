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
    
    # for testing
    fig = plt.figure(figsize=(8,6), dpi=128)
    plt.plot(freq, sig, ".", markersize=1.5)
    plt.ylabel('Signal (arb.)')
    plt.xlabel('Frequency (Hz)')
    plt.grid(True)
    #print(freq)
    
    return None

# Override
def exp():
    return None

def i_ii_line(x, *pars):
    a = pars[0]     # Frequency offset
    Au = pars[1]    # P state magnetic dipole coupling coefficient
    Bu = pars[2]    # P state electric quadrupole coupling coefficient
    Al = pars[3]    # D state magnetic dipole coupling coefficient
    Bl = pars[4]    # P state electric quadrupole coupling coefficient
    T = pars[5]     # Temperature
    # pars[6:14] are amplitudes for the 15 transition Gaussians
    
    # TODO: insert coefficients from deltaE_hf,i
    x0 = [1,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9,
          10,
          11,
          12,
          13,
          14,
          15]
    
    G = a + pars[6]*np.exp( -((x-x0[0])^2*M*C^2) / 4*K*T ) + pars[7]*np.exp( -((x-x0[1])^2*M*C^2) / 4*K*T ) + pars[8]*np.exp( -((x-x0[2])^2*M*C^2) / 4*K*T ) + pars[9]*np.exp( -((x-x0[3])^2*M*C^2) / 4*K*T ) + pars[10]*np.exp( -((x-x0[4])^2*M*C^2) / 4*K*T ) + pars[11]*np.exp( -((x-x0[5])^2*M*C^2) / 4*K*T ) + pars[12]*np.exp( -((x-x0[6])^2*M*C^2) / 4*K*T ) + pars[13]*np.exp( -((x-x0[7])^2*M*C^2) / 4*K*T ) + pars[14]*np.exp( -((x-x0[8])^2*M*C^2) / 4*K*T ) + pars[15]*np.exp( -((x-x0[9])^2*M*C^2) / 4*K*T ) + pars[16]*np.exp( -((x-x0[10])^2*M*C^2) / 4*K*T ) + pars[17]*np.exp( -((x-x0[11])^2*M*C^2) / 4*K*T ) + pars[18]*np.exp( -((x-x0[12])^2*M*C^2) / 4*K*T ) + pars[19]*np.exp( -((x-x0[13])^2*M*C^2) / 4*K*T ) + pars[20]*np.exp( -((x-x0[14])^2*M*C^2) / 4*K*T )
    
    return G

read_data()