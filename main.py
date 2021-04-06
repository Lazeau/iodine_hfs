# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 15:08:02 2021

@author: mlazo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from scipy.signal import find_peaks
from scipy.signal import argrelextrema

from constants import *

def read_data():
    try:
        f = input('Enter file number:\n')
        filename = "data/{}.DIODE_LIF".format(f)
        data = pd.read_csv(filename, header=34)
    except:
        raise NameError
    
    freq = np.asarray(data.iloc[:,0] + DAQ_OFFSET)  # GHz; must correct offset from DAQ code
    freq -= CORR_OFFSET                 # GHz
    sig = np.asarray(data.iloc[:,1])                # Counts (arb. units)
    
    # for testing
    # fig = plt.figure(figsize=(8,6), dpi=128)
    # plt.plot(freq, sig, ".", markersize=1.5)
    # plt.ylabel('Signal (arb.)')
    # plt.xlabel('Frequency (Hz)')
    # plt.grid(True)
    # print(freq)
    
    return freq, sig

def maxwellian(b, x, x0):
    return mag * np.exp( -((x-x0)**2 * M*C**2) / (4*K*T*x0) )

def i_ii_line(x, *pars):
    a = pars[0]     # Frequency offset ## This is the DC offset, not freq. should be small T.E.S.
    Au = pars[1]    # P state magnetic dipole coupling coefficient 
    Bu = pars[2]    # P state electric quadrupole coupling coefficient
    Al = pars[3]    # D state magnetic dipole coupling coefficient
    Bl = pars[4]    # D state electric quadrupole coupling coefficient
    dNu1 = pars[5]  # Predicted highest-intensity transition frequency, in GHz
    dNu2 = pars[6]  # Predicted second-highesty-intensity transition frequency, in GHz
    T = pars[7]     # Temperature
                    # pars[8:23] are amplitudes for the 15 transition Gaussians
    
    # Transitions listed from strongest to weakest relative intensity in parentheses
    x0 = [7.5*Au            + 0.25*Bu           - 10*Al             - 0.25*Bl,          # 100
          2*Au              - 0.3333333333*Bu   - 3.5*Al            + 0.2375*Bl,        # 76
          -1.628571429*dNu1 - 2.4*dNu2          + 0.39*Bu           - 0.2735714286*Bl,  # 56
          -2.742857143*dNu1 + 1.8*dNu2          + 0.9*Bu            - 0.7475*Bl,        # 40
          3.428571429*dNu1  + 1*dNu2            + 1.3333333333*Bu   - 1.28571429*Bl,    # 18
          1.457142857*dNu1  - 3.8*dNu2          - 1.77*Bu           + 1.176071429*Bl,   # 16
          3.771428571*dNu1  - 5.6*dNu2          - 2.64*Bu           + 2.066428571*Bl,   # 15
          -0.342857143*dNu1 - 2.4*dNu2          - 0.78*Bu           + 0.3335714286*Bl,  # 14
          6.6*dNu1          + 7.8*dNu2          - 3.12*Bu           + 2.86*Bl,          # 10
          -1.628571429*dNu1 - 1.4*dNu2          - 11.39*Bu          - 0.3485714286*Bl,  # 9
          0,
          0,
          0,
          0,
          0]
    
    G = a + maxwellian(pars[8], x, x0[0]) + maxwellian(pars[9], x, x0[1]) + maxwellian(pars[10], x, x0[2]) + maxwellian(pars[11], x, x0[3]) + maxwellian(pars[12], x, x0[4]) + maxwellian(pars[13], x, x0[5]) + maxwellian(pars[14], x, x0[6]) + maxwellian(pars[15], x, x0[7]) + maxwellian(pars[16], x, x0[8]) + maxwellian(pars[17], x, x0[9]) + maxwellian(pars[18], x, x0[10]) + maxwellian(pars[19], x, x0[11]) + maxwellian(pars[20], x, x0[12]) + maxwellian(pars[21], x, x0[13]) + maxwellian(pars[22], x, x0[14])
    
    return G

def get_peaks(freq, sig):
    x = np.linspace(0, sig.size, sig.size)
    peaks, _ = find_peaks(sig, height=1)
    
    plt.plot(x, sig, ".", markersize=1.5)
    plt.plot(peaks, sig[peaks], "x")
    
    print("Peak frequency values (GHz): ", freq[peaks])
    print("Peak signal values (arb.): ", sig[peaks])
    
    return None

# # # #
dNu1 = 1.62649497 # Highest-intensity transition frequency, in GHz
dNu2 = 0.72610497 # Second-highest-intensity transition frequency, in GHz

freq, sig = read_data()
sig_max = np.amax(sig)
freq_min = np.amin(freq)
freq_max = np.amax(freq)

fits = 100

# Not part of fitting routine, but present for reproduceability
# get_peaks(freq, sig)

sig = sig / sig_max # Normalize signal

# Fitting routine
for i in range(fits):
    # fit_pars = np.zeros(21)
    fit_pars = np.zeros(23) # added two spots for dNu1 and dNu2
    
    fit_pars[0] = np.random.random() * 0.1                      # Signal offset guess
    fit_pars[1:5] = np.random.uniform(freq_min, freq_max, 4)    # Coeff guesses ## should be 1:4
    fit_pars[5] = np.random.uniform(dNu1-1, dNu1+1)             # Highest-intensity transition  guess in +-1 GHz window around predicted value, dNu1
    fit_pars[6] = np.random.uniform(dNu2-1, dNu2+1)             # Second-highest-intensity transition guess in +-1 GHz window around predicted value, dNu2
    fit_pars[7] = np.random.random()                            # Temp guess
    fit_pars[8:23] = np.random.uniform(0, sig_max, 15)          # Amplitude guesses
    print(fit_pars)

    popt, pcov = so.curve_fit(i_ii_line, freq, sig, fit_pars)
    plt.plot(freq, i_ii_line(freq, *popt))
    plt.plot(freq,sig)
    plt.show()
    input('break')
    

