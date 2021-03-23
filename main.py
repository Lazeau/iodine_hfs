# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 15:08:02 2021

@author: mlazo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

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
    
    return freq, sig

def i_ii_line(x, *pars):
    a = pars[0]     # Frequency offset
    Au = pars[1]    # P state magnetic dipole coupling coefficient
    Bu = pars[2]    # P state electric quadrupole coupling coefficient
    Al = pars[3]    # D state magnetic dipole coupling coefficient
    Bl = pars[4]    # D state electric quadrupole coupling coefficient
    T = pars[5]     # Temperature
                    # pars[6:20] are amplitudes for the 15 transition Gaussians
    
    x0 = [-10*Au    + 0.6*Bu            + 12.5*Al   - 0.4910714286*Bl,
          -8.5*Au   + 5.28*Bu           + 12.5*Al   - 0.4910714286*Bl,
          -6*Au     - 0.02*Bu           + 12.5*Al   - 0.4910714286*Bl,
          -8.5*Au   + 5.28*Bu           + 10*Al     - 0.1964285714*Bl,
          -6*Au     - 0.02*Bu           + 10*Al     - 0.1964285714*Bl,
          -2.5*Au   - 0.3333333333*Bu   + 10*Al     - 0.1964285714*Bl,
          -6*Au     - 0.02*Bu           + 6.5*Al    + 0.1035714286*Bl,
          -2.5*Au   - 0.3333333333*Bu   + 6.5*Al    + 0.1035714286*Bl,
          2*Au      - 0.3333333333*Bu   + 6.5*Al    + 0.1035714286*Bl,
          -2.5*Au   - 0.3333333333*Bu   + 2*Al      + 0.2964285714*Bl,
          2*Au      - 0.3333333333*Bu   + 2*Al      + 0.2964285714*Bl,
          7.5*Au    + 1*Bu              + 2*Al      + 0.2964285714*Bl,
          2*Au      - 0.3333333333*Bu   - 3.5*Al    + 0.2375*Bl,
          7.5*Au    + 1*Bu              - 3.5*Al    + 0.2375*Bl,
          7.5*Au    + 1*Bu              - 10*Al     - 0.25*Bl]
    
    G = a + pars[6]*np.exp( -((x-x0[0])^2*M*C^2) / 4*K*T ) + pars[7]*np.exp( -((x-x0[1])^2*M*C^2) / 4*K*T ) + pars[8]*np.exp( -((x-x0[2])^2*M*C^2) / 4*K*T ) + pars[9]*np.exp( -((x-x0[3])^2*M*C^2) / 4*K*T ) + pars[10]*np.exp( -((x-x0[4])^2*M*C^2) / 4*K*T ) + pars[11]*np.exp( -((x-x0[5])^2*M*C^2) / 4*K*T ) + pars[12]*np.exp( -((x-x0[6])^2*M*C^2) / 4*K*T ) + pars[13]*np.exp( -((x-x0[7])^2*M*C^2) / 4*K*T ) + pars[14]*np.exp( -((x-x0[8])^2*M*C^2) / 4*K*T ) + pars[15]*np.exp( -((x-x0[9])^2*M*C^2) / 4*K*T ) + pars[16]*np.exp( -((x-x0[10])^2*M*C^2) / 4*K*T ) + pars[17]*np.exp( -((x-x0[11])^2*M*C^2) / 4*K*T ) + pars[18]*np.exp( -((x-x0[12])^2*M*C^2) / 4*K*T ) + pars[19]*np.exp( -((x-x0[13])^2*M*C^2) / 4*K*T ) + pars[20]*np.exp( -((x-x0[14])^2*M*C^2) / 4*K*T )
    
    return G

freq, sig = read_data()
sig_max = np.amax(sig)
freq_min = np.amin(freq)
freq_max = np.amax(freq)

fits = 100
fit_pars = np.zeros((21,fits))

for i in range(fits):
    fit_pars[0,i] = np.random.random() * freq_max               # Freq offset guess
    fit_pars[1:5,i] = np.random.uniform(freq_min, freq_max, 4)  # Coeff guesses
    fit_pars[5,i] = np.random.random()                          # Temp guess
    fit_pars[6:21,i] = np.random.uniform(0, sig_max, 15)        # Amp guesses
    
    print(fit_pars)
    
    popt, pcov = so.curve_fit(i_ii_line, freq, sig, fit_pars[i])
    plt.plot(freq, i_ii_line(freq, sig, popt))
    plt.plot(freq,sig)