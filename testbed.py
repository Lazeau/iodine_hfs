# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 13:17:06 2021

@author: matth
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from constants import *

def import_data(f):
    # f = None
    try:
        # f = input('Enter file number:\n')
        filename = "data/{}.DIODE_LIF".format(f)
        data = pd.read_csv(filename, header=34)
    except:
        raise NameError
    
    freq = np.asarray(data.iloc[:,0] + DAQ_OFFSET)  # GHz; must correct offset from DAQ code
    freq -= CORR_OFFSET                             # GHz
    sig = np.asarray(data.iloc[:,1])                # Counts (arb. units)
    
    return freq, sig, f

def find_fwhm(freq, sig):
    halfway_ind = int((sig.shape[0] / 2) - 1)
    index1 = np.argmin(np.abs(0.5 - sig[:halfway_ind]))
    index2 = np.argmin(np.abs(0.5 - sig[halfway_ind:]))
    
    plt.plot(freq[:halfway_ind], sig[:halfway_ind])
    plt.plot(freq[halfway_ind:], sig[halfway_ind:])
    plt.title('TEST FWHM')
    plt.show()
    
    f1 = freq[index1]
    f2 = freq[index2]
    fwhm = np.abs(f2 - f1)
    
    return fwhm

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
signal = pd.DataFrame(None)
frequency = pd.DataFrame(None)
shots = [112, 113, 124, 119, 120, 122]
laser_powers = np.asarray([2.2, 2.2, 1.7, 1.0, 0.7, 0.5]) # mW
laser_intensities = laser_powers / (np.pi * 0.001**2) # 1 mm radius laser beam

for i in range(len(shots)):
    freq, sig, shot = import_data(shots[i])
    sig = np.asarray(sig)
    signal = pd.concat([signal, pd.DataFrame(sig, columns=[shot])], axis=1)
    frequency = pd.concat([frequency, pd.DataFrame(freq, columns=[shot])], axis=1)
    
    plt.plot(freq, sig, label=shot)

plt.legend()
plt.title('Raw Signals')
plt.xlabel('\u0394f (GHz)')
plt.ylabel('Raw Signal (arb.)')
plt.show()

frequency = np.asarray(frequency)
signal = np.asarray(signal)
print(signal.shape)
print(np.max(signal))

# Normalize signal
sig_max = np.amax(signal)
signal = signal / sig_max
for i in range(signal.shape[1]):
    plt.plot(frequency[:,i], signal[:,i], label=shots[i])
plt.legend()
plt.title("Normalized Signals")
plt.xlabel('\u0394f (GHz)')
plt.ylabel('Normalized Signal (arb.)')
plt.show()

FWHM = np.empty(signal.shape[1])
for i in range(signal.shape[1]):
    FWHM[i] = find_fwhm(frequency[:,i], signal[:,i])
    
    # plt.plot(np.log(laser_intensities[i]), FWHM[i], '.', label=shots[i])


for i in range(signal.shape[1]):
    plt.plot(np.log(laser_intensities[i]), FWHM[i], '.', label=shots[i])

plt.legend()
plt.title('FWHM vs. log(I)')
plt.xlabel('log(I_L) (W/m^2)')
plt.ylabel('FWHM (arb.)')
plt.show()