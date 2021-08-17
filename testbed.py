# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 13:17:06 2021

@author: matth
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from constants import *

def import_data():
    f = None
    try:
        f = input('Enter file number:\n')
        filename = "data/{}.DIODE_LIF".format(f)
        data = pd.read_csv(filename, header=34)
    except:
        raise NameError
    
    freq = np.asarray(data.iloc[:,0] + DAQ_OFFSET)  # GHz; must correct offset from DAQ code
    freq -= CORR_OFFSET                             # GHz
    sig = np.asarray(data.iloc[:,1])                # Counts (arb. units)
    
    return freq, sig, f

def fwhm(freq, sig):
    halfway_ind = (sig.shape[0] / 2) - 1
    index1 = np.argnub(np.abs(0.5 - sig[:halfway_ind]))
    index2 = np.argmin(np.abs(0.5 - sig[halfway_indx:]))
    
    return None

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
signal = pd.DataFrame(None)
frequency = pd.DataFrame(None)
shots = []

importing = True
while(importing):
    freq, sig, shot = import_data()
    sig = np.asarray(sig)
    signal = pd.concat([signal, pd.DataFrame(sig, columns=[shot])], axis=1)
    frequency = pd.concat([frequency, pd.DataFrame(freq, columns=[shot])], axis=1)
    shots.append(shot)
    
    plt.plot(freq, sig, label=shot)
    
    option = input("import again?\n")
    if option == "t":
        continue
    else:
        importing = False

plt.legend()
plt.title("Raw Signals")
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
plt.show()

