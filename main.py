# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 15:08:02 2021

@author: mlazo
"""
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from scipy.signal import find_peaks

from constants import *

def import_data():
    try:
        f = input('Enter file number:\n')
        filename = "data/{}.DIODE_LIF".format(f)
        data = pd.read_csv(filename, header=34)
    except:
        raise NameError
    
    freq = np.asarray(data.iloc[:,0] + DAQ_OFFSET)  # GHz; must correct offset from DAQ code
    freq -= CORR_OFFSET                             # GHz
    sig = np.asarray(data.iloc[:,1])                # Counts (arb. units)
    
    return freq, sig

def export_data(data, cov):
    try:
        # For pseudorandom guesses, 21 parameters + RMSE
        f = input('Enter export file name for least-squares GUESSES:\n')
        filename = "data/{}.csv".format(f)
        
        df = pd.DataFrame(data)
        df.to_csv(filename, header=['sig_offset', 'dNu1 (GHz)', 'dNu2 (GHz)', 'B_U (GHz)', 'B_L (GHz)', 'T (eV)', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15', 'RMSE'])
        
        # For diagonalized covariance matrix, 21 parameters; no header
        f = input('Enter export file name for least-squares CONVARIANCE MATRIX:\n')
        filename = "data/{}.csv".format(f)
        
        df = pd.DataFrame(cov)
        df.to_csv(filename, header=['sig_offset', 'dNu1 (GHz)', 'dNu2 (GHz)', 'B_U (GHz)', 'B_L (GHz)', 'T (eV)', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15'])
    except:
        raise NameError
    
    return None

def import_excellent():
    try:
        f = input('Enter PARAMETER file name:\n')
        filename = "data/{}.csv".format(f)
        data = pd.read_csv(filename)
        
        f = input('Enter COVARIANCE file name:\n')
        filename = "data/{}.csv".format(f)
        cov = pd.read_csv(filename)
    except:
        raise NameError
    
    RMSE_CUTOFF = 0.007
    excellent_data = data[data['RMSE']<=RMSE_CUTOFF] # 0.007, 0.01
    confidence = cov[data['RMSE']<=RMSE_CUTOFF] # Extract same indices based on data RMSE
    # excellent = data[(data['RMSE']<0.0078) & (data['RMSE']>0.0072)] # investigation of RMSE distribution
    
    return excellent_data, confidence

def export_converged_values(vals, AU, AL, cov):
    try:
        f = input('Enter export file name for CONVERGED VALUES:\n')
        filename = 'data/{}.csv'.format(f)
        
        df = pd.DataFrame(vals)
        df.insert(4, 'A_U (GHz)', AU)
        df.insert(6, 'A_L (GHz)', AL)
        df.to_csv(filename, header=['fit_index', 'sig_offset', 'dNu1 (GHz)', 'dNu2 (GHz)', 'A_U (GHz)', 'B_U (GHz)', 'A_L (GHz)', 'B_L (GHz)', 'T (eV)', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15', 'RMSE'])
        
        
        f = input('Enter export file name for COVARIANCE ARRAY:\n')
        filename = 'data/{}.csv'.format(f)
        
        df = pd.DataFrame(cov)
        df.to_csv(filename)
    except:
        raise NameError
    
    return None

def maxwellian(b, x, x0, T):
    '''
    Convenience function for plotting a Maxwellian curve for our iodine plasma
    '''
    return b * np.exp( -(0.3186334108*((x-x0)**2)/T) )

def i_ii_line(x, *pars):
    '''
    Specify the full lineshape of the 15 hyperfine transitions of I II.
    
    Parameters
    ----------
    x : Numpy array, GHz
        Scanning frequency values
    *pars
        Array of 21 initial parameters to fit the curve
    
    Returns
    -------
    G
        Function containing 15 overlapping Gaussian curves to characterize the
        I II lineshape, in terms of the provided fit parameters.
    '''
    a = pars[0]     # DC offset
    dNu1 = pars[1]  # Highest-intensity transition
    dNu2 = pars[2]  # Second-highest-intensity transition
    Bu = pars[3]    # P state electric quadrupole coupling coefficient
    Bl = pars[4]    # D state electric quadrupole coupling coefficient
    T = pars[5]     # Temperature
                    # pars[6:20] are amplitudes for the 15 transition Gaussians
    
    # Transitions listed from strongest to weakest            Relative Intensity
    x0 = [dNu1,                                               # 100.000
          dNu2,                                               # 75.974
          -0.76*dNu1 + 1.6*dNu2  + 0.37*Bu - 0.2735714286*Bl, # 56.122
          -1.28*dNu1 + 1.8*dNu2  + 0.84*Bu - 0.6439285714*Bl, # 40.087
          -1.56*dNu1 + 1.6*dNu2  + 1.2*Bu  - 0.9664285714*Bl, # 27.551
          -1.6*dNu1  + 1*dNu2    + 1.3*Bu  - 1.128571429*Bl,  # 18.367
          0.68*dNu1  - 3.8*dNu2  - 1.61*Bu + 1.176071429*Bl,  # 16.035
          1.76*dNu1  - 5.6*dNu2  - 2.42*Bu + 2.066428571*Bl,  # 14.813
          -0.16*dNu1 - 2.4*dNu2  - 0.7*Bu  + 0.3335714286*Bl, # 13.994
          3.08*dNu1  - 7.8*dNu2  - 2.86*Bu + 2.86*Bl,         # 9.711
          -0.76*dNu1 - 1.4*dNu2  + 0.1*Bu  - 0.3485714286*Bl, # 9.184
          1.8*dNu1   - 8*dNu2    - 3.15*Bu + 2.153571429*Bl,  # 0.820
          0.64*dNu1  - 5.4*dNu2  - 1.8*Bu  + 0.9514285714*Bl, # 0.714
          3.2*dNu1   - 11*dNu2   - 4.4*Bu  + 3.516071429*Bl,  # 0.510
          4.84*dNu1  - 14.4*dNu2 - 5.28*Bu + 4.926428571*Bl]  # 0.012
    
    G = a + maxwellian(pars[6], x, x0[0], T) + maxwellian(pars[7], x, x0[1], T) + maxwellian(pars[8], x, x0[2], T) + maxwellian(pars[9], x, x0[3], T) + maxwellian(pars[10], x, x0[4], T) + maxwellian(pars[11], x, x0[5], T) + maxwellian(pars[12], x, x0[6], T) + maxwellian(pars[13], x, x0[7], T) + maxwellian(pars[14], x, x0[8], T) + maxwellian(pars[15], x, x0[9], T) + maxwellian(pars[16], x, x0[10], T) + maxwellian(pars[17], x, x0[11], T) + maxwellian(pars[18], x, x0[12], T) + maxwellian(pars[19], x, x0[13], T) + maxwellian(pars[20], x, x0[14], T)
    
    return G

def saturated_lineshape(x, *pars):
    THEOR_REL_INTENSITY = np.asarray([1, 0.75974, 0.56122, 0.40087, 0.27551, 0.18367, 0.16055, 0.14813, 0.13994, 0.09711, 0.09184, 0.00820, 0.00714, 0.00510, 0.00012])
    SAT = 0.5
    SAT_REL_INTENSITY = THEOR_REL_INTENSITY / (1 + (THEOR_REL_INTENSITY / SAT) )
    norm_sat_int = SAT_REL_INTENSITY / np.max(SAT_REL_INTENSITY)
    
    a = pars[0]     # DC offset
    Au = pars[1]    # P state magnetic dipole coupling coefficient
    Al = pars[2]    # D state magnetic dipole coupling coefficient
    Bu = pars[3]    # P state electric quadrupole coupling coefficient
    Bl = pars[4]    # D state electric quadrupole coupling coefficient
    T = pars[5]     # Temperature
    b = pars[6]     # Overall amplitude coefficient
                    # norm_sat_int contains the 15 theoretically predicted relative hf amplitudes
    
    # Transitions listed from strongest to weakest
    x0 = [7.5*Au  + 0.25*Bu  - 10*Al   - 0.25*Bl,
          2*Au    - 0.3*Bu   - 3.5*Al  + 0.2375*Bl,
          -2.5*Au - 0.3*Bu   + 2*Al    + 0.2964285714*Bl,
          -6*Au   - 0.02*Bu  + 6.5*Al  + 0.1035714286*Bl,
          -8.5*Au + 0.330*Bu + 10*Al   - 0.1964285714*Bl,
          -10*Au  + 0.6*Bu   + 12.5*Al - 0.4910714286*Bl,
          -2.5*Au - 0.3*Bu   + 6.5*Al  + 0.1035714286*Bl,
          2*Au    - 0.3*Bu   + 2*Al    + 0.2964285714*Bl,
          -6*Au   - 0.02*Bu  + 10*Al   - 0.1964285714*Bl,
          7.5*Au  + 0.25*Bu  - 3.5*Al  + 0.2375*Bl,
          -8.5*Au + 0.330*Bu + 12.5*Al - 0.4910714286*Bl,
          -2.5*Au - 0.3*Bu   + 10*Al   - 0.1964285714*Bl,
          -6*Au   - 0.02*Bu  + 12.5*Al - 0.4910714286*Bl,
          2*Au    - 0.3*Bu   + 6.5*Al  + 0.1035714286*Bl,
          7.5*Au  + 0.25*Bu  + 2*Al    + 0.2964285714*Bl]
    
    sat_line = a + b * ( maxwellian(norm_sat_int[0], x, x0[0], T) + maxwellian(norm_sat_int[1], x, x0[1], T) + maxwellian(norm_sat_int[2], x, x0[2], T) + maxwellian(norm_sat_int[3], x, x0[3], T) + maxwellian(norm_sat_int[4], x, x0[4], T) + maxwellian(norm_sat_int[5], x, x0[5], T) + maxwellian(norm_sat_int[6], x, x0[6], T) + maxwellian(norm_sat_int[7], x, x0[7], T) + maxwellian(norm_sat_int[8], x, x0[8], T) + maxwellian(norm_sat_int[9], x, x0[9], T) + maxwellian(norm_sat_int[10], x, x0[10], T) + maxwellian(norm_sat_int[11], x, x0[11], T) + maxwellian(norm_sat_int[12], x, x0[12], T) + maxwellian(norm_sat_int[13], x, x0[13], T) + maxwellian(norm_sat_int[14], x, x0[14], T) )
    
    return sat_line

def get_peaks(freq, sig):
    '''
    Determine peaks in the I II lineshape, display each peak with an 'X', and
    print the peak values to the console.
    
    Parameters
    ----------
    freq : array-like, GHz
        Array of scanning frequencies.
    sig : array-like, arb.
        Array of signal values.
    
    Returns
    -------
    None.
    '''
    x = np.linspace(0, sig.size, sig.size)
    peaks, _ = find_peaks(sig, height=1)
    
    plt.plot(x, sig, ".", markersize=1.5)
    plt.plot(peaks, sig[peaks], "x")
    
    print("Peak frequency values (GHz): ", freq[peaks])
    print("Peak signal values (arb.): ", sig[peaks])
    
    return None



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
freq, sig = import_data()
sig_max = np.amax(sig)
freq_min = np.amin(freq)
freq_max = np.amax(freq)

sig = sig / sig_max # Normalize signal

# get_peaks(freq, sig) # Used to determine dNu1 and dNu2 below
dNu1 = 1.62649497      # Highest-intensity transition frequency, in GHz
dNu2 = 0.72610497      # Second-highest-intensity transition frequency, in GHz

# Fitting routine
# fits = 10_000
# results = np.zeros((fits, 22))
# confidence = np.zeros((fits, 21))
# for i in range(fits):
#     fit_pars = np.zeros(21)
    
#     fit_pars[0] = np.random.random() * 0.1                      # Signal offset guess
#     fit_pars[1] = np.random.uniform(dNu1-0.5, dNu1+0.5)             # Highest-intensity transition guess in +-1 GHz window around predicted value, dNu1
#     fit_pars[2] = np.random.uniform(dNu2-0.5, dNu2+0.5)             # Second-highest-intensity transition guess in +-1 GHz window around predicted value, dNu2
#     fit_pars[3:5] = np.random.uniform(-2.5, 2.5, 2)               # Coeff guesses
#     fit_pars[5] = np.random.random() * 0.05                     # Temp guess
#     fit_pars[6:22] = np.random.uniform(0, 1, 15)                # Amplitude guesses
#     print("Fit ", i)
    
#     # Fit
#     LIMITS = ((0, dNu1-0.5, dNu2-0.5, -100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (1, dNu1+0.5, dNu2+0.5, 100, 100, 0.09, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
#     popt, pcov = so.curve_fit(i_ii_line, freq, sig, fit_pars, maxfev=10_000_000, bounds=LIMITS)
#     confidence[i,:] = np.sqrt(np.diagonal(pcov)).reshape((1,21))
#     print(confidence.shape)
    
#     # RMSE
#     RESIDUAL = sig - i_ii_line(freq, *popt)
#     rmse = (np.sum(RESIDUAL**2)/(RESIDUAL.size - 2))**(0.5)
#     print('RMSE:', rmse)
#     print('\npopt:', popt)
#     print('confidence:', confidence[i], '\n')
#     results[i,:] = np.asarray([*popt, rmse])
    
    # Plot the current fit (for testing)
    # fig = plt.figure(figsize=(8,6), dpi=128)
    # plt.plot(freq, sig, ".", markersize=1.5)
    # plt.plot(freq, i_ii_line(freq, *popt))
    # plt.ylabel('Signal (arb.)')
    # plt.xlabel('\u0394f (GHz)')
    # plt.show()


# Export to CSV
# export_data(results, confidence)

# # # # Select excellent fits and calculate coupling coefficients
ex_fits, ex_cov = import_excellent()
print("\nExcellent fits:\n", ex_fits)
num_ex = ex_fits.iloc[:,0].size

# Calculate the values from i_ii_line() for each of the excellent fits
ex_fit_vals = np.zeros((freq.shape[0], num_ex))
for i in range(num_ex):
    ex_fit_vals[:,i] = i_ii_line(freq, *ex_fits.iloc[i,1:22])

# Determine the "best" fit by averaging parameters and values from excellent fits
best_pars = np.asarray(np.mean(ex_fits.iloc[:,1:22]))
best_vals = np.asarray(np.mean(ex_fit_vals, axis=1))

dNu1 = best_pars[1]
dNu2 = best_pars[2]
B_U = best_pars[3]
B_L = best_pars[4]
T = best_pars[5]

# NOTE: below two techniques for finding A_U and A_L produce the same results if one manually averages over the ex_A_U and ex_A_L output to .csv
# A coeffs for each set of excellent params
ex_A_U = (14/25)*np.asarray(ex_fits.iloc[:,2]) - (8/5)*np.asarray(ex_fits.iloc[:,3]) - (31/50)*np.asarray(ex_fits.iloc[:,4]) + (13/25)*np.asarray(ex_fits.iloc[:,5])
ex_A_L = (8/25)*np.asarray(ex_fits.iloc[:,2]) - (6/5)*np.asarray(ex_fits.iloc[:,3]) - (11/25)*np.asarray(ex_fits.iloc[:,4]) + (73/200)*np.asarray(ex_fits.iloc[:,5])
# A coeffs from best_pars
A_U = (14/25)*dNu1 - (8/5)*dNu2 - (31/50)*B_U + (13/25)*B_L
A_L = (8/25)*dNu1 - (6/5)*dNu2 - (11/25)*B_U + (73/200)*B_L

# Standard deviations of each of the 21 parameters and RMSE for all excellent fits
devs = np.std(ex_fits, axis=0)

A_U_err = np.sqrt(( (14/25)*devs[2] )**2 + ( (8/5)*devs[3] )**2 + ( (31/50)*devs[4] )**2 + ( (13/25)*devs[5] )**2)
A_L_err = np.sqrt(( (8/25)*devs[2] )**2 + ( (6/5)*devs[3] )**2 + ( (11/25)*devs[4] )**2 + ( (73/200)*devs[5] )**2)

# Calculate all excellent fit transition predictions from sys. eqs. in terms of dNu1, dNu2, B_U, B_L
ex_dNu = [ex_fits.iloc[:,2],                                                                                       # 100.000 1
       ex_fits.iloc[:,3],                                                                                          # 75.974  2
       -0.76*ex_fits.iloc[:,2] + 1.6*ex_fits.iloc[:,3]  + 0.37*ex_fits.iloc[:,4] - 0.2735714286*ex_fits.iloc[:,5], # 56.122  3
       -1.28*ex_fits.iloc[:,2] + 1.8*ex_fits.iloc[:,3]  + 0.84*ex_fits.iloc[:,4] - 0.6439285714*ex_fits.iloc[:,5], # 40.087  4
       -1.56*ex_fits.iloc[:,2] + 1.6*ex_fits.iloc[:,3]  + 1.2*ex_fits.iloc[:,4]  - 0.9664285714*ex_fits.iloc[:,5], # 27.551  5
       -1.6*ex_fits.iloc[:,2]  + 1*ex_fits.iloc[:,3]    + 1.3*ex_fits.iloc[:,4]  - 1.128571429*ex_fits.iloc[:,5],  # 18.367  6
       0.68*ex_fits.iloc[:,2]  - 3.8*ex_fits.iloc[:,3]  - 1.61*ex_fits.iloc[:,4] + 1.176071429*ex_fits.iloc[:,5],  # 16.035  7
       1.76*ex_fits.iloc[:,2]  - 5.6*ex_fits.iloc[:,3]  - 2.42*ex_fits.iloc[:,4] + 2.066428571*ex_fits.iloc[:,5],  # 14.813  8
       -0.16*ex_fits.iloc[:,2] - 2.4*ex_fits.iloc[:,3]  - 0.7*ex_fits.iloc[:,4]  + 0.3335714286*ex_fits.iloc[:,5], # 13.994  9
       3.08*ex_fits.iloc[:,2]  - 7.8*ex_fits.iloc[:,3]  - 2.86*ex_fits.iloc[:,4] + 2.86*ex_fits.iloc[:,5],         # 9.711   10
       -0.76*ex_fits.iloc[:,2] - 1.4*ex_fits.iloc[:,3]  + 0.1*ex_fits.iloc[:,4]  - 0.3485714286*ex_fits.iloc[:,5], # 9.184   11
       1.8*ex_fits.iloc[:,2]   - 8*ex_fits.iloc[:,3]    - 3.15*ex_fits.iloc[:,4] + 2.153571429*ex_fits.iloc[:,5],  # 0.820   12
       0.64*ex_fits.iloc[:,2]  - 5.4*ex_fits.iloc[:,3]  - 1.8*ex_fits.iloc[:,4]  + 0.9514285714*ex_fits.iloc[:,5], # 0.714   13
       3.2*ex_fits.iloc[:,2]   - 11*ex_fits.iloc[:,3]   - 4.4*ex_fits.iloc[:,4]  + 3.516071429*ex_fits.iloc[:,5],  # 0.510   14
       4.84*ex_fits.iloc[:,2]  - 14.4*ex_fits.iloc[:,3] - 5.28*ex_fits.iloc[:,4] + 4.926428571*ex_fits.iloc[:,5]]  # 0.012   15
ex_dNu = np.asarray(ex_dNu)

# Get standard deviation of each transition from the N_excellent calculations above
trans_devs = np.zeros(15)
for i in range(15):
    trans_devs[i] = np.std(ex_dNu[i], axis=0)

# Standard deviations of each amplitude, extracted from best fits to lineshape
amp_devs = devs[7:22]

print("\n\nEstimated values:\n Signal Offset: {} +- {}\n A_U: {} +- {}\n B_U: {} +- {}\n A_L: {} +- {}\n B_L: {} +- {}\n T: {} +- {}\n dNu1: {} +- {}\n dNu2: {} +- {}\n b1: {} +- {}\n b2: {} +- {}\n b3: {} +- {}\n b4: {} +- {}\n b5: {} +- {}\n b6: {} +- {}\n b7: {} +- {}\n b8: {} +- {}\n b9: {} +- {}\n b10: {} +- {}\n b11: {} +- {}\n b12: {} +- {}\n b13: {} +- {}\n b14: {} +- {}\n b15: {} +- {}\n".format(best_pars[0], devs[1], A_U, A_U_err, B_U, devs[4], A_L, A_L_err, B_L, devs[5], T, devs[6], dNu1, devs[2], dNu2, devs[3], best_pars[6], devs[7], best_pars[7], devs[8], best_pars[8], devs[9], best_pars[9], devs[10], best_pars[10], devs[11], best_pars[11], devs[12], best_pars[12], devs[13], best_pars[13], devs[14], best_pars[14], devs[15], best_pars[15], devs[16], best_pars[16], devs[17], best_pars[17], devs[18], best_pars[18], devs[19], best_pars[19], devs[20], best_pars[20], devs[21]))

# Histogram of RMSE values
# plt.hist(ex_fits.iloc[:,22])

# Plot of measured lineshape and average of best fit values
# Commented lines were used to produce manuscript Fig. 3; may be needed again later
fig = plt.figure(figsize=(8,6), dpi=128)
# font = {'family' : 'normal',
#         'size'   : 16}
# plt.rc('font', **font)
plt.plot(freq, sig, '.', color='c', label='Measured lineshape')
plt.plot(freq, best_vals,  color='r', label='Fitted lineshape')
plt.legend()
plt.ylabel('Signal (arb.)')
plt.xlabel('\u0394f (GHz)')
# plt.xlim(-4, 4)

# Final equations for each hyperfine transition from each "best" parameter
dNu = [dNu1,                                               # 100.000
       dNu2,                                               # 75.974
       -0.76*dNu1 + 1.6*dNu2  + 0.37*B_U - 0.2735714286*B_L, # 56.122
       -1.28*dNu1 + 1.8*dNu2  + 0.84*B_U - 0.6439285714*B_L, # 40.087
       -1.56*dNu1 + 1.6*dNu2  + 1.2*B_U  - 0.9664285714*B_L, # 27.551
       -1.6*dNu1  + 1*dNu2    + 1.3*B_U  - 1.128571429*B_L,  # 18.367
       0.68*dNu1  - 3.8*dNu2  - 1.61*B_U + 1.176071429*B_L,  # 16.035
       1.76*dNu1  - 5.6*dNu2  - 2.42*B_U + 2.066428571*B_L,  # 14.813
       -0.16*dNu1 - 2.4*dNu2  - 0.7*B_U  + 0.3335714286*B_L, # 13.994
       3.08*dNu1  - 7.8*dNu2  - 2.86*B_U + 2.86*B_L,         # 9.711
       -0.76*dNu1 - 1.4*dNu2  + 0.1*B_U  - 0.3485714286*B_L, # 9.184
       1.8*dNu1   - 8*dNu2    - 3.15*B_U + 2.153571429*B_L,  # 0.820
       0.64*dNu1  - 5.4*dNu2  - 1.8*B_U  + 0.9514285714*B_L, # 0.714
       3.2*dNu1   - 11*dNu2   - 4.4*B_U  + 3.516071429*B_L,  # 0.510
       4.84*dNu1  - 14.4*dNu2 - 5.28*B_U + 4.926428571*B_L]  # 0.012
dNu = np.asarray(dNu)

# Plot transition locations with error bars
amps = best_pars[6:]
for i in range(dNu.shape[0]):
    # # # # For generating one transition location at a time for labeling
    # fig = plt.figure(figsize=(8,6), dpi=128)
    # plt.plot(freq, sig, '.', color='c', label='Measured lineshape')
    # plt.plot(freq, best_vals,  color='r', label='Fitted lineshape')
    # plt.ylabel('Signal (arb.)')
    # plt.xlabel('\u0394f (GHz)')
    # # # #
    
    plt.errorbar([dNu[i],dNu[i]], [amps[i],amps[i]], yerr=amp_devs[i], color='k', capsize=5, marker='.')
    plt.errorbar([dNu[i],dNu[i]], [0,0], xerr=trans_devs[i], color='k', capsize=5)
    plt.plot([dNu[i],dNu[i]], [0,amps[i]], color='k')
    
    # # # #
    # plt.show()
    # input('continue...')
    # # # #
plt.show()

# Export converged values for coupling coefficients, signal offset, temperature, amplitudes, and RMSE
ex_converged_vals = np.asarray(ex_fits.iloc[:,:])
export_converged_values(ex_converged_vals, ex_A_U, ex_A_L, ex_cov)

plt.savefig('data/final_result_errbars.png', format='png')



# Pull back in a set of converged values for testing
# test_pars = ex_converged_vals[0,1:22]
# print(test_pars)
# test_vals = i_ii_line(freq, *test_pars)
# plt.plot(freq, sig, label='signal')
# plt.plot(freq, test_vals, label='test line')
# plt.legend()
# plt.show()


# # # # Fitting with saturation with scipy.curve_fit()
#### ! Something erroneous in this implementation; fit produced not physical !
# for i in range(50):
#     fit_pars = np.zeros(7)
#     confidence = np.zeros(7)
    
#     fit_pars[0] = np.random.random() * 0.1              # Signal offset guess
#     fit_pars[1:5] = np.random.uniform(-2.5, 2.5, 4)     # Coeff guesses
#     fit_pars[5] = np.random.random() * 0.05             # Temp guess
#     fit_pars[6] = np.random.uniform(0, 1)               # Amplitude guesses
    
#     print("Fit ", i)
#     LIMITS = ((0, -100, -100, -100, -100, 0, 0), (1, 100, 100, 100, 100, 0.09, 1))
#     popt, pcov = so.curve_fit(saturated_lineshape, freq, sig, fit_pars, maxfev=10_000_000, bounds=LIMITS)
#     confidence = np.sqrt(np.diagonal(pcov))
#     vals = saturated_lineshape(freq, *fit_pars)
    
#     RESIDUAL = sig - saturated_lineshape(freq, *popt)
#     rmse = (np.sum(RESIDUAL**2)/(RESIDUAL.size - 2))**(0.5)
#     print('RMSE:', rmse, '\n')
    
#     plt.plot(freq, vals)
#     plt.title("Reviewer Technique")
#     plt.show()
# # # #