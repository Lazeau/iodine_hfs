# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 20:08:54 2021

@author: mlazo
"""

C = 299792458       # Speed of light, in m/s
K = 1.38064852e-23  # Boltzmann constant, in J/K
H = 6.62607015e-34  # Planck's constant, in J*s

M = 2.107224e-25    # Mass of singly ionized atomic iodine, in kg
# M = 126.9 amu

# Corrected frequency offsets in GHz
# 696.07 nm is DAQ wavelength, 696.0694 nm is corrected wavelength
DAQ_OFFSET = C / (696.07)
CORR_OFFSET = C / (696.0694)