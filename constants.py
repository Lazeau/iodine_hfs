# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 20:08:54 2021

@author: mlazo
"""

C = 299792458       # Speed of light, in m/s
K = 1.38064852e-23  # Boltzmann constant, in J/K

M = 2.107224-25     # Mass of singly ionized atomic iodine, in kg
                    # M = 126.9 amu (Steinberger 2018)

DAQ_OFFSET = C / (696.07e-9)
CORR_OFFSET = C / (696.0694e-9)