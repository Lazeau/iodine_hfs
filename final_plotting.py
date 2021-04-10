# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 15:21:19 2021

@author: mlazo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from constants import *

def import_excellent():
    try:
        f = input('Enter file name:\n')
        filename = "data/{}.csv".format(f)
        data = pd.read_csv(filename)
    except:
        raise NameError
    
    print(data)
    excellent = data[data['RMSE']<0.007].iloc[:,2:6]
    
    return excellent

# # # #
fits = import_excellent()
print(fits)