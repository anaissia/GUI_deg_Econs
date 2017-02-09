# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 04:46:53 2017

@author: anaissia
"""

import csv                         
import pandas as pd
import numpy as np
import math
from scipy.interpolate import interp1d
def read_V(SOC):
    data=pd.read_csv('V_vs_SOC_NF.csv') 
    SOC_d = [100-x for x in data['SOC'] ]
    V = [x for x in data['V'] ]
    func=interp1d(SOC_d,V)
    Vout=func(100-SOC)
    return Vout

