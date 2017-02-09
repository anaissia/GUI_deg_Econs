# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 01:03:10 2017

@author: anaissia
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#def write_crate_deg(path):

path='Brampton_23_East_Heavy Duty'
data=pd.read_csv(path+'/crate.csv',header=None) 
time_diff=data[0]
c_rate=data[1]
# Ning's model : 1 C is for 19 A/m2
cell_curr=[x*19 for x in c_rate]
time_add=[0]
for i in range(1,len(time_diff)):
    time_add.append(time_diff[i]+time_add[i-1])
    
f = open(path+'/input_curr.txt', 'w')
dimension=len(time_add)
#for i in range(0,3):
#    f.write(  	0.10451392 0.0	1 3.0 4.2 + Rtotal+  ) 
f.write(str(len(cell_curr)) + "\n")
for i,j in zip(time_add[0:dimension],cell_curr[0:dimension]):
        f.write(str(j)+" "+str(j)+" 	1 3.0 4.2 + Rtotal  \n")
f.close()
    


#plt.plot(c_rate)
#plt.show()
#av=[]
#for i in range(0,len(c_rate)-60):
#    average=np.average(c_rate[i:i+60])
#    for j in range(1,60):
#        av.append(average)
#    
#plt.plot(av)
#plt.show()

