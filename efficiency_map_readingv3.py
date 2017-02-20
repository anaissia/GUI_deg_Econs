# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:50:03 2016

@author: anaissia
"""

import csv                         
import pandas as pd
import numpy as np
import math
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d
def read_efficiency(RPM1, Torque):
    eff_map=pd.read_csv('efficiencymap_NF_V3.csv', header=0) 
    T=[2500, 1875, 1250, 625]#5]#5];
    RPM1=[9.55*x for x in RPM1] #from rad/s to RPM
    names = eff_map.columns.values 
    eff=[]
    data={}
    for i in range(0,len(T)):       
        for x in range(0,len(names)-1):
            if str(T[i])==names[x]:
                RPM=[y for y in eff_map[names[x]]]
                eff=[y for y in eff_map[names[x+1]]]
                data[T[i]]=RPM,eff
    x=[]
    y=[]
    
    for k in data.keys():
        actual=[]
        actual2=[]
        x= data[k]
        list1=x[0]
        list2=x[1]
        for i in range(0,len(list1)):
            if ~np.isnan(list1[i]) and ~np.isnan(list2[i]):
                actual.append(list1[i])
                actual2.append(list2[i])
                
        data[k]=actual,actual2
        
    
        

    
  
    # function x=efficiency, y=RPM
    f1=interp1d(data[T[0]][0],data[T[0]][1])  
    f2=interp1d(data[T[1]][0],data[T[1]][1]) 
    f3=interp1d(data[T[2]][0],data[T[2]][1]) 
    f4=interp1d(data[T[3]][0],data[T[3]][1]) 
        
    eff1=f1(RPM1)
    eff2=f2(RPM1)
    eff3=f3(RPM1)
    eff4=f4(RPM1)
    ratio=[]
    upperbound=[]
    tmin=300
    for i in range(0, len(Torque)):
        if Torque[i]<=T[3]:
            ratio.append((Torque[i])/(T[3]))
            upperbound.append(3)
        elif Torque[i]>T[3] and Torque[i]<=T[2]:
            ratio.append((Torque[i]-T[3])/(T[2]-T[3]))
            upperbound.append(2)
        elif Torque[i]>T[2] and Torque[i]<=T[1]:
            ratio.append((Torque[i]-T[2])/(T[1]-T[2]))   
            upperbound.append(1)
        elif Torque[i]>T[1] and Torque[i]<=T[0]:
            ratio.append((Torque[i]-T[1])/(T[0]-T[1]))   
            upperbound.append(0)
        elif Torque[i]>T[0]:
            ratio.append(1)
            upperbound.append(-1)
    
    new_data={}

    list_results1=[]
    list_results2=[]
    list_results3=[]
    list_results4=[]
    list_results5=[]
    for i in range(0,len(Torque)):
        if upperbound[i]==0:
            for x in range(0,len(eff1)):
                list_results1.append((eff1[x]-eff2[x])*ratio[i]+eff2[x])
            
            new_data[Torque[i]]=list_results1
        if upperbound[i]==1:
            for x in range(0,len(eff2)):
                list_results2.append((eff2[x]-eff3[x])*ratio[i]+eff3[x])
            
            new_data[Torque[i]]=list_results2           
        if upperbound[i]==2:
            for x in range(0,len(eff3)):
                list_results3.append((eff3[x]-eff4[x])*ratio[i]+eff4[x])
            
            new_data[Torque[i]]=list_results3             
        if upperbound[i]==3:
            for x in range(0,len(eff4)):
#                list_results4.append((eff4[x])*ratio[i]+0)
                list_results4.append((eff4[x]-80)*(ratio[i])+80)
#                list_results4.append(eff4[x])
            new_data[Torque[i]]=list_results4          
        if upperbound[i]==-1:
            for x in range(0,len(eff4)):
                list_results5.append(eff4[x]-eff4[x]*ratio[i])
            
            new_data[Torque[i]]=list_results5    
            
    return new_data
    

def read_efficiency_BYD(RPM1, T):
    eff_map=pd.read_csv('try BYD eff map.csv') 
    RPM1=[9.55*x for x in RPM1]
    names = eff_map.columns.values 
    motor_speed=[]
    torque=[]
    X=[]
    Y=[]
    effi=[]
    Z=[]
    for x in range(0,len(names)-1,2):
            motor_speed=[y*10 for y in eff_map[names[x]]]
            torque=[y for y in eff_map[names[x+1]]]
            effi=[float(names[x]) for y in eff_map[names[x]]]
            X.extend(motor_speed)
            Y.extend(torque)
            Z.extend(effi)
            if names[x]=='93':
                plt.scatter(motor_speed,torque)
                plt.show()

    for i in range(0,len(X)):
        if i<len(X):
            if np.isnan(X[i]):
                X[i]=0
                Y[i]=0
                Z[i]=0
            
    plt.scatter(X,Y,c=Z)
    plt.show()
    f=interp2d(X,Y,Z,kind='linear')
    eff={}
    for i in T: 
        eff[i]=f(RPM1,i)
    return eff

    
#T=[300,300]
#RPM=[290 for x in T]
#dic={}
#dic=read_efficiency_BYD(RPM,T) 
#print(dic)  
