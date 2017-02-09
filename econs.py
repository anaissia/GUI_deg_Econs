# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 06:52:11 2017

@author: anaissia
"""
## THIS DOCUMENT IS A CLEAN UP OF THE PREVIOUS CODES

## IT CALCULATES THE ENERGY CONSUMPTION OF AN EBUS FOR A GIVEN ROUTE
## IT CALCULATES THE SOC OF THE BUS
## IT CALCULATES THE C-RATE OF THE CURRENT APPLIED TO THE BATTERIES


import numpy as np
from scipy.integrate import simps
from numpy import trapz
from scipy import integrate
import pandas as pd
from pandas import Series, DataFrame
import math
import pickle
from datetime import datetime
import json
import scipy as sp
import matplotlib.pyplot as plt
import random
from import_data_v2 import import_time_data,import_speed_data,import_elevation_data,import_distance, import_grade
from efficiency_map_readingv3 import read_efficiency

def econs(path,capacity):
    """
    First define general variables
    """
    unit=1 #1 for km/hour or m/s defined in import_data or 2 for US (MPH)
    start_passenger=61
    eff_t=0.95 #transmission efficiency
    eff_convert=0.97 #converter efficiency
    reg_brak_eff=0.5 #how much kinetic energy is recovered
    gear_ratio=5.67 #http://www.binghamton.edu/es2/meetings-conferences/smart_energy_cid_meeting_march_april_2014/Electric_Bus_SBU-UTS_Mamalis.pdf
    radius=0.5715 #michelin 22.5 in diam to m, raidus of the wheels
    M=14864.222 #lbs to kg, CW= no passengers mass of the bus
    A=7.44 #101.5*130*0.0254**2.0 #Frontal area: width * height from in to m
    rho_air=1.225 #air density kg/m3
    Cdrag=0.65 #drag coefficient between 0.6-0.8 for a bus
    aux_load=10000 #auxilliary load in W see Reduction of Accessory Overdrive and Parasitic Loading on a Parallel Electric Hybrid City Bus
    T_max=2500 #maximum allowable torque in Nm depends on motor manufacturer
    P_max=153000 #maximum allowable power in W
    capacity=76#battery capacity in kWh
    """
    Then import all data necessary and store them in arrays
    """
    #date,time_diff=import_time_data(path)
    speed_long,speed_diff_long=import_speed_data(path,unit)
    time_data=pd.read_csv(path+'/time.csv') 
    grade_long=import_grade(path)
    distance_long,diff_long=import_distance(path)
    values = [x for x in time_data['Time'] ]
    time_long=[]
    for row in values:
            time_long.append(row)   #datetime value in the right date format 18/01/2016 00:00:03
    speed=[]
    time=[]
    grade=[]
    speed_diff=[]
    for x in range(0,len(speed_long)-1): 
        speed.append(speed_long[x])  
        speed_diff.append(speed_diff_long[x])
        time.append(time_long[x])
        grade.append(grade_long[x])
    time_diff=[]        
    for i in range(0,len(time)-1):
        time_diff.append(time[i+1]-time[i])
    #elevation,diff_elev,diff_dist_GPS=import_elevation_data(path)
    #distance,diff_distance=import_distance(path)
    
#idea of improvements: randomize passegner weight
        
   #mass of the vehicle including passenger mass (75 kg average)     
    mass=[M+start_passenger*75 for x in speed] 
    # calculate the road slope alpha can be zero if info not given
    
    
    alpha=[np.deg2rad(x) for x in grade]  
    #alpha=[0 for x in speed] 
    
    overal_distance= max(distance_long)/1000 #km
    
    
    """ Calculate all the forces acting on the vehicle"""
    
    c_roll=[0.008+0.00012*i*i for i in speed] #this is from NRC, but can use the road surface/rolling coefficient for smooth concrete  T.Gillespie 1992 
    g=9.81 #m/s2
    
    Fr=[] #rolling force
    Fg=[] #gravity force
    Fa=[0.5*rho_air*A*Cdrag*(x**2) for x in speed] #aerodynamic force 
    for i in range(0,len(alpha)):
        Fr.append(c_roll[i]*mass[i]*g*math.cos(alpha[i]))
        Fg.append(mass[i]*g*math.sin(alpha[i]))
        
    """ Total force and acceleration to find motor torque """
       
    acceleration=[]
    Fmotor_wheel=[]
    Ftot=[]
    T_wheel=[]
    T_motor=[]
    RPMw=[]
    RPM=[]
    for i in range(0, len(speed)-1):
        acceleration.append((speed_diff[i])/(time_diff[i]))
           

    # used calculation in "Simulation of an electric transportation system at The Ohio State University"
    for i in range(0, len(speed)-1):
        Ftot.append(Fa[i]+Fr[i]+Fg[i])
        Fmotor_wheel.append((mass[i]+0.1*mass[i])*acceleration[i]+Ftot[i]) #0.1 mass to account for inertia of internal components 
    
        T_wheel.append(Fmotor_wheel[i]*radius)
        T_motor.append(T_wheel[i]/(eff_t*gear_ratio)) 
        RPMw.append(speed[i]/radius)
        RPM.append(RPMw[i]*gear_ratio)
        if speed[i]==0 and speed[i+1]==0: #this makes sure that if the vehicle is idle there is no torque...
            T_motor[i]=0
    """The following for loop ensures that 
        the vehicle completes a feasible path, 
        eg the maximum allowable torque can never be reached
    """
    for i in range(0,len(T_motor)-1): 
        
    
        if T_motor[i]>T_max:
            T_motor[i]=T_max
            T_wheel[i]= (eff_t*gear_ratio)*T_motor[i]
            Fmotor_wheel[i]=T_wheel[i]/radius
            
            coeff1=(0.5*rho_air*A*Cdrag+0.00012*mass[i]*g*math.cos(alpha[i]) )
            coeff2=(mass[i]+0.1*mass[i])
            coeff3=(0.008*mass[i]*g*math.cos(alpha[i])+Fg[i]-(mass[i]+0.1*mass[i])*speed[i])-Fmotor_wheel[i]
            coeff=[coeff1,coeff2,coeff3]
            results=np.roots(coeff) #solves the quadratic equation
            # d=speed diff because time diff is = 1 second
            d=speed_diff[i]

            # recalculate the delta t (longer) and speed to complete drive cycle with a maximum torque
            speed[i+1]=results[1]
            speed_diff[i]=speed[i+1]-speed[i]
            
            delta_t=abs(d/speed_diff[i])

            if delta_t>1:
                time_diff[i]=delta_t
            else:
                time_diff[i]=1
              
            c_roll[i+1]=0.008+0.00012*speed[i+1]**2 
            acceleration[i]=speed_diff[i]/time_diff[i]
            RPM[i+1]=speed[i+1]*gear_ratio/radius
    
        if T_motor[i]<-T_max: #in case of braking, same method than when moving forward
            T_motor[i]=-T_max
            T_wheel[i]= (eff_t*gear_ratio)*T_motor[i]
            Fmotor_wheel[i]=T_wheel[i]/radius        
            coeff1=(0.5*rho_air*A*Cdrag+0.00012*mass[i]*g*math.cos(alpha[i]) )
            coeff2=(mass[i]+0.1*mass[i])
            coeff3=(0.008*mass[i]*g*math.cos(alpha[i])+Fg[i]-(mass[i]+0.1*mass[i])*speed[i])-Fmotor_wheel[i]
            coeff=[coeff1,coeff2,coeff3]
            results=np.roots(coeff)
            d=speed_diff[i]
            speed[i+1]=results[1]
            speed_diff[i]=speed[i+1]-speed[i]
            
            delta_t=abs(d/speed_diff[i])
            if delta_t>1:
                time_diff[i]=delta_t
            else:
                time_diff[i]=1
            c_roll[i+1]=0.008+0.00012*speed[i+1]**2 
    
            acceleration[i]=speed_diff[i]/time_diff[i]
            RPM[i+1]=speed[i+1]*gear_ratio/radius
            
    
    """ Recalculate the variable time (increased) """
    kk=1
    time[0]=0
    for i in range(0,len(time_diff)):
        time[kk]=time[kk-1]+time_diff[i]
        kk=kk+1
        
       
    def efficiency_map(RPM,TM):
        # returns the efficiency in % of the motor for a given torque and RPM. 
        RPM1=[RPM,RPM]
        Torque=[TM,TM]
        RPM={}
        eff=read_efficiency(RPM1,Torque)
        return eff[TM][0]
    def plotting_function(x,y,xaxis,yaxis,name,xlimv,xlimv2,ylimv,ylimv2):
        # creates consistent shapes and layout for graphs, I like this type of layout :-)
        ax = plt.subplot(111)
        ax.plot(x,y)
        plt.xlabel(xaxis,fontsize=18)
        plt.ylabel(yaxis,fontsize=18)
        if xlimv !=0 or xlimv2 !=0:        
            plt.xlim(xlimv,xlimv2)
        if ylimv !=0 or ylimv2 !=0:             
            plt.ylim(ylimv,ylimv2)
    
    
         
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.grid()    
        plt.savefig(path+'/'+name+'.pdf', bbox_inches='tight')    
        plt.show() 
    def plotting_scatter(x,y,xaxis,yaxis,name,xlimv,xlimv2,ylimv,ylimv2):
        #same than before but for scatter plots
        ax = plt.subplot(111)
        ax.scatter(x,y)
        plt.xlabel(xaxis,fontsize=18)
        plt.ylabel(yaxis,fontsize=18)
        if xlimv !=0 or xlimv2 !=0:        
            plt.xlim(xlimv,xlimv2)
        if ylimv !=0 or ylimv2 !=0:             
            plt.ylim(ylimv,ylimv2)
    
    
         
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.grid()    
        plt.savefig(path+'/'+name+'.pdf', bbox_inches='tight')    
        plt.show()         
        
    
    P_instantaneous_w=[] #power consumed instanteously at the wheel
    P_instantaneous=[] #power consumed instanteously at the electric motor
    eff=[]
    
    for x in range(0, len(speed)-1)    :
        
        if RPM[x]<2:
            eff.append(0.2) #the graph doesn't exist for very small RPMs
        else:
            eff.append(efficiency_map(RPM[x],abs(T_motor[x]))/100) #0.104719755 rpm to rad/s
        
        if acceleration[x]>=0 and T_motor[x]>0: # if the vehicle is moving forward
            
         P_instantaneous_w.append(T_motor[x]*RPM[x])
         
         P_instantaneous.append((P_instantaneous_w[x]/(eff[x]*eff_convert))+aux_load)  # get the power consumed by the motor
    
        else: #if the vehicle brakes
            P_instantaneous_w.append(T_motor[x]*RPM[x]*reg_brak_eff)    
            P_instantaneous.append((P_instantaneous_w[x]/(eff[x]*eff_convert))+aux_load)# get the power consumed by the motor with regenerative braking
        if P_instantaneous[x]>P_max: #just to make sure one last time that you never go beyond max power 
            P_instantaneous[x]=P_max
        if P_instantaneous[x]<-P_max:
            P_instantaneous[x]=-P_max
                    
    
    energy=trapz(P_instantaneous[1:len(P_instantaneous)-1],time[1:len(P_instantaneous)-1])/(1000*overal_distance*3600) #integral of power over time = energy  kwh/km
    plotting_function(time[0:len(P_instantaneous)],P_instantaneous,'Time (s)','P_instantaneous_w','P_instantaneous_w',0,max(time),0,0) 
    plotting_function(time[0:len(P_instantaneous)],RPM,'Time (s)','RPM','RPM',0,max(time),0,0) 
    plotting_function(time[0:len(P_instantaneous)],T_motor,'Time (s)','T_motor','T_motor',0,max(time),0,0) 
    plotting_function(time[0:len(P_instantaneous)],acceleration,'Time (s)','acceleration','acceleration',0,max(time),0,0) 
    plotting_function(time,Fa,'Time (s)','F drag (N)','fa',0,max(time),0,0)
    plotting_function(time,Fr,'Time (s)','F roll (N)','fr',0,max(time),0,0)
    plotting_function(time,Fg,'Time (s)','F gravity (N)','fg',0,max(time),0,0)
        
    print('Energy per kwh/km',energy)
    print('Total energy used ',energy*overal_distance)             
    energy_kwh= [time_diff[i]*(P_instantaneous[i])*0.000278*0.001 for i in range(1,len(P_instantaneous)-1)]
    ## ENERGY STORAGE CALCULATIONS
    SOC_i=0.95
    
    #SOC_discharge is the amount of capacity discharge at each delta t
    SOC_discharged=[energy_kwh[i]/capacity for i in range(0,len(energy_kwh))] 
    
    eff_SOC= 1 #need a converter efficiency of some sort
    SOC=[]
    soc_current=SOC_i
    
    
    t=[x/60 for x in time]
    
    """ Calculate the SOC of the battery in a really simple way, may need improvements""" 
    for i in range(0,len(SOC_discharged)):
        soc_current=soc_current-(1)*SOC_discharged[i]
        SOC.append(soc_current)
    SOC=[x*100 for x in SOC]    
    
    """ Estimate battery voltage for each SOC"""
    

econs('Brampton_23_East_Heavy Duty',200)