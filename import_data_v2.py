# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:52:56 2016

@author: anaissia
"""
import csv                         #
from datetime import datetime           #import the librairies   
import pandas as pd
def import_time_data(path):
# this function returns the date values as a time step, and the difference delta t in time (seconds)
    time_data=pd.read_csv(path+'/time.csv') 
    date=[]
    values = [x for x in time_data['Time'] ]
    for row in values:
        date.append(datetime.strptime(row, '%d/%m/%Y %H:%M:%S') )   #datetime value in the right date format 18/01/2016 00:00:03
    diff=[]    
    for i in range(0,len(date)-1):
        elapsed_t=date[i+1]-date[i]
        diff.append(elapsed_t.seconds)
        if diff[i]==0:
            diff[i]=0.1  #thats because the data may have issues. So either I delete it.... or do this. simpler.
            
    return date,diff

def import_speed_data(path,unit):
    #retyrbs soeed in m/s
    time_data=pd.read_csv(path+'/speed.csv') 
    if unit==1:
        speed = [x for x in time_data['Speed'] ] #km/h to m/s*0.2777
    else: 
                speed = [x*0.44704 for x in time_data['Speed'] ] # mph to m/s
    diff=[]            
    for i in range(0,len(speed)-1):
      diff.append(speed[i+1]-speed[i])                
    return speed,diff
def import_elevation_data(path):
    #returns elevation in m and diff
  elevation_data=pd.read_csv(path+'/elevation.csv') 
  elevation = [x for x in elevation_data['Elevation (m)'] ]
  dist_GPS=[x for x in elevation_data['Distance (m)'] ]
  diff=[]
  diff_dist_GPS=[]
  for i in range(0,len(elevation)-1):
      diff.append(elevation[i+1]-elevation[i])
      diff_dist_GPS.append(dist_GPS[i+1]-dist_GPS[i])      
  return elevation,diff,diff_dist_GPS
def import_distance(path):  
    #returns distance and diff in distance in m
  data=pd.read_csv(path+'/distance.csv') 
  distance = [x for x in data['Distance'] ]    #km to m
  diff=[]
  for i in range(0,len(distance)-1):
      diff.append(distance[i+1]-distance[i])
  return distance,diff
  
def import_grade(path):
# this function returns the date values as a time step, and the difference delta t in time (seconds)
    grade_data=pd.read_csv(path+'/grade.csv') 
    grade = [x for x in grade_data['Grade'] ]
            
    return grade
  
