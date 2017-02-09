#This code intends to run the dualfoil5.2.f Newman's model to include battery degradation
#with cycling
# V1: 08-04 Added change in concentration of lithium at negative electrode due to degradation
#V2: 08-05 Modified the calculation of the SEI growth rate and  volume-averaged capacity lost due to parasitic reaction
# so that its only acting on the discharge
#V3: Ning's input as of 12/08 can't get right current input
#V4: Validation for cycle 1 with Ning's experiment
#V5: trying to also change the thickness of electrode as SEI increases. Did not turn out great. Try to change the capaciy NEGATIVE
#V6: change also the overall resistnace to add the negative resistance.
#V7: now try model ebus cell
#V8:  Add a calculation for j_side myself
import numpy as np
import time, sys, os
import pandas as pd
from numpy import f2py
import math
import shutil
import scipy
import timeit
import matlab.engine

directory='C:/Users/anaissia/Documents/Degradation Model'
def IsNumeric(lst):
	#FOR EACH VALUE FROM 0 TO LENGTH OF LINE
	for x in (0, len(lst)):
		try:
			#TRUE IF NUMERICAL
			float(lst[x])
			return True	
		except ValueError:
			#FALSE IF NOT NUMERICAL
			return False			
   

def write_input_files(R,SOC_n):
    As=
    as1=
    as3=
    C_nom=
    ce_avg=
    Cp=
    cs1_max=
    cs3_max=
    diam=
    Ds1_ref=
    Ds3_ref=
    Ea_Ds1=
    Ea_Ds3=
    Ea_k1=
    Ea_k3=
    eps1s=
    eps3s=
    F=
    h=
    height=
    k1_ref=
    k3_ref=
    L=
    R=
    Rc=
    rho=
    Rs1=
    Rs3=
    SA_V=
    T_amb=
    T_ref=
    thick1=
    thick2=
    thick3=
    Vc=
    x1_soc0=
    x1_soc1=
    y3_soc0=
    y3_soc1=
    data_matlab=[As,\
        as1,\
        as3,\
        C_nom,\
        ce_avg,\
        Cp,\
        cs1_max,\
        cs3_max,\
        diam,\
        Ds1_ref,\
        Ds3_ref,\
        Ea_Ds1,\
        Ea_Ds3,\
        Ea_k1,\
        Ea_k3,\
        eps1s,\
        eps3s,\
        F,\
        h,\
        height,\
        k1_ref,\
        k3_ref,\
        L,\
        R,\
        Rc,\
        rho,\
        Rs1,\
        Rs3,\
        SA_V,\
        T_amb,\
        T_ref,\
        thick1,\
        thick2,\
        thick3,\
        Vc,\
        x1_soc0,\
        x1_soc1,\
        y3_soc0,\
        y3_soc1]
    return data_matlab
def calc_j_side():
    os.getcwd()
    os.chdir('%s' % directory)
    filename='potential.out'
    data=np.loadtxt(filename)
    Upara=0.38 #V
    j0=8*10**-7
    F=96485.3329 #faraday's constant
    R=8.3144598 # J/mol k
    overpotential=[]
    alpha_c=0.5
    jpara=[]
    time=[]
    for i in range(0,len(data)-1):
        overpotential.append(data[i][1]-Upara)
        jpara.append(-j0*math.exp(alpha_c*F*overpotential[i]/(R*data[i][2])))
        time.append(data[i][0])
    
    return jpara,time
def degradation(Rf,Rtotal,j_para,time):

      #get the rate of SEI growth
      # Ning et al. Cyclce life modeling of lithiunm ion batteries
  	rho=2.1*1000*1000 #kg/m3
#  	rho=rho*(10**-2)**3 #g/m3
  	kappa=3.79*(10**-7) #S/m
      
  	M=0.1#kg/mol
  	F=96485 #C/mol (C=A.s)
  	test=[]   
  	thickness_SEI=[]

               
  	thickness_SEI_tot=abs(scipy.integrate.simps(j_para,time)*M*60/(F*rho)) # because integral

  	Rcycle=thickness_SEI_tot/kappa

  	R= Rf+Rcycle #ohm.m2
  	Rtotal = Rtotal+Rcycle 	
      
  	return R,Rtotal
   
def decrease_concentration(Cs,j_para,time):
  	eps_e=0.440 #electrolyte volume frac
  	eps_fl=0.07 #conductive filler volume frac 
  	rs=2*10**-6 #radius particule m
  	eps_s=(1-eps_e-eps_fl) 
  	a_neg=3*eps_s/rs
  	F=96485 #C/mol (C=A.s) Faraday
      
  	# Calculate volume-averaged capacity lost due to parasitic reaction (C/m3)

#  	for i in range(0,len(j_para)-1):
      #     if data[i][1]<0 :
#              Q_1.append(-(data[i][1]*a_neg*(data[i+1][0]-data[i][0])*60)) #*(10**-4)
#                         Q_1.append(abs(scipy.integrate.simps(j_para,time)*60*a_neg))
#           if data[i][1]>0 :
#              Q_1.append((data[i][1]*a_neg*(data[i+1][0]-data[i][0])*60)) #*(10**-4)           
  	Q_tot=abs(scipy.integrate.simps(j_para,time)*60*a_neg)#  because integral
  	Ccycle=Cs-abs(Q_tot)/(eps_s*F)  
  	print('charge',abs(Q_tot)/(eps_s*F))
  	return Ccycle   
def run_cycles(N):
    	owd = os.getcwd()
    	os.chdir('%s' % directory)
    	print( 'Now Running Dualfoil for cycle: '+N+'\n')
    	os.system('dualfoil_jout_5.2_OCP_2.exe')
#    	with open("dualfoil_jout_5.2.f") as2sourcefile:
#             sourcecode = sourcefile.read()
#    	f2py.compile(sourcecode, modulename='dualfoil')     
#    	print( 'Done.'   )
     
N=2#1969#number of cycles min 2
#R='0.0087' #ohm.m2
R='0.0087' #ohm.m2

SOC_n='0.48' #(mol/m3),

Rtotal='0.030'
write_input_files('0.0087','0.78','0.030')
cycles=[1,430,822,1545,1968,9000]
start_time = time.time()
Cs_neg='25420' #mol/m3
Cs_neg_max= 30555  #mol/m3
run_cycles(str(1)) 
j_para,time=calc_j_side()
for i in range(1,N):        
    write_input_files(str(R),str(SOC_n),str(Rtotal))
    R,Rtotal=degradation(float(R),float(Rtotal),j_para,time)    
    Cs_neg=decrease_concentration(float(Cs_neg),j_para,time)      
    SOC_n=Cs_neg/Cs_neg_max
    print('SOC n', SOC_n)
      
    print('Resistance:'+str(Rtotal)+' Concentration: '+str(Cs_neg))
    if i in cycles:
        print(i)        
        run_cycles(str(i))
        file_in=directory+'/dualfoil5.out'
        file_out=directory+'/Results/dualfoil5_'+str(i)+'.out'
        shutil.copy2(file_in, file_out)
        file_in=directory+'/j_current.out'
        file_out=directory+'/Results/j_current_'+str(i)+'.out'
        shutil.copy2(file_in, file_out)    
        file_in=directory+'/halfcells.out'
        file_out=directory+'/Results/halfcells_'+str(i)+'.out'
        shutil.copy2(file_in, file_out) 
        file_in=directory+'/potential.out'
        file_out=directory+'/Results/potential'+str(i)+'.out'
        shutil.copy2(file_in, file_out)      
     
    if SOC_n <0.58:
        write_input_files(str(R),str(SOC_n),str(Rtotal))

        print('Final cycle is', i)  
        print('Resistance:'+str(Rtotal)+' Concentration: '+str(Cs_neg))        
        run_cycles(str(i))
        file_in=directory+'/dualfoil5.out'
        file_out=directory+'/Results/dualfoil5_'+str(i)+'.out'
        shutil.copy2(file_in, file_out)
        file_in=directory+'/j_current.out'
        file_out=directory+'/Results/j_current_'+str(i)+'.out'
        shutil.copy2(file_in, file_out)    
        file_in=directory+'/halfcells.out'
        file_out=directory+'/Results/halfcells_'+str(i)+'.out'
        shutil.copy2(file_in, file_out) 
        file_in=directory+'/potential.out'
        file_out=directory+'/Results/potential'+str(i)+'.out'
        shutil.copy2(file_in, file_out)      
        break