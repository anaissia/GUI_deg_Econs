import matlab.engine
import numpy as np
import scipy.integrate
import math
eng = matlab.engine.start_matlab()
def degradation():
#list of input
    Rf=8.7*10**-3
    Cnom=1.67
    Tref = 25 + 273.15
    Rc = 20.0*10**-4
    thick1 = 75.0*10**-6
    thick2 = 25.0*10**-6
    thick3 = 740.*10**-6
    As = 0.087
    Rs1 = 2.0*10**-6
    Rs3= 2.0*10**-6
    eps1s = 0.51
    eps3s = 0.48
    cs1max = 30550.0
    cs3max = 51555.0
    x1soc0 = 0.0115
    x1soc1 = 0.83
    y3soc0 = 0.95
    y3soc1 = 0.48
    Ds1ref = 3.8*10**-14
    Ds3ref = 1.0*10**-13
    EaDs1 = 35.0*10**3
    EaDs3 = 29.0*10**3
    k1ref = 1.0*10**-6
    k3ref = 10.*10**-6
    Eak1=20.0e3
    Eak3=58.0e3
    ceavg = 1.0*10**3   
    Cs_N= x1soc0*cs1max	
    Cs_N_disc= x1soc1*cs1max	
    cycle=2	
    N=[1,600,700]	
    for c in range(1, cycle):

	    print('Model running for cycle:',c,'/n')	
	#run Matlab model
	    result=eng.EXAMPLE_constant_charge(Cnom,\
	    Tref,\
	    Rc,\
	    thick1,\
	    thick2,\
	    thick3,\
	    As ,\
	    Rs1,\
	    Rs3,\
	    eps1s,\
	    eps3s,\
	    cs1max,\
	    cs3max,\
	    x1soc0 ,\
	    x1soc1,\
	    y3soc0 ,\
	    y3soc1 ,\
	    Ds1ref,\
	    Ds3ref ,\
	    EaDs1,\
	    EaDs3 ,\
	    k1ref,\
	    k3ref,\
	    Eak1,\
	    Eak3,\
	    ceavg)
	#Run degradation
	    j_p=result['j1']._data.tolist()
	    phase_diff_p=result['phase_diff']._data.tolist()
	    t_p=result['time']._data.tolist()	
	    temp_p=result['temperature']._data.tolist()	
	    j_para=calc_j_side(phase_diff_p,j_p,Rf,temp_p)
	    Rf=increase_resistance(Rf,j_para,t_p)
	    Rc=Rc+Rf	
	    Cs_N,diff_C=decrease_concentration(Cs_N,j_para,t_p)
	    #print(diff_C,"degradation C",Rf, "resistance deg") 
	    print(j_para)		
	    Cs_N=float(Cs_N)
	    Rf=float(Rf)
	    Rc=float(Rc)
	    x1soc0=float(Cs_N/cs1max)
	    x1soc1=float((Cs_N_disc-diff_C)/cs1max)

	


def calc_j_side(phase,j_p,Rf,temp):
  
    Upara=0.38 #V
    j0=0.8*10**-7
    F=96485.3329 #faraday's constant
    R=8.3144598 # J/mol k
    alpha_c=0.5
    overpotential=[]
    jpara=[]

    for i in range(0,len(phase)):
        overpotential.append(phase[i]-Upara-Rf*j_p[i])
        jpara.append(-j0*math.exp(alpha_c*F*overpotential[i]/(R*temp[i])))    
    return jpara

def increase_resistance(Rf_N,j_para,time):

      #get the rate of SEI growth
	#get the increase in resistance
	
      # Ning et al. Cyclce life modeling of lithiunm ion batteries
  	rho=2.1*1000 #kg/m3

  	kappa=3.79*(10**-7) #S/m
      
  	M=0.1#kg/mol
  	F=96485 #C/mol (C=A.s)

  	thickness_SEI=[]
	

               
  	thickness_SEI=abs(scipy.integrate.simps(j_para,time)*M/(F*rho)) # because integral

  	Rcycle=thickness_SEI/kappa

  	R_N1= Rf_N+Rcycle #ohm.m2

      
  	return R_N1
   
def decrease_concentration(Cs_N,j_para,time):
  	eps_e=0.440 #electrolyte volume frac
  	eps_fl=0.07 #conductive filler volume frac 
  	rs=2*10**-6 #radius particule m
  	eps_s=(1-eps_e-eps_fl) 
  	a_neg=3*eps_s/rs
  	F=96485 #C/mol (C=A.s) Faraday
      
  	# Calculate volume-averaged capacity lost due to parasitic reaction (C/m3)
        
  	Q_tot=abs(scipy.integrate.simps(j_para,time)*a_neg)
  	Cs_N1=Cs_N-(Q_tot)/(eps_s*F)
  	diff_C=  (Q_tot)/(eps_s*F)

  	return Cs_N1,diff_C   		



degradation()
