# -*- coding: utf-8 -*-
"""
Transient Radiator Model 
Created on Sun Nov  1 13:29:14 2020

@author: Hakan İbrahim Tol, PhD

References:
[1] Pedersen,Theis Heidmann; Hedegaard, RasmusElbæk; Kristensen, KristianFogh; 
    Gadgaard, Benjamin; Petersen, Steffen. The effect of including hydronic 
    radiator dynamics in model predictive control of space heating. 
    Energy and Buildings. Volume 183. pp 772 - 784.
[2] Brembilla, Christian. Transient model of a panel radiator. In: 14th Conference
    of International Building Performance Simulation Association. Hyderabad, India
    December 7-9, 2015. 
"""

import numpy as np
import matplotlib.pyplot as plt

""" Input for Operational Parameters """
mF=0.01         # [kg/s]    mass flow rate
Ts=55           # [°C]      radiator inlet temperature
Ta=20           # [°C]      air temperature - assumed to be constant

Ti=20.5          # [°C]      initial temperature of the radiator water mass

""" Radiator Properties """
# Panel Radiator Model: Lenhovda MP 25 500
h_rad=0.5       # [m]       radiator height
l_rad=1         # [m]       radiator length
n_r=1.286       # [-]       emprical radiator exponent

# Nominal 
Qn=276          # [W]       nominal heat output rate
LMTDn=30        # [°C]      nominal LMTD

# Mass properties
Mw=3.23         # [kg]      total water mass in the radiator
Mm=10.71        # [kg]      total metal mass of the radiator unit

""" Physical Properties """
Cw=4180         # [J/kgK]   specific heat capacity of water  
Cm=897          # [J/kgK]   specific heat capacity of metal

""" Numerical Parameters """
n=5             # [-]       number of elements
dx=l_rad/n      # [m]       uniform grid spacing (mesh step)
dt=1            # [s]       time step length
Nt=80*60         # [-]       number of time steps

Qi=Qn/n         # [W]       nominal heat output rate per element

Tinit=np.ones(n+1)*Ti # initial temperature distribution along the radiator length

""" Finite Difference Method """
Crad=(Mw*Cw+Mm*Cm)/n # [J/K] Heat capacity per element (water-metal)

D_t=Crad/dt
D_c=Cw*mF
D_q=Qi

# Initialize the solution space
T=np.zeros((n+1,Nt+1))
T[:,0]=Tinit    # initial condition
T[0,:]=Ts       # boundary condition - Dirichlet

for iT in range(1,Nt+1):
    for iS in range(1,n+1):
        T[iS,iT]=D_c/D_t*(T[iS-1,iT-1]-T[iS,iT-1])-D_q/D_t*((T[iS,iT-1]-Ta)/LMTDn)**n_r+T[iS,iT-1]

""" Plotting the Results """
# Outlet temperature from each element
t_plot=np.arange(0,80+dt/60,dt/60) # time between 0 and 80 min

for i in range(n+1):
    if i==0:
        plt.plot(t_plot,T[i,:],linewidth=1,label='Tinlet')
    else:
        plt.plot(t_plot,T[i,:],linewidth=1,label='N='+str(i))


plt.xlabel('Time [min]')
plt.ylabel('Outlet Temperature [°C]')
plt.legend()
plt.show()
