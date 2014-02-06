#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as pl

# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


#######################################
# function definitions
def tov_RHS(rad,p,rho,m):
    
    # RHS function
    
    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = -ggrav * m * rho / rad**2
        rhs[1] = 4.0*np.pi * rad**2 * rho
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def tov_integrate_FE(rad,dr,p,rho,m):

    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    new = old + dr*tov_RHS(rad,p,rho,m)
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

#######################################

npoints = np.array([1000, 10000, 100000])
radmax = 2.0e8 # 2000 km
euler_radius = np.zeros(3)
euler_mass = np.zeros(3)
RK2_radius = np.zeros(3)
RK2_mass = np.zeros(3)
RK3_radius = np.zeros(3)
RK3_mass = np.zeros(3)
RK4_radius = np.zeros(3)
RK4_mass = np.zeros(3)

for i in np.arange(3):

    print 'using',npoints[i],'points'

    # set up grid
    radius = np.linspace(0.0, radmax, num=npoints[i])
    dr = radius[1]-radius[0]
    
    # set up variables
    press = np.zeros(npoints[i])
    rho   = np.zeros(npoints[i])
    mass  = np.zeros(npoints[i])
    
    # set up central values
    rho[0]   = 1.0e10
    press[0] = polyK * rho[0]**polyG
    mass[0]  = 0.0
    
    # set up termination criterion
    press_min = 1.0e-10 * press[0]
    
    
    
    #######################################
    #  Forward Euler  #
    nsurf = 0
    for n in range(npoints[i]-1):
        
        (press[n+1],mass[n+1]) = tov_integrate_FE(radius[n],
                                                  dr,
                                                  press[n],
                                                  rho[n],mass[n])
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n
    
        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]
    
        # invert the EOS to get density
        rho[n+1] = (press[n+1] / polyK)**(1.0/polyG)
    
    
    euler_radius[i] = radius[nsurf]/1.0e5   # km
    euler_mass[i] = mass[nsurf]/msun        # solar mass
    
    
    
    ########################################
    #  RK2  #  
    nsurf = 0
    for n in np.arange(npoints[i] - 1):
        
        k1 = dr*tov_RHS(radius[n],press[n],rho[n],mass[n])
        
        k2 = dr*tov_RHS(radius[n]+dr/2.0,
                        press[n]+k1[0]/2.0,
                        ((press[n]+k1[0]/2.0) / polyK)**(1.0/polyG), 
                        mass[n] + k1[1]/2.0)
                        
        press[n+1] = press[n] + k2[0]
        mass[n+1] = mass[n] + k2[1]
        rho[n+1] = (press[n+1] / polyK)**(1.0/polyG)
        
        
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n
    
        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]
    
    
    RK2_radius[i] = radius[nsurf]/1.0e5   # km
    RK2_mass[i] = mass[nsurf]/msun        # solar mass
    
    
    
    ##############################################
    #  RK3  #
    nsurf = 0
    for n in np.arange(npoints[i] - 1):
        
        k1 = dr*tov_RHS(radius[n],press[n],rho[n],mass[n])
        
        k2 = dr*tov_RHS(radius[n] + dr/2.0,
                        press[n] + k1[0]/2.0,
                        ( (press[n] + k1[0]/2.0) / polyK )**(1.0/polyG), 
                        mass[n] + k1[1]/2.0)
                
        k3 = dr*tov_RHS(radius[n] + dr,
                        press[n] - k1[0] + 2.0*k2[0],
                        ( (press[n] - k1[0] + 2*k2[0]) / polyK )**(1.0/polyG),
                        mass[n] - k1[1] + 2.0*k2[1]) 
                        
        press[n+1] = press[n] + (1.0/6.0)*(k1[0] + 4.0*k2[0] + k3[0])
        mass[n+1] = mass[n] + (1.0/6.0)*(k1[1] + 4.0*k2[1] + k3[1])
        rho[n+1] = (press[n+1] / polyK)**(1.0/polyG)
        
        
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n
    
        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]
    
    
    RK3_radius[i] = radius[nsurf]/1.0e5   # km
    RK3_mass[i] = mass[nsurf]/msun        # solar mass
    
    
    
    ##############################################
    #  RK4  #
    nsurf = 0
    for n in np.arange(npoints[i] - 1):
        
        k1 = dr*tov_RHS(radius[n],press[n],rho[n],mass[n])
        
        k2 = dr*tov_RHS(radius[n] + dr/2.0,
                        press[n] + k1[0]/2.0,
                        ( (press[n] + k1[0]/2.0) / polyK )**(1.0/polyG), 
                        mass[n] + k1[1]/2.0)
                
        k3 = dr*tov_RHS(radius[n] + dr/2.0,
                        press[n] + k2[0]/2.0,
                        ( (press[n] + k2[0]/2.0) / polyK )**(1.0/polyG),
                        mass[n] + k2[1]/2.0) 
         
        k4 = dr*tov_RHS(radius[n] + dr,
                        press[n] + k3[0],
                        ( (press[n] + k3[0]) / polyK )**(1.0/polyG),
                        mass[n] + k3[1]) 
                        
        press[n+1] = press[n] + (1.0/6.0)*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])
        mass[n+1] = mass[n] + (1.0/6.0)*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])
        rho[n+1] = (press[n+1] / polyK)**(1.0/polyG)
        
        
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n
    
        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]
    
    
    RK4_radius[i] = radius[nsurf]/1.0e5   # km
    RK4_mass[i] = mass[nsurf]/msun        # solar mass
    
    
####################################################################### 

# Look at final values, find convergence rates    

eulerQ_mass = np.abs(euler_mass[2] - euler_mass[1]) \
            / np.abs(euler_mass[1] - euler_mass[0])
eulerQ_radius = np.abs(euler_radius[2] - euler_radius[1]) \
              / np.abs(euler_radius[1] - euler_radius[0])

RK2Q_mass = np.abs(RK2_mass[2] - RK2_mass[1]) \
          / np.abs(RK2_mass[1] - RK2_mass[0])
RK2Q_radius = np.abs(RK2_radius[2] - RK2_radius[1]) \
            / np.abs(RK2_radius[1] - RK2_radius[0])
            
RK3Q_mass = np.abs(RK3_mass[2] - RK3_mass[1]) \
          / np.abs(RK3_mass[1] - RK3_mass[0])
RK3Q_radius = np.abs(RK3_radius[2] - RK3_radius[1]) \
            / np.abs(RK3_radius[1] - RK3_radius[0])
            
RK4Q_mass = np.abs(RK4_mass[2] - RK4_mass[1]) \
          / np.abs(RK4_mass[1] - RK4_mass[0])
RK4Q_radius = np.abs(RK4_radius[2] - RK4_radius[1]) \
            / np.abs(RK4_radius[1] - RK4_radius[0])
     
Q1 = (1.0 / (npoints[2]-1) - 1.0 / (npoints[1]-1)) \
   / (1.0 / (npoints[1]-1) - 1.0 / (npoints[0]-1))       
Q2 = ((1.0 / (npoints[2]-1))**2 - (1.0 / (npoints[1]-1))**2) \
   / ((1.0 / (npoints[1]-1))**2 - (1.0 / (npoints[0]-1))**2)
Q3 = ((1.0 / (npoints[2]-1))**3 - (1.0 / (npoints[1]-1))**3) \
   / ((1.0 / (npoints[1]-1))**3 - (1.0 / (npoints[0]-1))**3)
Q4 = ((1.0 / (npoints[2]-1))**4 - (1.0 / (npoints[1]-1))**4) \
   / ((1.0 / (npoints[1]-1))**4 - (1.0 / (npoints[0]-1))**4)

print 'Euler'
print 'radius:',euler_radius
print 'mass:',euler_mass

print 'RK2'
print 'radius:',RK2_radius
print 'mass:',RK2_mass

print 'RK3'
print 'radius:',RK2_radius
print 'mass:',RK2_mass

print 'RK4'
print 'radius:',RK2_radius
print 'mass:',RK2_mass


print 'CONVERGENCES'

print 'order 1:',Q1
print 'order 2:',Q2
print 'order 3:',Q3
print 'order 4:',Q4

print 'mass convergence'
print 'Euler:',eulerQ_mass
print 'RK2:',RK2Q_mass
print 'RK3:',RK3Q_mass
print 'RK4:',RK4Q_mass

print 'radius convergence'
print 'Euler:',eulerQ_radius
print 'RK2:',RK2Q_radius
print 'RK3:',RK3Q_radius
print 'RK4:',RK4Q_radius


##################################################

# plot RK4 values for npoints = 100000, since they were the last done  

pl.clf()

fig, ax1 = pl.subplots(1, 1, sharex=True)
ax2 = ax1.twinx()

P, = ax1.plot(radius/1.0e5,np.log10(press),'g',linewidth=2)
RHO, = ax1.plot(radius/1.0e5,np.log10(rho),'b',linewidth=2)  
M, = ax2.plot(radius/1.0e5,mass/msun,'r',linewidth=2)

pl.xlim(radius[0]/1.0e5, radius[nsurf]/1.0e5)

ax1.set_xlabel('Radius (km)')
ax1.set_ylabel('log10(Density) (g/cm^3)')
ax2.set_ylabel('Mass (Solar)')
pl.title('Mass, Density, and Pressure vs Radius from RK4 Run')
pl.legend((M,P,RHO), ("Mass (Solar mass)","Pressure (dyne/cm^2)","Density (g/cm^3)"), loc=(0.45,0.4), frameon=False)

pl.plot()

pl.savefig('MPrho.pdf')