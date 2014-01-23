#!/usr/bin/env python

import numpy as np

a = 1.0          # AU, =1.496e6 km
T = 365.25635    # days
e = 0.0167
## switch values of the eccentricity 
## by commenting above value, and uncommenting below
#e = 0.99999 
##
omega = 2*np.pi / T  
b = a**2 * (1-e**2)
epsilon = 1.0e-10    # allowable error
    

# function we are trying to find roots of.
def f(x):
    return x - omega*t - e*np.sin(x)

   
E = np.zeros(1000) 
    # grid to put values of E calculated via root finding technique


## TIME t ##
t = 273.0

## Initial guesses ##
E[0] = omega*t
E[1] = (0.9)*E[0]


# Secant Method
for i in np.arange(1, 1000):
    E[i+1] = E[i] - f(E[i])*( E[i] - E[i-1] )/( f(E[i]) - f(E[i-1]) )
    error = np.abs( (E[i+1] - E[i]) / E[i] )
    print 'E=',E[i+1]
    print 'error is',error
    if error < epsilon:
        print 'calculation took',i,'iterations'
        print 'x =',a*np.cos(E[i+1]),'and y =',b*np.sin(E[i+1])
        break
        

if e == 0.99999:
# Newton's Method - may reduce number of iterations for e = 0.99999
    def fprime(x):
        return 1 - e*np.cos(x)
    
    for i in np.arange(1, 1000):
        E[i+1] = E[i] - f(E[i])/fprime(E[i])
        error = np.abs( (E[i+1] - E[i]) / E[i] )
        print 'E=',E[i+1]
        print 'error is',error
        if error < epsilon:
            print 'calculation took',i,'iterations'
            print 'x =',a*np.cos(E[i+1]),'and y =',b*np.sin(E[i+1])
            break


