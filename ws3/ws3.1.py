#!/usr/bin/env python

import numpy as np

h1 = np.pi / 10.0
h2 = np.pi / 20.0 #smaller step size

print 'Computing the integral of sin(x) from 0 to pi'

I = 2.0
print 'actual answer is I =',I 

#Midpoint rule
Q1 = 0.0
ai1 = 0.0
bi1 = h1
for i in np.arange(np.pi/h1):
    Q1 += h1*np.sin((ai1+bi1)/2)
    ai1 += h1
    bi1 += h1
print 'Midpoint rule, step size =',h1,'gives I =',Q1

Q2 = 0.0
ai2 = 0.0
bi2 = h2
for i in np.arange(np.pi/h2):
    Q2 += h2*np.sin((ai2+bi2)/2)
    ai2 += h2
    bi2 += h2
print 'Midpoint rule, step size =',h2,'gives I =',Q2

#relative error between different step sizes
print 'relative error =',(Q2-I)/(Q1-I)


#Trapezoidal rules
Q1 = 0.0
ai1 = 0.0
bi1 = h1
for i in np.arange(np.pi/h1):
    Q1 += 0.5*h1*(np.sin(ai1)+np.sin(bi1))
    ai1 += h1
    bi1 += h1
print 'Trapezoid rule, step size =',h1,'gives I =',Q1

Q2 = 0.0
ai2 = 0.0
bi2 = h2
for i in np.arange(np.pi/h2):
    Q2 += 0.5*h2*(np.sin(ai2)+np.sin(bi2))
    ai2 += h2
    bi2 += h2
print 'Trapezoid rule, step size =',h2,'gives I =',Q2

#relative error between different step sizes
print 'relative error =',(Q2-I)/(Q1-I)


#Simpson's rule
Q1 = 0.0
ai1 = 0.0
bi1 = h1
for i in np.arange(np.pi/h1):
    Q1 += (h1/6.0)*(np.sin(ai1)+4*np.sin((ai1+bi1)/2)+np.sin(bi1))
    ai1 += h1
    bi1 += h1
print 'Simpson\'s rule, step size =',h1,'gives I =',Q1

Q2 = 0.0
ai2 = 0.0
bi2 = h2
for i in np.arange(np.pi/h2):
    Q2 += (h2/6.0)*(np.sin(ai2)+4*np.sin((ai2+bi2)/2)+np.sin(bi2))
    ai2 += h2
    bi2 += h2
print 'Simpson\'s rule, step size =',h2,'gives I =',Q2

#relative error between different step sizes
print 'relative error =',(Q2-I)/(Q1-I)





print 'Computing the integral of x*sin(x) from 0 to pi'

I = np.pi # actual answer
print 'actual answer is I =',I 
#Midpoint rule
Q1 = 0.0
ai1 = 0.0
bi1 = h1
for i in np.arange(np.pi/h1):
    Q1 += 0.5*h1*(ai1+bi1) * np.sin((ai1+bi1)/2)
    ai1 += h1
    bi1 += h1
print 'Midpoint rule, step size =',h1,'gives I =',Q1

Q2 = 0.0
ai2 = 0.0
bi2 = h2
for i in np.arange(np.pi/h2):
    Q2 += 0.5*h2*(ai2+bi2) * np.sin((ai2+bi2)/2)
    ai2 += h2
    bi2 += h2
print 'Midpoint rule, step size =',h2,'gives I =',Q2

#relative error between different step sizes
print 'relative error =',(Q2-I)/(Q1-I)


#Trapezoidal rules
Q1 = 0.0
ai1 = 0.0
bi1 = h1
for i in np.arange(np.pi/h1):
    Q1 += 0.5*h1*(ai1*np.sin(ai1) + bi1*np.sin(bi1))
    ai1 += h1
    bi1 += h1
print 'Trapezoid rule, step size =',h1,'gives I =',Q1

Q2 = 0.0
ai2 = 0.0
bi2 = h2
for i in np.arange(np.pi/h2):
    Q2 += 0.5*h2*(ai2*np.sin(ai2) + bi2*np.sin(bi2))
    ai2 += h2
    bi2 += h2
print 'Trapezoid rule, step size =',h2,'gives I =',Q2

#relative error between different step sizes
print 'relative error =',(Q2-I)/(Q1-I)


#Simpson's rule
Q1 = 0.0
ai1 = 0.0
bi1 = h1
for i in np.arange(np.pi/h1):
    Q1 += (h1/6.0)*(ai1*np.sin(ai1) + 2*(ai1+bi1)*np.sin((ai1+bi1)/2) + bi1*np.sin(bi1))
    ai1 += h1
    bi1 += h1
print 'Simpson\'s rule, step size =',h1,'gives I =',Q1

Q2 = 0.0
ai2 = 0.0
bi2 = h2
for i in np.arange(np.pi/h2):
    Q2 += (h2/6.0)*(ai2*np.sin(ai2) + 2*(ai2+bi2)*np.sin((ai2+bi2)/2) + bi2*np.sin(bi2))
    ai2 += h2
    bi2 += h2
print 'Simpson\'s rule, step size =',h2,'gives I =',Q2

#relative error between different step sizes
print 'relative error =',(Q2-I)/(Q1-I)

