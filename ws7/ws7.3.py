#!/usr/bin/env python

import numpy as np
import random
import matplotlib.pyplot as pl

np.random.seed(1)

N = np.round(10**np.linspace(1, 5.5, num=300))

# the array of the MC integral for different sample sizes
I = np.zeros(N.size)

def f(x):
    return x**2 + 1
    
# integration boundaries
a = 2.0
b = 3.0

# constant values used in hit-or-miss method
A1 = (b-a)*f(b)
A2 = (b-a)*(f(b)-f(a))


for i in np.arange(N.size):
    
    print 'doing size',N[i]
    
    # random values of x and y, taken from the uniform
    # distribution from 0 to 1 and altered to fit the range
    # x = [2, 3] and y = [5, 10]
    x = np.random.rand(N[i]) + 2.0
    y = 5.0 * np.random.rand(N[i]) + 5.0

    # Returns an array of true (value 1), when the y coordinate
    # is under the curve in the upper rectangle, and false, (value 0)
    # when the y coordinate is not. Summing the resulting array
    # will then count the number of points under the curve.
    num_under_curve = (f(x) >= y)
    n = float(num_under_curve.sum())
    
    # formula for MC integration using the constants, the counts 
    # calculated above, and the total number of points
    I[i] = A2*(n/N[i]) + A1-A2


# calculate actual integral, get array of errors
actual_I = 22.0 / 3.0
err = actual_I - I

# PLOTTING
#
# error vs log10(N) plot

pl.clf()

errorplot, = pl.plot(np.log10(N),err,'b.')  
    
pl.xlim(min(np.log10(N)),max(np.log10(N)))
pl.ylim(min(err)*1.05,max(err)*1.05)

pl.title('Error of Hit-or-Miss MC Integration vs Sample Size')
pl.xlabel("log10(Sample Size)")
pl.ylabel("Error")

pl.plot()

pl.savefig("MCint_err_plot.pdf")

# log10(error) vs log10(N) plot

#pl.clf()
#
#errorplot, = pl.plot(np.log10(N),np.log10(np.abs(err)),'b.')  
#    
#pl.xlim(min(np.log10(N)),max(np.log10(N)))
#pl.ylim(-4.0,0.05)
#
#pl.title('Error of Hit-or-Miss MC Integration vs Sample Size')
#pl.xlabel("log10(Sample Size)")
#pl.ylabel("log10(Error)")
#
#pl.plot()
#
#pl.savefig("MCint_logerr_plot.pdf")