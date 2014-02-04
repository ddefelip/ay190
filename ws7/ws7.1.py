#!/usr/bin/env python

import numpy as np
import random
import matplotlib.pyplot as pl

# seed the random number generator
np.random.seed(1)

# using this line for values of N gets 250 logarithmically evenely spaced
# numbers from 10 to about 300,000, which is about what it takes before 
# the computation gets slow
N = np.round(10**np.linspace(1, 5.5, num=250))

# counter array
n = np.zeros(250)


for i in np.arange(N.size):
    
    # get random x and y variables
    # N[i] sets the number of points
    x = np.random.rand(N[i])
    y = np.random.rand(N[i])
    
    dist = x**2 + y**2
        # technically distance squared, but we choose the 
        # radius of the circle to be 1 so we don't need to 
        # take the square root

    for point in dist:
        
        if point <= 1.0:
            n[i] += 1
                # increment counter
            
    print 'did N = ',N[i]


# value of pi according to MC method
pi = 4.0*n/N

# error from actual value
err = np.pi - pi


# PLOTTING
#
# error vs log10(N) plot

pl.clf()

errorplot, = pl.plot(np.log10(N),err,'b.')  
    
pl.xlim(min(np.log10(N)),max(np.log10(N)))
pl.ylim(-1.05,1.05)

pl.title('Error of Monte Carlo Experiment vs Sample Size')
pl.xlabel("log10(Sample Size)")
pl.ylabel("Error")

pl.plot()

pl.savefig("pi_err_plot.pdf")


# log10(error) vs log10(N) plot

#pl.clf()
#
#errorplot, = pl.plot(np.log10(N),np.log10(np.abs(err)),'b.')  
#    
#pl.xlim(min(np.log10(N)),max(np.log10(N)))
#pl.ylim(min(np.log10(np.abs(err)))-0.05,max(np.log10(np.abs(err)))+0.05)
#
#pl.title('Error of Monte Carlo experiment vs Sample Size')
#pl.xlabel("log10(Sample Size)")
#pl.ylabel("log10(Error)")
#
#pl.plot()
#
#pl.savefig("pi_logerr_plot.pdf")
