#!/usr/bin/env python

import numpy as np
import random
import matplotlib.pyplot as pl

np.random.seed(1)

# trying numbers of people from 1 to 30
num_of_people = np.arange(1, 31)

# sample sizes
N = np.round(10**np.linspace(2, 4.75, num=25))

# 2D array of probabilities, with the rows being a specific sample size
# and the columns being a specific number of people
prob = np.zeros((N.size, num_of_people.size))


for i in np.arange(N.size):
    
    print 'doing N =',N[i]
    
    tempprob = np.zeros(num_of_people.size)

    for j in num_of_people:
        
        success = 0.0
        
        for n in np.arange(N[i]):
        
            people = np.random.randint(365, size=j)
            
            # The set() function creates an ordered array
            # without duplicates of whatever you give it.
            # So, if doing that to the array of people reduces
            # the number, there must be duplicates. 
            if len(set(people)) != len(people):
                success += 1.0
        
        tempprob[j-1] = success / N[i]
    
    prob[i] = tempprob



days = np.arange(1.0, 366.0)
actual_prob = np.zeros(num_of_people.size)

# actual way to calculate probabilities from wikipedia
for j in np.arange(num_of_people.size):
    actual_prob[j] = 1.0  - np.prod(days[365-num_of_people[j]:]) / 365.0**num_of_people[j]

# index 22 is 23 people, which is the actual answer,
# which is the error we are interested in.
err = actual_prob[22] - prob[:,22]



# PLOTTING
#
# error vs log10(N) plot

pl.clf()

errorplot, = pl.plot(np.log10(N),err,'b.')  
    
pl.xlim(min(np.log10(N)),max(np.log10(N)))
#pl.ylim(-1.05,1.05)

pl.title('Error of Birthday Problem vs Sample Size')
pl.xlabel("log10(Sample Size)")
pl.ylabel("Error")

pl.plot()

pl.savefig("birth_err_plot.pdf")


# log10(error) vs log10(N) plot

#pl.clf()
#
#errorplot, = pl.plot(np.log10(N),np.log10(np.abs(err)),'b.')  
#    
#pl.xlim(min(np.log10(N)),max(np.log10(N)))
#pl.ylim(min(np.log10(np.abs(err)))-0.05,max(np.log10(np.abs(err)))+0.05)
#
#pl.title('Error of Birthday Problem vs Sample Size')
#pl.xlabel("log10(Sample Size)")
#pl.ylabel("log10(Error)")
#
#pl.plot()
#
#pl.savefig("birth_logerr_plot.pdf") 





