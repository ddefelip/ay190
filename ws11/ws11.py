import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

mpl.clf()

def apply_bcs(x, y):
    # apply boundary conditions
    # sets the boundary grid points of y (0th and last) to be the 
    # last interior points of x (1st and 2nd to last)
    y[0] = x[1]
    y[-1] = x[-2]

def analytic(x, x0, sigma):
    # gaussian function
    return np.exp(-(x-x0)**2 / (2*sigma**2))


# set up the grid so dx = 0.1
x = np.linspace(0, 100, num=1001)

# parameters
dx = x[1]-x[0]
n = len(x)
v = 0.1

##### 
# Change cfl value to observe differing values of v*dt/dx
cfl = 0.5
#####

dt = cfl*(dx/v)

# for initial data
sigma = np.sqrt(15.0) 
x0 = 30.0

#set up initial conditions
yupwind = analytic(x, x0, sigma)
yftcs = analytic(x, x0, sigma)
ylax = analytic(x, x0, sigma)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,yupwind,'x-') # numerical data - upwind
mpl.plot(x,yftcs,'x-') # numerical data - FTCS
mpl.plot(x,ylax,'x-') # numerical data - Lax-Friedrich
mpl.plot(x,analytic(x, x0, sigma),'r-') # analytic data
mpl.show()

yold_ftcs = yftcs
yold_up = yupwind
yold_lax = ylax
ntmax = 1000
t = np.linspace(0., ntmax*dt, num=ntmax+1)
err = np.zeros(ntmax+1)
err[0] = 0.
for it in range(ntmax):
    # save previous data
    yold_ftcs = yftcs
    yold_up = yupwind
    yold_lax = ylax

    # get new data for UPWIND scheme, unstable FTCS, and Lax-Friedrich
    for j in range(1, n-1):
        yupwind[j] = yold_up[j] - \
        (v*dt/dx)*(yold_up[j] - yold_up[j-1])
        yftcs[j] = yold_ftcs[j] - \
        (v*dt/(2.*dx))*(yold_ftcs[j+1] - yold_ftcs[j-1])
        ylax[j] = 0.5*(yold_lax[j+1] + yold_lax[j-1]) - \
        (v*dt/(2.*dx))*(yold_lax[j+1] - yold_lax[j-1])
    # after updates, apply boundary conditions
    apply_bcs(yold_up, yupwind)
    apply_bcs(yold_ftcs, yftcs)
    apply_bcs(yold_lax, ylax)

    # get analytic result for time t
    yana = analytic(x, x0+(t[it]*v), sigma)
    # compute error estimage

    err[it] = max(yupwind) - max(yana)
    print "it = ",it,err[it]
    mpl.clf()
    # plot numerical results
    mpl.plot(x,yupwind,'x-')
    mpl.plot(x,yftcs,'x-')
    mpl.plot(x,ylax,'x-')
    # plot analytic results
    mpl.plot(x,yana,'r-')
    mpl.draw()

mpl.show()
mpl.clf()
sigma1, = mpl.plot(t,abs(err),'r')
mpl.title('Upwind Scheme Error vs time')
mpl.xlabel('Time')
mpl.ylabel('Difference between peaks')
