#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

pl.clf()

time = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
mag = np.array([0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561, 0.468, 0.302])
npoints=9
xpoints = np.linspace(0.0, 1.0, num=1001)
    # discretize x = [0, 1]

#LAGRANGIAN INTERPOLATION
#Lagrangian Interpolation Product function
def L(x, num_of_points, jth_point):
    L = 1.0
    for k in np.arange(num_of_points):
        if k != jth_point:
            L *= (x - time[k]) / (time[jth_point] - time[k])
    return L

Lx = np.array([L(xpoints, npoints, j) for j in np.arange(npoints)]) 
    #apply L to all of the points
Lx = np.swapaxes(Lx, 0, 1)
    #change the format of the resulting array so that np.sum works correctly
p = np.sum(mag*Lx, axis=1)
    #sum along horizontal axes


#PIECEWISE LINEAR

fprime = (mag[1:] - mag[:-1]) / (time[1:] - time[:-1])
    #define coefficient of (x-x_i) using numpy arrays
p_linear = mag[0] + fprime[0]*(xpoints[:200] - time[0])
    #first and last intervals are larger than the others,
    #so deal with them separately
for i in np.arange(2, 8):
    p_linear = np.concatenate((p_linear, mag[i-1]+fprime[i-1]*(xpoints[100*i:100*i+100]-time[i-1])))
p_linear = np.concatenate((p_linear, mag[7] + fprime[7]*(xpoints[800:] - time[7])))


#PIECEWISE QUADRATIC

fa = mag[:-2] / ((time[:-2] - time[1:8])*(time[:-2] - time[2:]))
fb = mag[1:8] / ((time[1:8] - time[:-2])*(time[1:8] - time[2:]))
fc = mag[2:] / ((time[2:] - time[:-2])*(time[2:] - time[1:8]))
    #coefficients defined ahead of time using existing arrays

p_quad = ( fa[0]*(xpoints[:200]-time[1])*(xpoints[:200]-time[2]) + 
           fb[0]*(xpoints[:200]-time[0])*(xpoints[:200]-time[2]) +
           fc[0]*(xpoints[:200]-time[0])*(xpoints[:200]-time[1]) )
    #first and last intervals are larger, deal with them separately           
for i in np.arange(1, 7):
    p_quad = np.concatenate((p_quad, ( fa[i]*(xpoints[100*(i+1):100*(i+1)+100]-time[i+1])*(xpoints[100*(i+1):100*(i+1)+100]-time[i+2]) + 
                                       fb[i]*(xpoints[100*(i+1):100*(i+1)+100]-time[i])*(xpoints[100*(i+1):100*(i+1)+100]-time[i+2]) +
                                       fc[i]*(xpoints[100*(i+1):100*(i+1)+100]-time[i])*(xpoints[100*(i+1):100*(i+1)+100]-time[i+1]) ) ))
                                       
p_quad = np.concatenate((p_quad, ( fa[6]*(xpoints[800:]-time[7])*(xpoints[800:]-time[8]) + 
                                   fb[6]*(xpoints[800:]-time[6])*(xpoints[800:]-time[8]) +
                                   fc[6]*(xpoints[800:]-time[6])*(xpoints[800:]-time[7]) ) ))


#PLOTTING

lagrange, = pl.plot(xpoints,p,"r",linewidth=2)
linear, = pl.plot(xpoints,p_linear,"g",linewidth=2)
quadratic, = pl.plot(xpoints,p_quad,"c",linewidth=2)
data, = pl.plot(time,mag,"bo")
pl.plot()

pl.xlim(min(time),max(time))
pl.ylim(0, 3)

pl.title('Linear, Quadratic, and Lagrange Interpolations')
pl.xlabel("time (days)")
pl.ylabel("magnitude")

pl.legend((lagrange,linear,quadratic,data), ("Lagrange p(t)","Piecewise linear p(t)","Piecewise quadratic p(t)","data"), loc=(0.05,0.65), frameon=False)

pl.savefig("linear_quadratic_lagrange_interpolation.pdf")

pl.plot()
