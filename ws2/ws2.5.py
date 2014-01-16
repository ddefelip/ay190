#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as inter

pl.clf()

time = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
mag = np.array([0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561, 0.468, 0.302])
npoints=9
xpoints = np.linspace(0.0, 1.0, num=1001)
    # discretize x = [0, 1]



#PIECEWISE CUBIC HERMITE INTERPOLATION

fprime_forward = np.array([(mag[1]-mag[0]) / (time[1]-time[0])])
fprime_central = (mag[2:]-mag[:-2]) / (time[2:]-time[:-2])
fprime_backward = np.array([(mag[8]-mag[7]) / (time[8]-time[7])])
fprime = np.concatenate((fprime_forward, fprime_central, fprime_backward))


z = (xpoints[:200] - time[0]) / (time[1] - time[0])
for i in np.arange(1, 7):
    z = np.concatenate((z,(xpoints[100*(i+1):100*(i+1)+100] - time[i]) / (time[i+1] - time[i]) ))
z = np.concatenate((z, (xpoints[800:] - time[7]) / (time[8] - time[7])))

psi0 = 2*z**3 - 3*z**2 + 1
psi0_1z = 2*(1-z)**3 - 3*(1-z)**2 + 1
psi1 = z**3 - 2*z**2 + z
psi1_1z = (1-z)**3 - 2*(1-z)**2 + (1-z)

H = ( mag[0]*psi0[:200] + mag[1]*psi0_1z[:200] + 
      fprime[0]*(time[1]-time[0])*psi1[:200] - fprime[1]*(time[1]-time[0])*psi1_1z[:200] )

for i in np.arange(1, 7):
    H = np.concatenate((H, ( mag[i]*psi0[100*(i+1):100*(i+1)+100] + 
                             mag[i+1]*psi0_1z[100*(i+1):100*(i+1)+100] + 
                             fprime[i]*(time[i+1]-time[i])*psi1[100*(i+1):100*(i+1)+100] - 
                             fprime[i+1]*(time[i+1]-time[i])*psi1_1z[100*(i+1):100*(i+1)+100] ) ))
                            
H = np.concatenate((H, (mag[7]*psi0[800:] + mag[8]*psi0_1z[800:] + 
                        fprime[7]*(time[8]-time[7])*psi1[800:] - 
                        fprime[8]*(time[8]-time[7])*psi1_1z[800:] ) ))




# CUBIC SPLINE INTERPOLATION

#methods in scipy.interpolate
sp = inter.splrep(time, mag, s=0)
S = inter.splev(xpoints, sp)


#PLOTTING

hermite, = pl.plot(xpoints,H,"r",linewidth=2)
spline, = pl.plot(xpoints,S,"g",linewidth=2)
data, = pl.plot(time,mag,"bo")
pl.plot()

pl.xlim(min(time),max(time))
pl.ylim(0, max(1.05*S))

pl.title('Cubic Spline and Hermite Interpolations')
pl.xlabel("time (days)")
pl.ylabel("magnitude")

pl.legend((hermite,spline,data), ("Cubic hermite p(t)","Cubic spline p(t)","data"),loc=(0.05,0.75),frameon=False)
pl.savefig("hermite_spline_interpolation.pdf")

pl.plot()