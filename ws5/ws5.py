#!/usr/bin/env python

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as pl


pl.clf()


# Part a

# read data, get log variables
data = ascii.read("m_sigma_table.dat",readme="m_sigma_ReadMe.dat")  
logsigma = np.array(np.log10(data["sigma*"]))
logM = np.array(data["logM"])

points, = pl.plot(logsigma,logM,"r.")

pl.xlim(min(logsigma)/1.1,max(logsigma)*1.1)
#pl.ylim(min(logM)/1.05, max(logM)*1.05)

pl.title('Black Hole Mass vs Stellar Velocity Dispersion')
pl.xlabel("log(sigma)")
pl.ylabel("log(M)")

#pl.savefig("ws5-a.pdf")




print 'Part b'

# ignoring errors means S in III.6.8 in the notes is logM.size.

# calculate fit parameters, where y = a1 + a2*x, using formula III.6.9.
a1 = ( np.sum(logM)*np.sum(logsigma**2) - np.sum(logsigma)*np.sum(logsigma*logM) )\
     / (logM.size*np.sum(logsigma**2) - np.sum(logsigma)**2)
     
a2 = ( logM.size*np.sum(logsigma*logM) - np.sum(logsigma)*np.sum(logM) )\
     / (logM.size*np.sum(logsigma**2) - np.sum(logsigma)**2)

# calculate uncertainties on the fit parameters using formula III.6.12.
sigma_a1 = np.sqrt(np.sum(logsigma**2) / (logM.size*np.sum(logsigma**2) - np.sum(logsigma)**2))
sigma_a2 = np.sqrt(logM.size / (logM.size*np.sum(logsigma**2) - np.sum(logsigma)**2))

# plot data points (no errors) and regression line together
xpoints = np.linspace(1.4, 2.6, num=1000)
regression, = pl.plot(xpoints,a1+a2*xpoints,'b',linewidth=2)


#pl.legend((points,regression), ("data points","regression line"), loc=(0.05,0.75), frameon=False)
#pl.plot()
#pl.savefig("ws5-b.pdf")

# print the fir parameters
print 'Fit parameters:'
print 'a1 =',a1,'+/-',sigma_a1
print 'a2 =',a2,'+/-',sigma_a2






print 'Part c'

# approximate the derivative dy/dx at all points by using the slope 
# of the fit line with no uncertainties, given by a2, for use in III.6.10.
dydx = a2

# correctly transform the errors on sigma to the errors on log(sigma):
# e_logsigma = 1/ln(10) * e_sigma/sigma
e_logsigma = (1/np.log(10.0)) * np.array(data["e_sigma*"]) / np.array(data["sigma*"])
e_logM = np.array(data["e_logM"])

# add errors in quadrature to get total error to use in formulas
# for fit parameters
e_total = np.sqrt(e_logM**2 + dydx**2 * e_logsigma**2)

Sigma_x = np.sum(logsigma/e_total**2)
Sigma_y = np.sum(logM/e_total**2)
Sigma_x2 = np.sum((logsigma/e_total)**2)
Sigma_y2 = np.sum((logM/e_total)**2)
Sigma_xy = np.sum((logsigma*logM)/e_total**2)
S = np.sum(1.0/e_total**2)
    # note, S is not necessarily equal to 1 here

# calculate new fit parameters, and errors on those parameters
a1_new = (Sigma_y*Sigma_x2 - Sigma_x*Sigma_xy) / (S*Sigma_x2 - Sigma_x**2)
a2_new = (S*Sigma_xy - Sigma_x*Sigma_y) / (S*Sigma_x2 - Sigma_x**2)

sigma_a1_new = np.sqrt(Sigma_x2 / (S*Sigma_x2 - Sigma_x**2))
sigma_a2_new = np.sqrt(S / (S*Sigma_x2 - Sigma_x**2))

regression_new, = pl.plot(xpoints,a1_new+a2_new*xpoints,'g',linewidth=2)
pl.errorbar(logsigma,logM, yerr=e_logM, xerr=e_logsigma, fmt='ro')

pl.legend((points,regression,regression_new), ("data points with error bars","regression without errors","regression with errors"), loc=(0.03,0.78), frameon=False)
pl.plot()
pl.savefig("ws5-c.pdf")

# print the fit parameters
print 'Fit parameters:'
print 'a1 =',a1_new,'+/-',sigma_a1_new
print 'a2 =',a2_new,'+/-',sigma_a2_new

