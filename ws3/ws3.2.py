#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from scipy import special as sp

pl.clf()

kT = 20.0 #MeV     
hbar = 1.05457266e-27*6.24150934e5 #MeV * s
c = 2.99792458e10 #cm*s^(-1)


print 'Gauss-Laguerre Quadrature'

I = np.zeros(10)

for n in np.arange(2, 21, 2):
    # trying different numbers of nodes, from 2 to 20, to confirm convergence
    
    
    [laguerre_roots,laguerre_weights] = sp.l_roots(n, 0)
        # line of code from worksheet
        # provudes easily accessible laguerre weights and roots
        # with n nodes.
    
    laguerre_array = laguerre_weights * laguerre_roots**2  / (1 + np.exp(-laguerre_roots))
        # array of weights times the function f(x_i), where x_u are the roots,
        # and f(x)*x^c*exp(-x) = x^2 / (exp(x) + 1). (No poles, so set c = 0).
    
    I[n/2.0 - 1] = (8*np.pi*(kT)**3)/(2*np.pi*hbar*c)**3*np.sum(laguerre_array)
        # now, simply sum the array!
    
    print 'nodes =',n,', number density =',I[n/2.0 - 1]





print 'Gauss-Legendre Quadrature'

[legendre_roots,legendre_weights] = sp.p_roots(20, 0)

x = np.linspace(0, 200, num=41) / 20.0
    # Energy bins of 5 MeV correspond to x bins of E/20 = 0.25
    # because kT = 20 MeV and x = E/kT.
    # Going out to E = 200 MeV to make sure all relevant terms 
    # are covered.

dndE = np.zeros(40)

for i in np.arange(x.size - 1):
    
    
    legendre_roots_transformed = legendre_roots/8.0 + (x[i]+x[i+1])/2.0 
        # need to do this because Gauss-Legendre quadrature needs 
        # to be done from -1 to 1, and this transformation allows that.
    
    legendre_array = legendre_weights * (legendre_roots_transformed)**2 / (np.exp(legendre_roots_transformed)+1)
        # no exp(x) term in the numerator here, because W(x) = 1 for 
        # Gauss-legendre quadrature.
    
    Isub = (8*np.pi*(kT)**3)/(2*np.pi*hbar*c)**3 * (0.125) * np.sum(legendre_array)
        # the 1/8 comes from the transformation. Since the x bins are all the 
        # same size, (x[i+1] - x[i])/2 = 1/8 = 0.125 always.
    
    print 'fraction of bin E = [%g, %g] is %f' % (20*x[i], 20*x[i+1], Isub / I[9])
        
    
    dndE[i] = Isub / 5.0
    


print 5.0*np.sum(dndE)
print 'relative error from part a =',(5.0*np.sum(dndE) - I[9]) / I[9]
    # make sure that this is close to the most accurate integral 
    # calculated in part a, because the infinite sum 
    # should be identical
  
  
#Plot!
energy, = pl.plot(x[:-1]*20.0,dndE,"r",linewidth=2)

pl.xlim(min(x[:-1]*20.0),max(x[:-1]*20.0))
pl.ylim(0, 1.05*max(dndE))

pl.title('Energy Distribution dn/dE of Electrons')
pl.xlabel("E (MeV)")
pl.ylabel("dn/dE (cm^3 MeV^(-1))")

pl.savefig("gauss-legendre.pdf")

pl.plot()