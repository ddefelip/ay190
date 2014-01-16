#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

#clears any previous figures on the plot
pl.clf()

def f(x):
    return x**3 - 5*x**2 + x
    
def fprime(x):
    return 3*x**2 - 10*x + 1   
   
h1 = 0.01
h2 = 0.005

h1forward = (f(np.arange(-2.0+h1, 6.0+h1, h1)) - f(np.arange(-2.0, 6.0, h1)))/h1
h1central = (f(np.arange(-2.0+h1, 6.0+h1, h1)) - f(np.arange(-2.0-h1, 6.0-h1, h1)))/(2*h1)
h2forward = (f(np.arange(-2.0+h2, 6.0+h2, h2)) - f(np.arange(-2.0, 6.0, h2)))/h2
h2central = (f(np.arange(-2.0+h2, 6.0+h2, h2)) - f(np.arange(-2.0-h2, 6.0-h2, h2)))/(2*h2)

x1 = np.arange(-2.0, 6.0, h1)
x2 = np.arange(-2.0, 6.0, h2)
#errors of the different numerical derivatives
y1forward = h1forward-fprime(np.arange(-2.0, 6.0, h1))
y1central = h1central-fprime(np.arange(-2.0, 6.0, h1))
y2forward = h2forward-fprime(np.arange(-2.0, 6.0, h2))
y2central = h2central-fprime(np.arange(-2.0, 6.0, h2))


p1forward, = pl.plot(x1,y1forward,"r",linewidth=2)
p2forward, = pl.plot(x2,y2forward,"b",linewidth=2)
#p1central, = pl.plot(x1,y1central,"r",linewidth=2)
#p2central, = pl.plot(x2,y2central,"b",linewidth=2)

pl.xlim(min(x1),max(x1))
pl.ylim(min(y1forward*1.05),max(y1forward*1.05))
#pl.ylim(0, 1.2e-4)

pl.title('Forward Differencing')
#pl.title('Central Differencing')
pl.xlabel("x")
pl.ylabel("f(x;h)-f(x)")

pl.legend((p1forward,p2forward), ("h1=0.01","h2=0.005"), loc=(0.05,0.75), frameon=False)
#pl.legend((p1central,p2central), ("h1=0.01","h2=0.005"), loc=(0.05,0.45), frameon=False)

pl.savefig("forward_differencing.pdf")
#pl.savefig("central_differencing.pdf")
pl.plot()

