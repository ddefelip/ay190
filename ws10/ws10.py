#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as pl

# set up grid
xmin = 0.0
xmax = 1.0
npoints = 1000
# set up grid
x = np.linspace(0.0, 1.0, num=npoints)
# dx based on x[1] and x[0]
dx = x[1] - x[0]

# boundary values
A = 0.0 # inner boundary
B = 0.1 # outer boundary

def calc_rhs(u,xx):
    # rhs routine
    # rhs[0] is rhs for y
    # rhs[1] is rhs for u
    rhs = np.zeros(2)
    rhs[0] = u
    rhs[1] = 12.0*xx - 4.0

    return rhs

def integrate_FE(z,x):
    # forward-Euler integrator
    
    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    yy = np.zeros((npoints,2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    for i in range(npoints-1):
        yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i])

    return yy
    
    
def integrate_RK2(z,x):
    
    yy = np.zeros((npoints,2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    for i in range(npoints-1):
        
        k1 = dx*calc_rhs(yy[i,1],x[i])
        
        k2 = dx*calc_rhs(yy[i,1] + k1[1]/2.0,
                         x[i] + dx/2.0)
                        
        yy[i+1,:] = yy[i,:] + k2 

    return yy    


print 'Euler'
# get initial guess for derivative
z0 = -1100000.0
z1 = 10000000.0
yy0 = integrate_FE(z0,x)
yy1 = integrate_FE(z1,x)
phi0 = yy0[npoints-1,0] - B
phi1 = yy1[npoints-1,0] - B
dphidz = (phi1 - phi0) / (z1 - z0) # dphi/dz

i = 0
itmax = 100
err = 1.0e99
criterion = 1.0e-12

z0 = z1
phi0 = phi1
while (err > 1.0e-12 and i < itmax):
    z1 = z0 - phi0 / dphidz # secand update
    yy = integrate_FE(z1,x)
    phi1 = yy[npoints-1,0] - B
    dphidz = (phi1 - phi0) / (z1 - z0) # dphi/dz numerical
    err = np.abs(phi1) # your error measure
    z0 = z1
    phi0 = phi1
    i = i+1

    print i,z1,phi1


print 'RK2'
# get initial guess for derivative
z0 = -1100000.0
z1 = 10000000.0
yy0 = integrate_RK2(z0,x)
yy1 = integrate_RK2(z1,x)
phi0 = yy0[npoints-1,0] - B
phi1 = yy1[npoints-1,0] - B
dphidz = (phi1 - phi0) / (z1 - z0) # dphi/dz

i = 0
itmax = 100
err = 1.0e99
criterion = 1.0e-12

z0 = z1
phi0 = phi1
while (err > 1.0e-12 and i < itmax):
    z1 = z0 - phi0 / dphidz # secand update
    yy_RK2 = integrate_RK2(z1,x)
    phi1 = yy[npoints-1,0] - B
    dphidz = (phi1 - phi0) / (z1 - z0) # dphi/dz numerical
    err = np.abs(phi1) # your error measure
    z0 = z1
    phi0 = phi1
    i = i+1

    print i,z1,phi1


def soln(xx):
    return 2.0*xx**3 - 2*xx**2 + 0.1*xx

# store values for testing convergence
if npoints == 10:
    euler10 = yy[:,0]
    RK210 = yy_RK2[:,0]
    soln10 = soln(x)

if npoints == 1000:
    euler1000 = yy[:,0]
    RK21000 = yy_RK2[:,0]
    soln1000 = soln(x)


# uncomment these next lines to calculate Q factor
# after running the code for npoints = 10 and 1000 separately:

#Q_euler = max(np.abs(euler1000 - soln1000)) / max(np.abs(euler10 - soln10))
#Q_RK2 = max(np.abs(RK21000 - soln1000)) / max(np.abs(RK210 - soln10))
#print 'Euler: Q =',Q_euler
#print 'RK2: Q =',Q_RK2



# Plotting

pl.clf()

euler, = pl.plot(x,yy[:,0],"r.")
RK2, = pl.plot(x,yy_RK2[:,0],"c.")
solution, = pl.plot(x,soln(x),"k-")

pl.legend((euler, RK2, solution), 
          ('Euler Integrator', 'RK2 Integrator', 'Actual Solution'),
          loc=(0.43, 0.68), frameon=False)

pl.title('Euler and RK2 integrators with 1000 points')
pl.xlabel('x')
pl.ylabel('y(x)')

pl.plot()
#pl.savefig('1000points.pdf')

