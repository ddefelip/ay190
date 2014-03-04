import sys,math
import numpy as np
import matplotlib.pyplot as mpl

# functions
def apply_bcs(y):
    # apply boundary conditions
    # sets the boundary grid points of y (0th and last) to be the 
    # last interior points of y (1st and 2nd to last)
    y[0] = y[1]
    y[-1] = y[-2]

def analytic(x, LL):
    # sin function
    return 0.125*np.sin(2*np.pi*x/LL)



L = 100.
# set up the grid so dx = 0.1
x = np.linspace(0, L, num=1001)
dx = x[1]-x[0]
n = len(x)

##### 
# Change cfl value to observe differing values of v*dt/dx
cfl = 1.
#####

# use max value of function to determine a constant dt 
dt = cfl*(dx/max(np.abs(analytic(x, L))))



#set up initial conditions
y = analytic(x, L)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data - upwind
mpl.plot(x,analytic(x, L),'r-') # analytic initial condition
mpl.show()

yold = np.copy(y)
ntmax = 1000
t = np.linspace(0., ntmax*dt, num=ntmax+1)

for it in range(ntmax):
    # save previous data
    yold = np.copy(y)

    # get new data for upwind scheme
    # note: it matters what direction the velocity is in at a given
    # point because that determines what direction "upwind" is
    for j in range(1, n-1):
        if yold[j] > 0:
            # positive velocity: get data from the left
            y[j] = yold[j] - (yold[j]*dt/dx)*(yold[j] - yold[j-1])
        else:
            # negative velocity: get data from the right
            y[j] = yold[j] - (yold[j]*dt/dx)*(yold[j+1] - yold[j])

            
    # after updates, apply boundary conditions
    apply_bcs(y)
    
    print 'time =',t[it]
    
    # plot data and make a movie!
    mpl.clf()
    data, = mpl.plot(x,y,'x-')
    initial, = mpl.plot(x,analytic(x,L),'r-')
    mpl.title('Upwind Scheme')
    mpl.xlabel('x')
    mpl.ylabel('Psi')
    mpl.legend((data,initial),
               ('Result at t = %g' % t[it], 'Initial Condition (t = 0)'),
                frameon=False,
                loc=(0.52,0.60))
    if t[it] % 40 == 0 or t[it] % 140 == 0:
        mpl.savefig('t=%g.pdf' % t[it])
    mpl.draw()


mpl.show()
mpl.clf()

