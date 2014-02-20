import numpy as np
import matplotlib.pyplot as pl
import scipy as sp

ggrav = 6.67e-8

zone, mass, rad, temp, rho, infall, efrac, ang =\
       np.loadtxt('./presupernova.dat', unpack=True)

# column 1 is mass - g
# column 2 is radius - cm
# column 3 is temperature - K
# column 4 is density - g/cm^3
# column 5 is infall velocity
# column 6 is electron fraction
# column 7 is angular velocity (this star is non-rotating)

pl.clf()

### log-log plot of density vs radius
#pl.plot(np.log10(rad), np.log10(rho), 'r')
#pl.title('Density vs Radius - Log-Log')
#pl.xlabel('Log10(radius) - cm')
#pl.ylabel('Log10(density) - g/cm^3')
#pl.xlim(min(np.log10(rad)), 9)
#pl.ylim(4, max(np.log10(rho))*1.05)
#pl.savefig('density-loglog.pdf')



# edit rad and rho to suit needs
rad = np.concatenate(([0.], rad[rad <= 10.**9]))
rho = np.concatenate((rho[:1], rho[:(rad.size-1)]))

# get total mass within rad = 10^9 cm
tot_mass = mass[rad.size-1]



# use linear interpolation
# p(x) = rho(x_i) + rhoprime*(x-x_i)

newrad = np.linspace(0., max(rad), num=6001)

rhoprime = (rho[1:] - rho[:-1]) / (rad[1:] - rad[:-1])

newrho = np.zeros(newrad.size)

for i in range(newrad.size - 1):
    for j in range(rad.size):
        if newrad[i] < rad[j]:
            
            newrho[i] = rho[j-1] + rhoprime[j-1]*(newrad[i] - rad[j-1])
            break

# last value is simply copied because it is defined to be at the
# same max value of the radius array       
newrho[newrad.size-1] = rho[rad.size-1]


#### Euler Forward Integration Functions ###
def RHS(rad,rho,z):
    
    # RHS function
    rhs = np.zeros(2)
    if rad == 0:
        rhs[0] = 0.
        rhs[1] = 4.*np.pi*ggrav*rho ### trying z/r = 1?
    else:
        rhs[0] = z
        rhs[1] = 4.*np.pi*ggrav*rho - 2.*z/rad

    return rhs

def integrate_FE(rad,dr,rho,phi,z):

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = phi
    old[1] = z

    # forward Euler integrator
    new = old + dr*RHS(rad,rho,z)

    return (new[0], new[1])

#######################################

dr = newrad[1] - newrad[0]
zarray = np.zeros(newrad.size)
phiarray = np.zeros(newrad.size)

# do this initially, and correct phi later with an additive constant
phiarray[0] = 0.
zarray[0] = 4.*np.pi*ggrav*newrho[0]



###### uncomment (and remove [n] from newrho below) 
###### to check that it works 
## set rho to some constant
#newrho = 42.0
## boudary conditions for sphere of constant density
#phi_outer_test = -(4./3.)*np.pi*ggrav*newrho*newrad[-1]**2
#phi_0_test = -2.*np.pi*ggrav*newrho*newrad[-1]**2
#################################

#  Forward Euler  #
for n in range(newrad.size-1):
        
    (phiarray[n+1],zarray[n+1]) = integrate_FE(newrad[n],
                                               dr,
                                               newrho[n],
                                               phiarray[n],
                                               zarray[n])


#### uncomment to check that it works ###### 
# adjust solution so it matches boudary condition there
#phiarray += phi_0_test
## print relative error
#print 'error at boundary'
#print (phi_outer_test - phiarray[-1]) / phi_outer_test
#################################

# correct phi to satisfy outer boundary condition:
phiarray += -ggrav*tot_mass/max(newrad) - phiarray[-1]


# plot potential
pl.plot(newrad,phiarray,'b')
pl.title('Star\'s Gravitaional Potential vs Radius')
pl.xlabel('Radius - cm')
pl.ylabel('Potential - erg')
#pl.savefig('grav_pot.pdf')



