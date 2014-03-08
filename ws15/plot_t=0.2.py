import numpy as np
import matplotlib.pyplot as pl

pl.clf()

(itSPH,rSPH,rhoSPH,pressSPH,vSPH) = \
                 np.loadtxt('./out_0160.dat',unpack=True)
                 
(itexact,rexact,rhoexact,pressexact,vexact) = \
                 np.loadtxt('./output',unpack=True,skiprows=2)
                 

# Density plot
densitySPH, = pl.plot(rSPH,rhoSPH,'b+')
densityexact, = pl.plot(rexact,rhoexact,'r+')
pl.title('Density Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Density')
pl.legend((densitySPH,densityexact),
          ('SPH','Exact'),frameon=False,loc=(0.6,0.7))
pl.savefig('density.pdf')

pl.clf()

# Velocity plot
velocitySPH, = pl.plot(rSPH,vSPH,'b+')
velocityexact, = pl.plot(rexact,vexact,'r+')
pl.title('Velocity Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Velocity')
pl.legend((velocitySPH,velocityexact),
          ('SPH','Exact'),frameon=False,loc=(0.45,0.3))
pl.savefig('velocity.pdf')

pl.clf()


# Pressure plot       
pressureSPH, = pl.plot(rSPH,pressSPH,'b+')
pressureexact, = pl.plot(rexact,pressexact,'r+')
pl.title('Pressure Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Pressure')
pl.legend((pressureSPH,pressureexact),
          ('SPH','Exact'),frameon=False,loc=(0.6,0.7))
pl.savefig('pressure.pdf')


