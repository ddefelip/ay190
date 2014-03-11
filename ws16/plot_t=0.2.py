import numpy as np
import matplotlib.pyplot as pl

pl.clf()

(itSPH,rSPH,rhoSPH,pressSPH,vSPH) = \
                 np.loadtxt('./t=0.2_SPH.dat',unpack=True)
                 
(itexact,rexact,rhoexact,pressexact,vexact) = \
                 np.loadtxt('./t=0.2_exact.dat',unpack=True,skiprows=2)
         
         
(rPC,rhoPC,pressPC,vPC) = \
                 np.loadtxt('./t=0.2_pc.dat',unpack=True)
                 
(rMC,rhoMC,pressMC,vMC) = \
                 np.loadtxt('./t=0.2_mc.dat',unpack=True)
                 
(rMINMOD,rhoMINMOD,pressMINMOD,vMINMOD) = \
                 np.loadtxt('./t=0.2_minmod.dat',unpack=True)


# Density plot comparing reconstruction methods
densityPC, = pl.plot(rPC,rhoPC,'r+')
densityMINMOD, = pl.plot(rMINMOD,rhoMINMOD,'g+')
densityMC, = pl.plot(rMC,rhoMC,'b+')
pl.title('Density Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Density')
pl.legend((densityPC,densityMINMOD,densityMC),
          ('Piecewise Constant','Minmod-limited','Monotonized Central'),
           frameon=False,loc=(0.5,0.7))
pl.ylim(0.20,1.05)
pl.savefig('reconstruction_densities.pdf')
# Best one is clearly Monotonized Central, so use that going forward

pl.clf()

# Compare to results from ws15

# Density plot
densitySPH, = pl.plot(rSPH,rhoSPH,'r+')
densityexact, = pl.plot(rexact,rhoexact,'g+')
densityMC, = pl.plot(rMC,rhoMC,'b+')
pl.title('Density Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Density')
pl.legend((densitySPH,densityexact,densityMC),
          ('SPH','Exact','MC2'),frameon=False,loc=(0.6,0.7))
pl.savefig('density_with_ws15.pdf')


pl.clf()

# Velocity plot
velocitySPH, = pl.plot(rSPH,vSPH,'r+')
velocityexact, = pl.plot(rexact,vexact,'g+')
velocityMC, = pl.plot(rMC,vMC,'b+')
pl.title('Velocity Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Velocity')
pl.legend((velocitySPH,velocityexact,velocityMC),
          ('SPH','Exact','MC2'),frameon=False,loc=(0.45,0.3))
pl.savefig('velocity_with_ws15.pdf')



pl.clf()

# Pressure plot       
pressureSPH, = pl.plot(rSPH,pressSPH,'r+')
pressureexact, = pl.plot(rexact,pressexact,'g+')
pressureMC, = pl.plot(rMC,pressMC,'b+')
pl.title('Pressure Profiles: t=0.2')
pl.xlabel('Position')
pl.ylabel('Pressure')
pl.legend((pressureSPH,pressureexact,pressureMC),
          ('SPH','Exact','MC2'),frameon=False,loc=(0.6,0.7))
pl.savefig('pressure_with_ws15.pdf')


