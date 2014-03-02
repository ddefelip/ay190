#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

# system parameters - use commented values to compute 
# trajectories for stars file on cgs scale. Leave them
# to just look at arcsecond scale (also comment out 
# "seconds_per_year" for t1).
initial_data_file = "sun_earth.asc"#"sgrAstar.asc"
distance_unit_to_cm = 1.#*0.04*8000*1.5e13
time_unit_to_s = 1.#*seconds_per_year
mass_unit_to_g = 1.#*msun
Nsteps = 200
t0 = 0.
t1 = 10. * seconds_per_year
dt = (t1-t0)/Nsteps

# either 2 or 14 for the files for this worksheet
num_objects = 2

final_data_file = "final_stars_positions.asc"



def NbodyRHS(u,mass,time):
    
    rhs = np.zeros((num_objects, 6))
    # rhs has same dimensions as u
    # but vx, vy, vz replace x, y, z, and vx', vy', vz' replace vx, vy, vz

    for i in range(num_objects):
        
        # get altered mass array (without the mass of the iteration
        # we are on)
        a = False
        massi = np.zeros(num_objects-1)
        for m in range(num_objects):
            if m == i:
                a = True 
            elif not a:
                massi[m] = mass[m]   
            else:
                massi[m-1] = mass[m]
                
        # altered distance array (without the coordinates of the
        # iteration we are on so we don't divide by 0)    
        disti = np.zeros((num_objects-1, 3))
        
        # velocities
        rhs[i,:3] = u[i,3:]   
        
        # fill in distance array
        for k in range(3):
            disti[:,k] = np.concatenate(( u[i,k]-u[:i,k], u[i,k]-u[i+1:,k] ))
      
        # compute RHS velocity derivatives
        for k in range(3):
            
            if num_objects == 2:
                rhs[i,k+3] = -ggrav*np.sum(massi*disti[0,k]/np.sum(disti**2)**1.5)
            else:
                for j in range(num_objects-1):
                    rhs[i,k+3] += -ggrav*massi[j]*disti[j,k]/np.sum(disti[j]**2)**1.5
   
    return rhs



def NbodyRK4(u,mass,time,dt):
    
    k1 = dt*NbodyRHS(u,mass,time)
    k2 = dt*NbodyRHS(u+k1/2.,mass,time+dt/2.)
    k3 = dt*NbodyRHS(u+k2/2.,mass,time+dt/2.)
    k4 = dt*NbodyRHS(u+k3,mass,time+dt)
    
    unew = u + (k1 + 2.*k2 + 2.*k3 + k4)/6.
    
    return unew



def TotalEnergy(u,mass,time):
    
    Ekin = 0.
    Epot = 0.
    
    for i in range(num_objects):
        
        Ekin += 0.5*mass[i]*np.sum(u[i,3:]**2)
        
        # don't want to double count pairs of objects, so
        # the j array gets smaller after each iteration
        for j in range(i+1, num_objects):
            Epot += -ggrav*mass[i]*mass[j] / np.sum((u[i,:3]-u[j,:3])**2)**0.5

    return Ekin + Epot
    
    
# main program
plt.ion()

(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)


# convert from units in initial data file to cgs
# comment out if you want to keep scales from stars data file
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s 
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
# change '2.5' to get larger or smaller scales (I went as high as 10)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()
energy = np.zeros(Nsteps)

for it in np.arange(0, Nsteps):
    
    time = t0 + it * dt
    u = NbodyRK4(u,mass,time,dt)
    
    if it % max(1,Nsteps/100) == 0:
        
        print "it = %d, time = %g years, energy = %g" % \
              (it, time / seconds_per_year,
               TotalEnergy(u,mass,time))
               
        plt.clf()
        fig = plt.gcf()
        ax = mpl3d.Axes3D(fig)
        ax.scatter(u[:,0],u[:,1],u[:,2])
        ax.set_xlim((-rmax,rmax))
        ax.set_ylim((-rmax,rmax))
        ax.set_zlim((-rmax,rmax))
        plt.draw()
      
    energy[it] = TotalEnergy(u,mass,time)

# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u, header=file_header)

# Plot energy
plt.clf()
plt.plot(np.arange(0,10,(t1-t0)/(200*seconds_per_year)),energy,'r')
plt.title('Total Energy vs Time')
plt.xlabel('Time (years)')
plt.ylabel('Total Energy (erg)')

