"""
@Author: Louis Sharma
@Date: 2023-08-19
@Description: a python class for generating a simulation of an 2D ideal gas 
comprised of n particles that have the same mass. 
"""

import numpy as np
import matplotlib.pyplot as plt

class IdealGas:

    def __init__(self, nparticles, mass, radius,L, duration, nsteps):
        self.nparticles = nparticles #number of particles 
        self.mass = mass #mass of particles (kg)
        self.radius = radius #radius of particles (m)
        self.L = L #length of the box (m)
        self.duration = duration #duration (s)
        self.nsteps = nsteps #number of steps
        self.t = np.linspace(0,duration,nsteps) #time array 
        self.dt = duration/nsteps #timestep
        #initialize positions and velocities
        self.r = np.stack((np.random.uniform(2*radius,L-2*radius, nparticles), 
                           np.random.uniform(2*radius,L-2*radius, nparticles))).T
        
        v0 = 3 #m/s
        theta = np.random.uniform(0,1,nparticles)*2*np.pi #random angles (uniformly distributed)
        self.v = v0*np.stack((np.cos(theta),np.sin(theta))).T #initial velocities

        
        
    def dij(self,r):

        dist_matrix = np.zeros((self.nparticles, self.nparticles))
        for i in range(self.nparticles):
            for j in range(i, self.nparticles):
                dist_matrix[i,j] = np.linalg.norm(r[i]-r[j])

        return dist_matrix

    

    def check_collisions(self):

        """Checks if there has been a collision between 2 particles. If so, 
        updates the velocities taking into account conservation of energy
        and momentum."""
        
        r_next = self.r + self.v*self.dt #next theoretical position of the particles

        #particle collisions
        dist_pairs = self.dij(r_next)
        rows,cols = np.triu_indices(self.nparticles, k=1)
        collisions = np.where(dist_pairs[rows,cols]<=2*self.radius)[0]
        i,j = rows[collisions], cols[collisions]
        vdiff = self.v[i] - self.v[j]
        rdiff = self.r[i] - self.r[j]
        self.v[i] = self.v[i] -rdiff*(np.sum(rdiff*vdiff,axis=1)/(np.sum(rdiff**2,axis=1)))[:,np.newaxis]
        self.v[j] = self.v[j] +rdiff*(np.sum(rdiff*vdiff,axis=1)/(np.sum(rdiff**2,axis=1)))[:,np.newaxis]


        #wall collisions, invert orthogonal velocity component
        self.v[r_next[:,0]<self.radius, 0] *=-1
        self.v[r_next[:,0]>self.L-self.radius,0] *=-1
        self.v[r_next[:,1]<self.radius,1] *=-1
        self.v[r_next[:,1]>self.L-self.radius,1] *=-1


    def step(self): 
        self.check_collisions() #check for collisions
        self.r += self.v*self.dt #step the positions


    def MaxwellBoltzmann(self, v0, v):

        KE_avg = 1/2*self.mass*np.sum(v0**2) * 1/self.nparticles #average kinetic energy
        kT =KE_avg #temperature
        sigma_sq = kT/self.mass #variance of the distribution
        f = np.exp(-v**2/(2*sigma_sq)) * v/sigma_sq #Maxwell-Boltzmann distribution
        return f

    


    def energy(self, velocities):
        #takes as input the velocity of the particles at a certain time. 
        average_energies = [1/2*self.mass*np.sum(v**2) for v in velocities] #average kinetic energy (=total E/N)
        return average_energies


    def animate(self):

        positions = np.zeros((self.nparticles,2,len(self.t))) #empty array to store positions
        velocities = np.zeros((self.nparticles, len(self.t))) #empty array to store velocity norm

        for n,t in enumerate(self.t): #iterate through all timesteps

            self.step() #step 
            positions[:,:,n] = self.r #append to positions
            velocities[:,n] = np.linalg.norm(self.v, axis=1) #append to velocities

        return positions,velocities




    




