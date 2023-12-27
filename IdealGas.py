"""
@Author: Louis Sharma
@Date: 2023-08-19
@Description: a python class for generating a simulation of an 2D ideal gas. 
"""

import numpy as np
from itertools import product

class IdealGas:

    def __init__(self, N, mass, radius,L, v0, duration, nsteps):

        self.N = N #number of particles 
        self.mass = mass #mass of particles (kg)
        self.radius = radius #radius of particles (m)
        self.L = L #length of the box (m)
        self.duration = duration #duration (s)
        self.nsteps = nsteps #number of steps
        self.dt = duration/nsteps #timestep (s)
        self.v0 = v0 #initial velocity (m/s)

        #intialise random positions and velocities
        grid_size = int(np.ceil(np.sqrt(N))) # Create enough positions in grid for the particles
        spacing = L/ grid_size 
        x = np.linspace(radius + spacing/2, L - radius - spacing/2, grid_size) 
        pos = list(product(x, x))
        
        self.r = np.array(pos[:N]) # Take the positions needed

        # Retrieve specified number of particle positions
        theta = np.random.uniform(0, 2*np.pi, size=N)
        vx,vy = self.v0*np.cos(theta), self.v0*np.sin(theta)
        self.v = np.stack((vx,vy), axis=1)


    def check_collisions(self):
        """Checks for particle-particle collisions and particle-wall collisions.
        If collisions are found, update the velocities of the particles that are colliding, accounting for 
        the conservation of energy and momentum."""
        r_next = self.r + self.v*self.dt #positions at the next timestep
        #check for collisions with the walls
        self.v[r_next[:,0] < self.radius,0] *= -1 #collisions on the left wall
        self.v[r_next[:,0] > self.L-self.radius,0] *= -1 #collisions on the right wall
        self.v[r_next[:,1] < self.radius,1] *= -1 #collisions on the bottom wall
        self.v[r_next[:,1] > self.L-self.radius,1] *= -1 #collisions on the top wall 


        #check for collisions between particles
        for i in range(self.N):
            for j in range(i+1, self.N):

                if np.linalg.norm(r_next[i] - r_next[j]) < 2*self.radius:

                    rdiff = self.r[i] - self.r[j] #vector between particle 1 and particle 2
                    vdiff = self.v[i] - self.v[j]
                    self.v[i] = self.v[i] - rdiff.dot(vdiff)/rdiff.dot(rdiff)*rdiff #update velocity of particle i
                    self.v[j] = self.v[j] + rdiff.dot(vdiff)/rdiff.dot(rdiff)*rdiff #update velocity of particle j

    def step(self):
        """Computes the positions at the next timestep."""
        self.check_collisions()
        self.r += self.v*self.dt



    def animate(self):
        """Evolves the ideal gas for the specified number of steps."""

        positions = np.zeros((self.nsteps, self.N,2)) #empty array to store positions
        speeds = np.zeros((self.nsteps, self.N)) #empty array to store velocity norms

        for n in range(self.nsteps): #iterate through all timesteps

            positions[n,:,:] = self.r #append to positions
            speeds[n,:] = np.linalg.norm(self.v, axis=1) #append to velocities
            self.step() #step 

        return positions,speeds
    

    def MaxwellBoltzmann(self,v):

        KE_avg = 1/2*self.mass*np.sum(self.v0**2) #average kinetic energy
        kT =KE_avg #temperature
        sigma_sq = kT/self.mass #spread of the distribution
        f = np.exp(-v**2/(2*sigma_sq)) * v/sigma_sq #Maxwell-Boltzmann distribution
        return f




    




