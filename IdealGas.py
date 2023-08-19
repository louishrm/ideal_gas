"""
@Author: Louis Sharma
@Date: 2023-08-19
@Description: a python class for generating a simulation of an 2D ideal gas 
comprised of n particles that have the same mass. 
"""

import numpy as np
import matplotlib.pyplot as plt

class IdealGas:

    def __init__(self, nparticles, mass, radius, x1,x2, y1,y2, t):   
        self.nparticles = nparticles #number of particles 
        self.mass = mass #mass of particles kg
        self.radius = radius #radius of particles m 
        self.x1 = x1 #boundaries of the container, all in m
        self.x2 = x2 
        self.y1 = y1 
        self.y2 = y2
        self.t = t #array of times for the animation in seconds
        self.fps = 60 #frames per second
        self.dt = t[1]-t[0] #timestep in seconds

        #initialize positions and velocities
        self.r = np.stack((np.random.uniform(x1+2*radius,x2-2*radius, nparticles), 
                           np.random.uniform(y1+2*radius,y2-2*radius, nparticles))).reshape((nparticles,2))
        
        self.v = np.stack((np.random.uniform(-15,15, nparticles), 
                           np.random.uniform(-15,15, nparticles))).reshape((nparticles,2))


        
        
    def dij(self,r):
        """Computes the distance between all pairs of particles. They are in the upper triangle 
        of a distance matrix."""

        dist_matrix = np.zeros((self.nparticles, self.nparticles)) #initialize distance matrix

        for i in range(self.nparticles): #loop through all distinct pairs i and j
            for j in range(i, self.nparticles):

                dist_matrix[i,j] = np.linalg.norm(r[i]-r[j]) #append the norm of the difference vector

        return dist_matrix

    

    def check_collisions(self):

        """Checks if there has been a collision between 2 particles. If so, 
        updates the velocities taking into account conservation of energy
        and momentum."""
        
        r_next = self.r + self.v*self.dt #next theoretical position of the particles

        #particle collisions
        dist_pairs = self.dij(r_next) #get the distance matrix at the next timestep
        rows,cols = np.triu_indices(self.nparticles, k=1) #get the indices of the upper triangle
        collisions = np.where(dist_pairs[rows,cols]<=2*self.radius)[0] #find those where the particles are colliding
        i,j = rows[collisions], cols[collisions]

        #update velocities
        vdiff = self.v[i] - self.v[j]
        rdiff = self.r[i] - self.r[j]
        self.v[i] = self.v[i] -rdiff*(np.sum(rdiff*vdiff,axis=1)/(np.sum(rdiff**2,axis=1)))[:,np.newaxis]
        self.v[j] = self.v[j] +rdiff*(np.sum(rdiff*vdiff,axis=1)/(np.sum(rdiff**2,axis=1)))[:,np.newaxis]


        #wall collisions, invert orthogonal velocity component
        self.v[r_next[:,0]<self.x1+self.radius, 0] *=-1
        self.v[r_next[:,0]>self.x2-self.radius,0] *=-1
        self.v[r_next[:,1]<self.y1+self.radius,1] *=-1
        self.v[r_next[:,1]>self.y2-self.radius,1] *=-1


    def step(self): 
        """Steps the positions of the particles forward in time."""

        self.check_collisions() #check for collisions
        self.r += self.v*self.dt #step the positions by integrating the velocity


    def MaxwellBoltzmann(self, v0, v):
        """2D Maxwell-Boltzmann distribution for an ideal gas."""

        KE_avg = 1/2*self.mass*np.sum(v0**2) * 1/self.nparticles #average kinetic energy
        kT =KE_avg #temperature
        sigma_sq = kT/self.mass #variance of the distribution
        f = np.exp(-v**2/(2*sigma_sq)) * v/sigma_sq #Maxwell-Boltzmann distribution
        return f

    


    def energy(self, velocities):
        """Computes the energy of the system at a certain time. 
        Used to check that the energy is correctly conserved."""

        average_energies = [1/self.nparticles*0.5*self.mass*np.sum(v**2) for v in velocities] #average kinetic energy (=total E/N)
        return average_energies


    def animate(self):
        """Animates the simulation. Returns 2 arrays, positions and velocities recording the 
        positions and velocities of the particles at each timestep."""

        positions = np.zeros((self.nparticles,2,len(self.t))) #empty array to store positions
        velocities = np.zeros((self.nparticles, len(self.t))) #empty array to store velocity norm

        for n,t in enumerate(self.t): #iterate through all timesteps
            self.step() #step to get next state
            positions[:,:,n] = self.r #append to positions
            velocities[:,n] = np.linalg.norm(self.v, axis=1) #append to velocities

        return positions,velocities




    


