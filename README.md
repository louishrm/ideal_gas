# Ideal gas simulation in Python

The purpose of this project is to write a Python code that animate sa 2D ideal gas of $N$ particles evolving in a square box. 
We check that the speeds of the particles at equilibrium follow the **Maxwell-Boltzmann distribution**. 
Check out [this video] on my [channel](https://www.youtube.com/channel/UCrD42YxJdVEKl-9-iER80hA) which explains the code in detail. 


## The Maxwell-Boltzmann distribution

In the 2D ideal gas, the velocities of the particles follow the **Maxwell-Boltzmann distribution**: 

$$f(v) = \frac{m v}{k_B T} e^{-\frac{mv^2}{2k_BT}},$$

where $m$ is the mass of the partiles and the average energy per particle is $\langle E \rangle=k_B T.$ 

## Assumptions
* The particles are represented by circles that all have the same radius and mass. 
* The particles start at random positions and with random velocities. 
* **Important:** when we initialize our particle positions, we need to make sure that the particles aren't overlapping, otheriwse this can lead to numerical errors. The best way to do this is to arrange them on a grid
* The number of particles in the container remains fixed throughout the simulation. 
* The particles can collide *elastically* with each other and the boundaries of the container


## Algorithm 
The idea behind the algorithm to run the simulation is simple: 

* Fix the number of particles, the size of the container and length of simulation.
* Initialize their positions and velocities at time $t=0$. 
* At each time step: check if some of the particles are colliding with one another/the wall. If so, update their velocities accounting for the conservation of energy + momentum.
* Once the velocities have been updated, we get the positions of the particles at the next step: $\vec r(t+\Delta t) = \vec r(t) + \vec v(t) \Delta t$.

## Collisions 
### Particle-particle collisions
We use the following formula for updating the velocities of the particles that are colliding:

$$ \vec v_i' = \vec v_i - \frac{(\vec v_i - \vec v_j) \cdot (\vec r_i - \vec r_j)}{(\vec r_i - \vec r_j) \cdot (\vec r_i - \vec r_j)} (\vec r_i - \vec r_j) $$
$$ \vec v_j' = \vec v_j - \frac{(\vec v_j - \vec v_i) \cdot (\vec r_j - \vec r_i)}{(\vec r_j - \vec r_i) \cdot (\vec r_i - \vec r_j)} (\vec r_j - \vec r_i) $$ 

Source: [wikipedia](https://en.wikipedia.org/wiki/Elastic_collision)

### Particle-wall collisions
When a particle collides with one of the walls, we invert the component of its velocity that is orthogonal to the wall. 

