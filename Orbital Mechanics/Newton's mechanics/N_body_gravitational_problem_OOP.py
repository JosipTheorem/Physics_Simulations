import numpy as np
#import autograd.numpy as np
#from autograd import grad
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import time
from numba import njit


class Body:
    
    instances = []  ## All bodies are here
    body_count = 0
    positions = []  ## Not using this
    Masses = []

    def __init__(self, mass, position, velocity):
        
        self.mass = mass
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        Body.positions.append(np.array(position)) 
        Body.Masses.append(mass) 
        Body.instances.append(self)  ## Saving body in list of all bodies
        Body.body_count +=1  ## Not using this

        
        #self.Potential_grad = grad(self.Potential) ## For numerical calculation of force from potential, not using

    @classmethod
    def preprocess(cls):
        """Convert lists to NumPy arrays for efficient computation."""
        cls.Masses = np.array(cls.Masses)
        cls.positions = np.array(cls.positions)

    
    def display(self):
        return (self.mass, self.position, self.velocity)
        

    @property
    def Kinetic(self):
        return np.dot(self.velocity, self.velocity) * 0.5 * self.mass
        

    def Potential(self, position): ##Only used for numerics!
        potential_energy = 0.0
        
        for other_body in Body.instances:
            if other_body is not self:  # Exclude the self-interaction term
                r = np.linalg.norm(position - other_body.position)  # Distance between the bodies
                potential_energy -= G * other_body.mass * self.mass / r
        
        return potential_energy
        

    @property 
    def Force_numerical(self):  ## For numerics
        return -self.Potential_grad(self.position)
        
    
    def acceleration(self):
        force = np.array([0.0, 0.0, 0.0])
        for other in Body.instances:
            if other is not self:
                r_vector = self.position - other.position
                r = np.linalg.norm(r_vector)
                force += -G * self.mass * other.mass * r_vector / r**3
        return force/self.mass


    def acceleration_mybe1(self, positions):
        """ Calculate acceleration due to gravitational force from all other bodies """
        # Vectorized calculation of forces from all bodies
        mask = ~(np.all(positions == self.position, axis=1))
        
        diff = positions[mask] - self.position
        
        diff_mass = Body.Masses[mask, np.newaxis] * diff  # Scale each vector by the corresponding mass

        distances = np.linalg.norm(diff, axis=1)  # Calculate distances between bodies
        
        acc = G  * (diff_mass.T / (distances**3 + eps))  # Gravitational force
        
        return np.sum(acc, axis=1)   # Total force divided by mass gives acceleration

    def acceleration_mybe2(self, positions):
        """ Calculate acceleration due to gravitational force from all other bodies """
        # Vectorized calculation of forces from all bodies
        diff = positions - self.position
        
        diff_mass = Body.Masses[:, np.newaxis] * diff  # Scale each vector by the corresponding mass
        
        distances = np.linalg.norm(diff, axis=1)  # Calculate distances between bodies
        
        acc = G  * (diff_mass.T / (distances**3 + eps))  # Gravitational force
        
        acc[:, Body.instances.index(self)] = 0  # Set self-interaction forces to 0

        
        return np.sum(acc, axis=1)   # Total force divided by mass gives acceleration
    


@njit
def fun2(t,y):
    num_bodies = len(masses)
    
    # Reshape y into a 2D array for easier manipulation
    y_reshaped = y.reshape((num_bodies, 6))  # Each row is [x, vx, y, vy, z, vz]
        
    # Extract positions and velocities
    positions = y_reshaped[:, r_es]  # Columns: x, y, z
    velocities = y_reshaped[:, v_es]  # Columns: vx, vy, vz
    
    # Pairwise differences in positions
    diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    
    # Compute pairwise distances manually (Euclidean distance)
    distances = np.sqrt(np.sum(diff**2, axis=2)) + eps  # Adding epsilon to avoid division by zero
        
    # Mask self-interaction
    np.fill_diagonal(distances, np.inf)  # Avoid division by zero for self-interactions
        
    # Compute forces using the gravitational formula
    acceleration = -np.sum(G * masses[np.newaxis, :, np.newaxis] * diff / (distances[:, :, np.newaxis]**3), axis=1)
    
    # Create dydt array (velocity and acceleration components)
    dydt = np.zeros_like(y)
    dydt[0::6] = velocities[:, 0]  # dx/dt = velocity x-component
    dydt[1::6] = acceleration[:, 0]  # dv/dt = acceleration x-component
    dydt[2::6] = velocities[:, 1]  # dy/dt = velocity y-component
    dydt[3::6] = acceleration[:, 1]  # dv/dt = acceleration y-component
    dydt[4::6] = velocities[:, 2]  # dz/dt = velocity z-component
    dydt[5::6] = acceleration[:, 2]  # dv/dt = acceleration z-component

    return dydt


def calculate_energy(Pos,Vel):

    energies = np.empty(steps)
    pos = np.empty([num_bodies*3,3])

    for i in range(steps):
        
        kinetic_i = 0
        # Extract current positions and velocities from the solution at time i
        idx = 0
        j = 0
        for body in Body.instances:
            
            vel = np.array([Vel[idx][i], Vel[idx+1][i], Vel[idx+2][i]])
            kinetic_i = kinetic_i + body.mass*np.dot(vel,vel)/2

            pos[j] = np.array([Pos[idx][i], Pos[idx+1][i], Pos[idx+2][i]])
            
            idx += 3
            j += 1

        # Calculate total kinetic energy
        

        # Calculate total potential energy
        potential_i = 0
        for j in range(num_bodies):
            for k in range(num_bodies):
                if j < k:  # Ensure each pair is counted only once
                    r_vector = pos[j] - pos[k]
                    r = np.linalg.norm(r_vector)
                    potential_i += -G * masses[j] * masses[k] / r

        
        # Total energy at time i
        total_energy = kinetic_i + potential_i
        energies[i] = total_energy

    return energies



if __name__ == '__main__':
    print('yes')

    G=1.0
    eps=1.e-10
    r_e = 1
    m_e = 1

    Earth = Body(mass=m_e, position=[1.0, 0, 0], velocity=[0, 500, 0])
    Mercury = Body(mass=0.055 * m_e, position=[0.39, 0, 0], velocity=[0, 875, 0])
    Venus = Body(mass=0.815 * m_e, position=[0.72, 0, 0], velocity=[0, 622, 0])
    Mars = Body(mass=0.107 * m_e, position=[1.52, 0, 0], velocity=[0, 407, 0])
    Jupiter = Body(mass=318 * m_e, position=[5.2, 0, 0], velocity=[0, 223, 0])
    Saturn = Body(mass=95.2 * m_e, position=[9.58, 0, 0], velocity=[0, 165, 0])
    Uranus = Body(mass=14.5 * m_e, position=[19.2, 0, 0], velocity=[0, 115, 0])
    Neptune = Body(mass=17.1 * m_e, position=[30.1, 0, 0], velocity=[0, 88, 0])
    Sun = Body(mass = m_e * 300000, position = [0,0,0], velocity = [0,0,0])

    Body.preprocess()
    num_bodies = len(Body.instances)
    masses = np.array(Body.Masses)
    r_es = np.array([0, 2, 4])
    v_es = np.array([1, 3, 5])


    initial_conditions = []
    for body in Body.instances:
        initial_conditions.extend([body.position[0], body.velocity[0],
                                body.position[1], body.velocity[1],
                                body.position[2], body.velocity[2]])

    initial_conditions = np.array(initial_conditions)



    start = time.perf_counter()
    tmax = 1
    steps = 10000
    t = np.linspace(0, tmax, steps)
    # Perform numerical integration
    solution = solve_ivp(fun2, [0, tmax], initial_conditions, method='DOP853', t_eval=t, rtol=1e-6, atol=1e-9)
    #solution = solve_ivp(func, [0, tmax], initial_conditions, method='DOP853', t_eval=t)
    end = time.perf_counter()
    print('Calculation time:',end - start,'s')


    Pos = np.empty([num_bodies * 3, steps])
    Vel = np.empty([num_bodies * 3, steps])

    for i in range(0,num_bodies * 6,2):

        Pos[int(i/2)] = solution.y[i]
        Vel[int(i/2)] = solution.y[i+1]


    energies = calculate_energy(Pos,Vel)

    print('Energy E_end-E_initial:', energies[len(energies)-1] - energies[0])

    plt.plot(t,energies)
    plt.ylim(energies[0]+100, energies[0]-100)
    plt.xlabel('t')
    plt.ylabel('E')
    plt.show()


    n_bodies = len(Body.instances)
    body_solutions = {}
    for i in range(n_bodies):
        body_solutions[f'body_{i+1}'] = {
            'x': solution.y[i*6],
            'vx': solution.y[i*6+1],
            'y': solution.y[i*6+2],
            'vy': solution.y[i*6+3],
            'z': solution.y[i*6+4],
            'vz': solution.y[:, i*6+5]
        }

    # Create the figure and 3D axis
    fig10 = plt.figure(10, figsize=(12, 12))
    ax = fig10.add_subplot(111, projection='3d')

    # Loop through each body in the body_solutions dictionary and plot its trajectory
    for i, body in enumerate(body_solutions.values()):
        x = body['x']
        y = body['y']
        z = body['z']
        
        # Use a different color for each body (can use a colormap or just a predefined set of colors)
        color = plt.cm.get_cmap('tab10')(i % 10)  # Use 'tab10' colormap to cycle colors
        
        # Plot the trajectory in 3D
        ax.plot3D(x, y, z, label=f'Body {i+1}', color=color)

    # Set labels for the axes
    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)
    ax.set_zlabel(r'$z$', fontsize=16)
    ax.set_zlim(-1,1)
    #ax.set_xlim(-1000,-600)
    #ax.set_ylim(500,600)


    # Add a legend
    plt.legend(fontsize=16)

    # Show the plot
    plt.show()