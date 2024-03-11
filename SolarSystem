# Load packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pylab as p
from IPython import display
from tqdm import tqdm


# Constants
N = 10  # Number of Planets
t = 0 # Initial time
tEnd = 1500 * 24 * 60**2  # End time
dt = 1 * 24 * 60**2  # Time step
Num_ts = int(tEnd / dt)  # Number of time steps
G = 6.67428e-11 # Gravitational constant
AU = 149.6e6 * 1000  # Astronomical unit [m]
SM = 1.99e30  # Solar mass [kg]
CMv = 29.789e3  # Earths orbital velocity [km/s]
softening = 0.0 * AU  # Softening length - Can fix 'close encounters'
plotRealTime = False # To animate or not to animate. Really slows speed

def Acceleration(position, mass, G, softening):
    '''
        Calculate the Acceleration given position, mass, G and softening
        Creates an array to store acceleration values in x, y and z
        Loop through all pair separations to calculate denominator in Newton Equation
        and store them in the right coordinate.
        Returns acceleration
    '''

    N = position.shape[0]
    a = np.zeros((N,3))

    for i in range(N):
        for j in range(N):
            if i != j:
                # Find all separations
                dx = position[j,0] - position[i,0]
                dy = position[j,1] - position[i,1]
                dz = position[j,2] - position[i,2]

                denom = (dx**2 + dy**2 + dz**2 + softening**2)**(-1.5)

                a[i,0] +=  G * (dx * denom) * mass[j] 
                a[i,1] +=  G * (dy * denom) * mass[j]
                a[i,2] +=  G * (dz * denom) * mass[j]

    return a

def Energy(position, velocity, mass, G):
    '''
        Returns Potential and Kinetic energy
    '''
     
    # Kinetic Energy
    kin = 0.5 * np.sum( mass * velocity**2 )

    # Vectorize
    x = position[:,0:1]
    y = position[:,1:2]
    z = position[:,2:3]
    
    dx = x.T - x
    dy = y.T - y
    dz = z.T - z

    # Denominator term
    denom = np.sqrt(dx**2 + dy**2 + dz**2)
    denom[denom>0] = denom[denom>0]
    
    # Potential Energy
    # Only include the diagonal triangle.
    # This avoids double counting. Otherwise Energy will be off
    # Diagonal above main Diagonal elements (,1)
    pot = G * np.sum(np.triu(-(mass*mass.T)/denom,1))
    
    return kin, pot

def initialize():
    '''
    Initialize configuration.
    Returns mass and initial position and velocity of the objects of interest
    '''

    # Initial masses are set to the masses of the planets
    mass = np.ones((N, 1))
    mass[0], mass[1], mass[2] = 1 * SM, 0.33e24, 4.87e24
    mass[3], mass[4], mass[5] = 5.97e24, 0.642e24, 1898e24
    mass[6], mass[7], mass[8] = 568e24, 86.8e24, 102e24
    mass[9] = 0.0130e24

    # Initial positions of the planets
    position = np.ones((N, 3))
    position[0], position[1], position[2] = (0, 0, 0), (0.47 * AU, 0, 0), (0.728 * AU, 0, 0)
    position[3], position[4], position[5] = (1.02 * AU, 0, 0), (1.67 * AU, 0, 0), (5.45 * AU, 0, 0)
    position[6], position[7], position[8] = (10.0 * AU, 0, 0), (20.1 * AU, 0, 0), (30.3 * AU, 0, 0)
    position[9] = (49.9 * AU, 0, 0)

    # Initial velocities of the planets
    velocity = np.ones((N, 3))
    velocity[0], velocity[1], velocity[2] = (0, 0, 0), (0, 38.86e3, 0), (0, 34.79e3, 0)
    velocity[3], velocity[4], velocity[5] = (0, 29.29e3, 0), (0, 21.97e3, 0), (0, 12.44e3, 0)
    velocity[6], velocity[7], velocity[8] = (0, 9.09e3, 0), (0, 6.80e3, 0), (0, 5.37e3, 0)
    velocity[9] = (0, 3.71e3, 0)

    return mass, position, velocity

def main():
    '''
    Main function for running and visualizing the simulation
    '''

    mass, position, velocity = initialize()
    a = Acceleration(position, mass, G, softening)
    kin, pot = Energy(position, velocity, mass, G)

    # Save positions, energies and time for plotting purposes
    position_save = np.zeros((N, 3, Num_ts + 1)) # Save position for every particle in the x, y and z direction for every single time step
    position_save[:, :, 0] = position
    kin_save = np.zeros(Num_ts + 1) # Save energies every time step to track total Energy is conserved at the end
    kin_save[0] = kin
    pot_save = np.zeros(Num_ts + 1)
    pot_save[0] = pot
    time_save = np.arange(Num_ts + 1) * dt

    # Prepare figures
    fig = plt.figure(figsize=(8, 8), dpi=80)
    grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
    ax1 = plt.subplot(grid[0:2, 0])
    ax2 = plt.subplot(grid[2, 0])

    # Plot starting positions
    plt.sca(ax1)
    plt.plot(position[1] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[2] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[3] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[4] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[5] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[6] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[7] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[8] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[9] * AU, 0, 'o', color='black', ms=3)
    plt.plot(position[0], 0, 'o', color='black', ms=3)
    ax1.set(xlim=(-51 * AU, 51 * AU), ylim=(-51 * AU, 51 * AU))
    ax1.set_aspect('equal', 'box')

    for i in tqdm(range(0, Num_ts)):
        # Perform Leap-Frog Integration (KDK)
        velocity += (dt / 2) * a # Kick
        position += dt * velocity # Drift
        a = Acceleration(position, mass, G, softening) # Update acceleation
        velocity += (dt / 2) * a # Kick

        t += dt # Update time
        
        kin, pot = Energy(position, velocity, mass, G) # Check energy
        kin_save[i + 1] = kin
        pot_save[i + 1] = pot

        position_save[:, :, i + 1] = position

        if plotRealTime or (i == Num_ts - 1):
            plt.sca(ax1)
            plt.cla()
            plt.draw()
            xx = position_save[:, 0, :]
            yy = position_save[:, 1, :]
            # Scatter position as trail
            plt.scatter(xx, yy, s=1, color=[.7, .7, 1])
            plt.scatter(position[:, 0], position[:, 1], s=13, color='blue')
    
            display.clear_output(wait=True)
            display.display(plt.gcf())
    
            plt.sca(ax2)
            plt.cla()
            # Energy plot
            plt.scatter(time_save, kin_save, color='red', s=1, label='Kin')
            plt.scatter(time_save, pot_save, color='blue', s=1, label='Pot')
            plt.scatter(time_save, kin_save + pot_save, color='black', s=1, label='$E_{tot}$')
            plt.legend()
            plt.ylabel('Energy [J]')


main()
