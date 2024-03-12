# N-Body Simulation

# Overview
This Python script simulates the dynamics of N planets interacting gravitationally over time. It employs the Leap-Frog integration scheme (KDK) to update the positions and velocities of the planets. The script visualizes the evolution of the system's energy and the trajectories of the planets in real-time or after the simulation has completed.

# Features
- Simulation of N-body dynamics using gravitational forces.
- Visualization of planetary trajectories and energy conservation throughout the simulation.
- Option to animate the simulation in real-time.

# Requirements
- Python 3.x environment.
- Packages: numpy, matplotlib, tqdm.

# Usage
1) Ensure Python and the required packages are installed on your system.
2) Open a terminal or command prompt.
3) Navigate to the directory containing the Python script (n_body_simulation.py).
4) Run the script using the command python n_body_simulation.py.
5) The script will execute the simulation and display visualizations of planetary trajectories and energy evolution.

# Configuration
- The script allows configuration of parameters such as the number of planets, simulation duration, time step, and initial conditions of the planets.
- Users can adjust the softening length to handle close encounters between planets.

# Visualization
The script generates two plots:
  1) Planetary trajectories in a 2D space, with an option to animate in real-time.
  2) Energy evolution over the course of the simulation.

# Notes
- The simulation uses the Leap-Frog integration scheme (KDK) for numerical stability.
- The script may take longer to execute for larger numbers of planets or longer simulation durations.
- Gas planets are so far away from the Sun that terrestrial planets become as small as points. Zoom in is requires.
- Users can customize the script to incorporate additional features or analysis as needed.
