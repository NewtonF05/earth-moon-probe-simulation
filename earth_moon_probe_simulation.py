#!/usr/bin/env python
# coding: utf-8

"""
Earth-Moon-Probe Orbital Dynamics Simulation
============================================

Author: Newton Fernihough

A numerical simulation of two coupled orbital-mechanics problems, solved
with SciPy's adaptive Runge-Kutta integrator (``solve_ivp``).

Two interactive problems are available from a menu:

    1. The two-body Earth-Moon problem.
       The Moon is integrated for roughly one lunar month under Earth's
       gravity alone, and its trajectory is plotted in the Earth-centred
       frame.

    2. The three-body Earth-Moon-Probe problem.
       A spacecraft is launched into low lunar orbit (200 km above the
       Moon's surface) while the Moon continues its orbit around the
       Earth. The probe is subject to gravity from both Earth and Moon,
       and its trajectory is plotted alongside the Moon's motion.

Both problems are integrated with tight relative and absolute tolerances
(1e-8 and 1e-11 respectively) to preserve orbital energy over the
integration window.

Run the file and select a problem from the prompt; enter ``q`` to quit.
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define physical constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
m_earth = 5.972e24  # Mass of Earth (kg)
m_moon = 7.348e22   # Mass of Moon (kg)
d_earth_moon = 3.844e8  # Average Earth-Moon distance (m)
v_moon = 1.022e3  # Moon's orbital speed (m/s)
d_earth = 1.2742e7  # Diameter of Earth (m)
d_moon = 3.4748e6  # Diameter of Moon (m)

choice = '0'
while choice != 'q':
    choice = input('Enter a choice: "1" for lunar orbit, "2" for Earth-Moon-Probe, or "q" to quit: ')
    print('You entered:', choice)
    
    if choice == '1':
        print('Selected: Lunar orbit simulation')

        # Function to compute derivatives for Moon's orbit
        def moon_orbit_derivatives(t, state, G, m_earth):
            xm, ym, vxm, vym = state
            rm = np.sqrt(xm**2 + ym**2)
            
            # Acceleration due to Earth's gravity
            axm = -G * m_earth * xm / rm**3
            aym = -G * m_earth * ym / rm**3
            
            return [vxm, vym, axm, aym]

        # Initial conditions: [x position, y position, x velocity, y velocity]
        initial_state_moon = [d_earth_moon, 0, 0, v_moon]

        # Time parameters for the simulation
        t_start = 0
        t_end = 2628288  # ~1 month in seconds
        t_points = np.linspace(t_start, t_end, 1001)

        # Solving the ODE
        moon_solution = solve_ivp(moon_orbit_derivatives, (t_start, t_end), initial_state_moon, t_eval=t_points, args=(G, m_earth), rtol=1e-8,atol=1e-11)

        # Extracting the results for plotting
        x_traj_moon = moon_solution.y[0, :]
        y_traj_moon = moon_solution.y[1, :]

        # Plotting Moon's orbit
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
        ax.set_aspect('equal')
        ax.plot(x_traj_moon, y_traj_moon, label='Moon')
        ax.plot(0, 0, 'ro', label='Earth')  # Earth at the origin
        ax.set_xlabel('x coordinate (m)')
        ax.set_ylabel('y coordinate (m)')
        ax.set_title('Moon Orbit around Earth')
        ax.legend()
        plt.show()

    elif choice == '2':
        print('Selected: Earth-Moon-Probe system')

        # Function to compute derivatives for the Earth-Moon-Probe system
        def probe_moon_derivatives(t, state, G, m_earth, m_moon):
            xm, ym, vxm, vym, xp, yp, vxp, vyp = state
            
            # Distances
            rm = np.sqrt(xm**2 + ym**2)
            rp = np.sqrt(xp**2 + yp**2)
            rpm = np.sqrt((xp - xm)**2 + (yp - ym)**2)

            # Acceleration components for the Moon
            axm = -G * m_earth * xm / rm**3
            aym = -G * m_earth * ym / rm**3

            # Acceleration components for the Probe
            axp = (-G * m_earth * xp / rp**3) - (G * m_moon * (xp - xm) / rpm**3)
            ayp = (-G * m_earth * yp / rp**3) - (G * m_moon * (yp - ym) / rpm**3)
            
            return [vxm, vym, axm, aym, vxp, vyp, axp, ayp]

        # Initial conditions for Moon and Probe
        initial_state_moon = [d_earth_moon, 0, 0, v_moon]

        probe_altitude = 200e3 + d_moon / 2  # 200 km above Moon's surface
        probe_speed = 1.825e3  # Speed for Moon orbit
        initial_state_probe = [d_earth_moon + probe_altitude, 0, 0, v_moon + probe_speed]

        # Combine initial states
        initial_state = initial_state_moon + initial_state_probe

        # Time parameters for the simulation
        t_start = 0
        t_end = 100000  # Simulation duration in seconds
        t_points = np.linspace(t_start, t_end, 1001)

        # Solving the ODE for the combined system
        system_solution = solve_ivp(probe_moon_derivatives, (t_start, t_end), initial_state, t_eval=t_points, args=(G, m_earth, m_moon), rtol=1e-8,atol=1e-11)

        # Extract Moon and Probe trajectories for plotting
        x_traj_moon = system_solution.y[1, :]
        y_traj_moon = system_solution.y[0, :]
        x_traj_probe = system_solution.y[5, :]
        y_traj_probe = system_solution.y[4, :]

        # Plotting the trajectories
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(x_traj_moon, y_traj_moon, label='Moon', color='blue')
        ax.plot(x_traj_probe, y_traj_probe, '--', label='Probe', color='red')
        ax.set_xlabel("x coordinate (m)")
        ax.set_ylabel("y coordinate (m)")
        ax.set_title("Probe and Moon Trajectories")
        ax.set_aspect('equal')
        ax.legend()
        plt.show()

    elif choice != 'q':
        print('Invalid choice, please try again')

print('Exiting the program. Goodbye!')
