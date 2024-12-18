import numpy as np
import matplotlib.pyplot as plt

class Rocket:
    def __init__(self, masses, mass_locations, mass_tolerance, mass_location_tolerance, aerodynamic_dimensions, thrust):
        self.n = np.size(masses) # Number of point masses provided
        self.thrust = thrust

        # Randomly define masses and locations of masses
        self.masses   = masses + np.random.uniform(low=-mass_tolerance, high=mass_tolerance, size=self.n) # [kg]
        self.x_masses = mass_locations + np.random.uniform(low=-mass_location_tolerance, high=mass_location_tolerance, size=self.n) # [m]
        self.mass_tolerance = mass_tolerance
        self.location_tolerance = mass_location_tolerance
        self.total_mass = np.sum(self.masses) # [kg]

        # Aerodynamic Properties of Rocket
        self.dimensions = aerodynamic_dimensions
        self.drag_coefficient = 0.5

        # Get CP, CNalpha, S_ref, and aero_moment
        barrowman_outputs = self.barrowman_eqns(self.dimensions)
        self.x_cp         = barrowman_outputs[0] # [m]
        self.CNalpha      = barrowman_outputs[1]
        self.S_ref        = barrowman_outputs[2] # [m^2]
        self.aero_moment  = barrowman_outputs[3] # [N*m]

        # Get CG
        self.x_cg = self.get_x_cg() # [m]

        self.static_margin = self.x_cg - self.x_cp

        # Get moment of inertia about CG
        self.inertia = self.get_inertia() # [kg*m^2]

        # Stability Derivatives
        self.CNq = self.CNalpha * self.static_margin
        self.CMq = self.CNalpha * self.aero_moment + 2*self.CNalpha*self.x_cg*self.x_cp - self.CNalpha*self.x_cg**2

    def barrowman_eqns(self, x_data): # Outputs (Cp, CNalpha, S_ref, aero_moment)
        # Ref "The Theoretical Prediction of the Center of Pressure", James Barrowman
        # Assume small angles of attack (< 10 deg)

        location_tolerance = self.location_tolerance
        
        # Geometry w/ uncertainty
        nose_L = x_data[0] + np.random.uniform(low=-location_tolerance, high=location_tolerance)
        body_L = x_data[1] + np.random.uniform(low=-location_tolerance, high=location_tolerance)
        body_D = x_data[2] + np.random.uniform(low=-location_tolerance, high=location_tolerance)
        nozzle_L = x_data[3] + np.random.uniform(low=-location_tolerance, high=location_tolerance)
        fin_cr = x_data[4] + np.random.uniform(low=-location_tolerance, high=location_tolerance)
        fin_ct = x_data[5] + np.random.uniform(low=-location_tolerance, high=location_tolerance)
        fin_D = (x_data[6] + np.random.uniform(low=-location_tolerance, high=location_tolerance)) * np.sqrt(2) # Fin-to-fin diameter
        fin_span = (fin_D-body_D) / 2
        rocket_length = nose_L + body_L + nozzle_L

        # Nose
        CNalpha_nose = 2 # Coefficient of Normal Force W.R.T. AoA for an ogive nosecone
        xbar_nose = rocket_length - 0.466*nose_L # CP location for an ogive nosecone, as measured from the engine

        # Body
        CNalpha_body = 0 # Cylindrical body, axisymmetric
        xbar_body = 0

        # Fins
        x_r = m = (fin_cr-fin_ct)/2
        n = 4 # Number of fins
        body_R = body_D/2 # Body radius
        x_b = body_L + nose_L # Distance from nose tip to front edge of fin root
        xf = rocket_length - x_b # As measured from the engine

        CNalpha_fins_noInterference = (4*n*(fin_span/body_D)**2) / (1 + np.sqrt(1 + (2*fin_span/(fin_cr+fin_ct))**2))
        fin_interference = 1 + body_R/(fin_span+body_R)
        CNalpha_fins = CNalpha_fins_noInterference * fin_interference
        delta_xf = x_r/3 * (fin_cr + 2*fin_ct) / (fin_cr + fin_ct) + 1/6 * ((fin_cr + fin_ct) - (fin_cr*fin_ct)/((fin_cr + fin_ct)))
        xbar_fins = xf + delta_xf # CP location for fins

        # Full Rocket
        CNalpha = CNalpha_nose + CNalpha_body + CNalpha_fins 
        xbar = ((CNalpha_nose*xbar_nose) + (CNalpha_body*xbar_body) + (CNalpha_fins*xbar_fins)) / CNalpha # CP location for rocket

        # More Parameters
        S_ref = np.pi * body_R**2 # Reference surface area
        aero_moment = (xbar_fins**2 * CNalpha_fins + xbar_body**2 * CNalpha_body + xbar_nose**2 * CNalpha_nose) / CNalpha # Aerodynamic moment about the CP

        return (xbar, CNalpha, S_ref, aero_moment)

    def dryden_turbulence(self, dt, V, h, gust, gust_intensity):
        # t, dt given in [s]
        # h given in [m]
        # gust is vector given in [m/s]
        # gust is given in [m/s]

        V_fps = V * 3.28084 # [ft/s]
        h_ft = h * 3.28084 # [ft]

        if h_ft <= 1000:
            L_u = (h_ft / (0.177 + 0.000823*h_ft)**1.2) / 3.28084 # [m]
            L_v = L_u # [m]
            L_w = h # [m]

            sigma_w = gust_intensity # [m/s]
            sigma_u = sigma_w / (0.177 + 0.000823*h_ft)**0.4 # [m/s]
            sigma_v = sigma_u # [m/s]
        elif h_ft >= 2000:
            L_u = 1750 / 3.28084 # [m]
            L_v = 1750 / 3.28084 # [m]
            L_w = 1750 / 3.28084 # [m]

            sigma_u = gust_intensity # [m/s]
            sigma_v = gust_intensity # [m/s]
            sigma_w = gust_intensity # [m/s]
        else:
            # Linearly Interprolate Turbulence Scale Length
            L_u_1000 = (1000 / (0.177 + 0.000823*1000)**1.2) # [ft]
            L_v_1000 = L_u_1000 # [ft]
            L_w_1000 = 1000 # [ft]

            L_u_2000 = 1750 # [ft]
            L_v_2000 = 1750 # [ft]
            L_w_2000 = 1750 # [ft]

            L_u = np.interp(h_ft, [1000, 2000], [L_u_1000, L_u_2000]) / 3.28084 # [m]
            L_v = np.interp(h_ft, [1000, 2000], [L_v_1000, L_v_2000]) / 3.28084 # [m]
            L_w = np.interp(h_ft, [1000, 2000], [L_w_1000, L_w_2000]) / 3.28084 # [m]

            # Linearly Interprolate Turbulence Intensity
            sigma_w_1000 = gust_intensity # [m/s]
            sigma_u_1000 = sigma_w_1000 / (0.177 + 0.000823*1000)**0.4 # [m/s]
            sigma_v_1000 = sigma_u_1000 # [m/s]

            sigma_u_2000 = gust_intensity # [m/s]
            sigma_v_2000 = gust_intensity # [m/s]
            sigma_w_2000 = gust_intensity # [m/s]

            sigma_u = np.interp(h_ft, [1000, 2000], [sigma_u_1000, sigma_u_2000]) # [m/s]
            sigma_v = np.interp(h_ft, [1000, 2000], [sigma_v_1000, sigma_v_2000]) # [m/s]
            sigma_w = np.interp(h_ft, [1000, 2000], [sigma_w_1000, sigma_w_2000]) # [m/s]

        L = np.array([L_u, L_v, L_w]) # [ft]
        sigma = np.array([sigma_u, sigma_v, sigma_w]) # [m/s]
        noise = np.random.normal(0, 1, size=[1, 3])

        gust_next = np.exp(-V_fps * dt / L) * gust + sigma * np.sqrt(1 - np.exp(-2 * V_fps * dt / L)) * noise # [m/s]

        gust_rotated = np.array([gust_next[0][2], gust_next[0][1], gust_next[0][0]])
        return gust_rotated
    
    def get_x_cg(self): # get location of center of gravity
        m = self.masses
        x = self.x_masses

        x_cg = m @ x / np.sum(m)
        return x_cg
    
    def get_inertia(self): # get moment of inertia about center of gravity
        m = self.masses
        x = self.x_masses
        x_cg = self.x_cg
        r = x - x_cg

        I = np.sum(m * r**2)
        return I

    def state_dot(self, state, dt, gust_intensity):
        # Define states
        x      = state[0]
        z      = state[1]
        vx     = state[2]
        vz     = state[3]
        theta  = state[4]
        q      = state[5]
        gust_u = state[6]
        gust_v = state[7]
        gust_w = state[8]

        gust_state = np.array([gust_u, gust_v, gust_w])

        T = self.thrust
        M = self.total_mass
        I = self.inertia
        S = self.S_ref
        CNalpha = self.CNalpha
        CNq = self.CNq
        Cd = self.drag_coefficient
        sm = self.static_margin

        rho_sl = 1.225 # [kg/m^3] - Air Density at Sea Level
        beta = 1 / 9042 # scale length
        rho = rho_sl * np.exp(-beta*z) # [kg/m^3]


        u = vx * np.cos(theta) + vz * np.sin(theta)  - gust_u # velocity in x-body axis
        w = vx * np.sin(theta) - vz * np.cos(theta)  - gust_w # velocity in z-body axis

        V = np.sqrt(u**2 + w**2) # Freestream velocity

        psi = np.arctan2(vz, vx) # flight path angle
        alpha = np.arctan2(w, u) # angle of attack

        q_bar = 0.5 * rho * V**2 # dynamic pressure
        N = q_bar * S * (CNalpha*alpha + CNq*q) # normal force
        L = N # lift
        D = q_bar * S * Cd
        g = 9.81 # [m/s^2]

        # Derivatives
        x_dot     = vx
        z_dot     = vz
        vx_dot    = ( T*np.cos(theta) - D*np.cos(psi) - L*np.sin(psi) ) / M
        vz_dot    = ( T*np.sin(theta) - D*np.sin(psi) + L*np.cos(psi) - M*g) / M
        theta_dot = q
        q_dot     = ( -L*np.cos(alpha)*sm -D*np.sin(alpha)*sm ) / I

        gust_state_next = self.dryden_turbulence(dt, V, z, gust_state, gust_intensity)

        state_dot = np.array([
            x_dot,
            z_dot,
            vx_dot,
            vz_dot,
            theta_dot,
            q_dot,

            alpha,

            gust_state_next[0],
            gust_state_next[1],
            gust_state_next[2]
        ])

        return state_dot
    
    def integration_sim(self, dt, num_iterations):
        alpha_list = np.array([])
        for i in range(num_iterations):
            if i == 0:
                state = np.zeros(9)
                state[4] = np.deg2rad(90)
                state[3] = 0.001
                state[1] = 0.01
                state_matrix = np.zeros((num_iterations, 7)) # State Matrix with the states as columns and times as rows

            dx = self.state_dot(state=state, dt=0.01, gust_intensity=4.5)
            state[0:6] = dx[0:6]*dt + state[0:6] # Integrate the first 6 states
            state[6] = dx[6] # Storage for alpha
            state[7:9] = dx[7:9] # Last 3 states are gusts, not to be integrated
            state_matrix[i, :] = state[0:7] 

        return state_matrix

def main():
    path = "mass_placement.csv"
    data = np.loadtxt(path, delimiter=",")

    x_data    = data[:, 0] # [m]
    mass_data = data[:, 1] # [kg]

    path = "rocket_data.csv"
    dimensions = np.loadtxt(path, delimiter=",")

    # Tolerances are +/- unless otherwise defined
    x_tolerance    = 2/1000 # [m]
    mass_tolerance = 0.1    # [kg]

    TXE2_thrust = 15.5 * 1000 # [N]

    dt = 0.01
    num_iterations = 6000
    start_time = 0
    end_time = num_iterations * dt
    time = np.arange(start_time, end_time, dt)

    N = 500
    x_list = z_list = vx_list = vz_list = theta_list = q_list = alpha_list = np.array([])

    runs = np.empty(shape=(N, num_iterations, 7)) # states as columns, time as rows, every matrix layer is a run

    fig1, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1)


    for i in range(N):
        halcyon = Rocket(masses=mass_data, mass_locations=x_data, 
                         mass_tolerance=mass_tolerance, mass_location_tolerance=x_tolerance, 
                         aerodynamic_dimensions=dimensions, thrust=TXE2_thrust)

        halcyon_state = halcyon.integration_sim(dt=dt, num_iterations=num_iterations)

        x_list = np.append(x_list, halcyon_state[-1][0])
        z_list = np.append(z_list, halcyon_state[-1][1])
        vx_list = np.append(vx_list, halcyon_state[-1][2])
        vz_list = np.append(vz_list, halcyon_state[-1][3])
        theta_list = np.append(theta_list, halcyon_state[-1][4])
        q_list = np.append(q_list, halcyon_state[-1][5])
        alpha_list = np.append(alpha_list, halcyon_state[-1][6])

        runs[i] = halcyon_state

    
    mean = np.mean(runs, axis=0)
    std = np.std(runs, axis=0)

    ax1.plot(time, mean[:, 0] + 2 * std[:, 0], 'b')
    ax1.plot(time, mean[:, 0] - 2 * std[:, 0], 'b')
    ax1.plot(time, mean[:, 0], 'r')

    ax2.plot(time, mean[:, 1] + 2 * std[:, 1], 'b')
    ax2.plot(time, mean[:, 1] - 2 * std[:, 1], 'b')
    ax2.plot(time, mean[:, 1], 'r')

    ax3.plot(time, mean[:, 2] + 2 * std[:, 2], 'b')
    ax3.plot(time, mean[:, 2] - 2 * std[:, 2], 'b')
    ax3.plot(time, mean[:, 2], 'r')

    ax4.plot(time, mean[:, 3] + 2 * std[:, 3], 'b')
    ax4.plot(time, mean[:, 3] - 2 * std[:, 3], 'b')
    ax4.plot(time, mean[:, 3], 'r')

    ax5.plot(time, mean[:, 4] + 2 * std[:, 4], 'b')
    ax5.plot(time, mean[:, 4] - 2 * std[:, 4], 'b')
    ax5.plot(time, mean[:, 4], 'r')

    ax6.plot(time, mean[:, 5] + 2 * std[:, 5], 'b')
    ax6.plot(time, mean[:, 5] - 2 * std[:, 5], 'b')
    ax6.plot(time, mean[:, 5], 'r')

    ax7.plot(time, mean[:, 6] + 2 * std[:, 6], 'b')
    ax7.plot(time, mean[:, 6] - 2 * std[:, 6], 'b')
    ax7.plot(time, mean[:, 6], 'r')

    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('$x$ (m)')

    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('$z$ (m)')

    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('$v_x$ (m/s)')

    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('$v_z$ (m/s)')

    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('$\\theta$ (rad)')

    ax6.set_xlabel('Time (s)')
    ax6.set_ylabel('q (rad/s)')

    ax7.set_xlabel('Time (s)')
    ax7.set_ylabel('$\\alpha$ (rad)')

    fig1.suptitle('Rocket States vs. Time for %d Halcyon Objects' % N)

    fig2, (ax11, ax22, ax33, ax44, ax55, ax66, ax77) = plt.subplots(7, 1)

    ax11.hist(x_list, bins=50)
    ax11.set_xlabel('distance from launch (m)')
    ax11.set_ylabel('count')

    ax22.hist(z_list, bins=50)
    ax22.set_xlabel('altitude (m)')
    ax22.set_ylabel('count')

    ax33.hist(vx_list, bins=50)
    ax33.set_xlabel('horizotnal velocity (m/s)')
    ax33.set_ylabel('count')

    ax44.hist(vz_list, bins=50)
    ax44.set_xlabel('vertical velocity (m/s)')
    ax44.set_ylabel('count')

    ax55.hist(theta_list, bins=50)
    ax55.set_xlabel('pitch angle (rad)')
    ax55.set_ylabel('count')

    ax66.hist(q_list, bins=50)
    ax66.set_xlabel('pitch rate (rad/s)')
    ax66.set_ylabel('count')

    ax77.hist(alpha_list, bins=50)
    ax77.set_xlabel('angle of attack (rad)')
    ax77.set_ylabel('count')

    fig2.suptitle('Distribution of Final Rocket States')
    plt.show()
        

if __name__ == "__main__":
    main()