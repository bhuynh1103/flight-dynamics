import numpy as np

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

    def barrowman_eqns(self, x_data): # Outputs (Cp, CNalpha, S_ref, aero_moment)
        # Ref "The Theoretical Prediction of the Center of Pressure", James Barrowman
        # Assume small angles of attack (< 10 deg)

        location_tolerance = self.location_tolerance
        
        # Geometry w/ uncertainty
        # TODO: make the locations of the components a function of the object parameters
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
            L_u = h / (0.177 + 0.000823*h)**1.2 # [ft]
            L_v = L_u # [ft]
            L_w = h_ft # [ft]

            sigma_w = gust_intensity # [m/s]
            sigma_u = sigma_w / (0.177 + 0.000823*h_ft)**0.4 # [m/s]
            sigma_v = sigma_u # [m/s]
        elif h_ft >= 2000:
            L_u = 1750 # [ft]
            L_v = 1750 # [ft]
            L_w = 1750 # [ft]

            sigma_u = gust_intensity
            sigma_v = gust_intensity
            sigma_w = gust_intensity
        else:
            # Linearly Interprolate Turbulence Scale Length
            L_u_1000 = 1000 / (0.177 + 0.000823*1000)**1.2
            L_v_1000 = L_u_1000
            L_w_1000 = 1000

            L_u_2000 = 1750
            L_v_2000 = 1750
            L_w_2000 = 1750

            L_u = np.interp(h_ft, [1000, 2000], [L_u_1000, L_u_2000]) # [ft]
            L_v = np.interp(h_ft, [1000, 2000], [L_v_1000, L_v_2000]) # [ft]
            L_w = np.interp(h_ft, [1000, 2000], [L_w_1000, L_w_2000]) # [ft]

            # Linearly Interprolate Turbulence Intensity
            sigma_w_1000 = gust_intensity # [m/s]
            sigma_u_1000 = sigma_w_1000 / (0.177 + 0.000823*1000)**0.4 # [m/s]
            sigma_v_1000 = sigma_u_1000

            sigma_u_2000 = gust_intensity # [m/s]
            sigma_v_2000 = gust_intensity # [m/s]
            sigma_w_2000 = gust_intensity # [m/s]

            sigma_u = np.interp(h_ft, [1000, 2000], [sigma_u_1000, sigma_u_2000]) # [m/s]
            sigma_v = np.interp(h_ft, [1000, 2000], [sigma_v_1000, sigma_v_2000]) # [m/s]
            sigma_w = np.interp(h_ft, [1000, 2000], [sigma_w_1000, sigma_w_2000]) # [m/s]

        L = np.array([L_u, L_v, L_w]) # [ft]
        sigma = np.array([sigma_u, sigma_v, sigma_w]) # [m/s]
        noise = np.random.normal(0, 1, size=[1, 3])

        gust_next = (1 - V_fps * dt / L) * gust + sigma * np.sqrt(2 * V_fps * dt / L) * noise # [m/s]

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
        Cd = self.drag_coefficient
        sm = self.static_margin

        rho_sl = 1.225 # [kg/m^3] - Air Density at Sea Level
        beta = 1 / 9042 # scale length
        rho = rho_sl * np.exp(-beta*z) # [kg/m^3]


        u = vx * np.cos(theta) + vz * np.sin(theta) - gust_u # velocity in x-body axis
        w = vx * np.sin(theta) - vz * np.cos(theta) - gust_w # velocity in z-body axis

        V = np.sqrt(u**2 + w**2) # Freestream velocity

        psi = np.arctan2(vz, vx) # flight path angle
        alpha = np.arctan2(w, u) # angle of attack

        q_bar = 0.5 * rho * V**2 # dynamic pressure
        L = q_bar * S * CNalpha * alpha # lift
        D = q_bar * S * Cd

        # Derivatives
        x_dot     = vx
        z_dot     = vz
        vx_dot    = ( T*np.cos(theta) - D*np.cos(psi) - L*np.sin(psi) ) / M
        vz_dot    = ( T*np.sin(theta) - D*np.sin(psi) + L*np.cos(psi) ) / M
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
            gust_state_next[0],
            gust_state_next[1],
            gust_state_next[2]
        ])

        return state_dot
    
    def integration_sim(self, dt, num_iterations):
        # from state: vx, vz, q
        # constants: mass, thrust
        # calculate: lift, drag, static margin
        # misc.: aero_moment

        for i in range(num_iterations):
            if i == 0:
                state = np.zeros(9)
                state[4] = np.deg2rad(90)
                state[3] = 0.001
                state[1] = 0.01

            dx = self.state_dot(state=state, dt=0.01, gust_intensity=4.5) * dt
            state = dx*dt + state
            return state

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
    for i in range(3):
        halcyon = Rocket(masses=mass_data, mass_locations=x_data, 
                         mass_tolerance=mass_tolerance, mass_location_tolerance=x_tolerance, 
                         aerodynamic_dimensions=dimensions, thrust=TXE2_thrust)
        print(f"Halcyon x_cg = {halcyon.x_cg}")
        print(f"Halcyon x_cp = {halcyon.x_cp}")
        print(f"Halcyon inertia = {halcyon.inertia}")
        print()

        print(halcyon.integration_sim(dt=0.01, num_iterations=100))
        # if i == 1:
        #     state = np.zeros(1, 9)
        # dx = halcyon.state_dot(state=state, dt=0.01, gust_intensity=4.5)
        # state = state + dx*dt
        

if __name__ == "__main__":
    main()