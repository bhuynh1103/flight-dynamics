import numpy as np

class Rocket:
    def __init__(self, masses, mass_locations, mass_tolerance, location_tolerance, thrust):
        self.masses = masses
        self.x_masses = mass_locations
        self.mass_tolerance = mass_tolerance
        self.location_tolerance = location_tolerance

        # Find CP, CNalpha, S_ref, and aero_moment
        barrowman_outputs = self.barrowman_eqns()
        self.x_cp = barrowman_outputs[0]
        self.CNalpha = barrowman_outputs[1]
        self.S_ref = barrowman_outputs[2]
        self.aero_moment = barrowman_outputs[3]

        # Find CG
        self.x_cg = self.get_x_cg()

    def barrowman_eqns(self): # Outputs (Cp, CNalpha, S_ref, aero_moment)
        # Ref "The Theoretical Prediction of the Center of Pressure", James Barrowman
        # Assume small angles of attack (< 10 deg)

        # Geometry w/ uncertainty
        # TODO: make the locations of the components a function of the object parameters
        nose_L = 1.2 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)
        body_L = 6 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)
        body_D = 0.4 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)
        nozzle_L = 0.8 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)
        fin_cr = 0.8 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)
        fin_ct = 0.7 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)
        fin_station = 1.2 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance) # Distance from the fins to end of body
        fin_D = (0.8 + np.random.uniform(low=-self.location_tolerance, high=self.location_tolerance)) * np.sqrt(2) # Fin-to-fin diameter
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

    def get_x_cg(self):
        N = np.size(self.masses)
        m_tolerance = self.mass_tolerance
        x_tolerance = self.location_tolerance

        m = self.masses + np.random.uniform(low=-m_tolerance, high=m_tolerance, size=N)
        x = self.x_masses + np.random.uniform(low=-x_tolerance, high=x_tolerance, size=N)
        x_cg = m @ x / np.sum(m)

        return x_cg

def main():
    path = "mass_placement.csv"
    data = np.loadtxt(path, delimiter=",")

    x_data    = data[:, 0] # [m]
    mass_data = data[:, 1] # [kg]

    # Tolerances are +/- unless otherwise defined
    x_tolerance    = 2/1000 # [m]
    mass_tolerance = 0.1    # [kg]

    TXE2_thrust = 15.5 * 1000 # [N]

    for i in range(3):
        halcyon = Rocket(masses=mass_data, mass_locations=x_data, 
                         mass_tolerance=mass_tolerance, location_tolerance=x_tolerance, thrust=TXE2_thrust)
        print(f"Halcyon x_cg = {halcyon.x_cg}")
        print(f"Halcyon x_cp = {halcyon.x_cp}")
        print()

if __name__ == "__main__":
    main()