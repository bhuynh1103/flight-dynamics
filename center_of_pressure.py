import numpy as np
import matplotlib.pyplot as plt

file = "rocket_data.csv"
data = np.loadtxt(file, delimiter=",") # Geometric measurements of rocket features, m

# Define tolerances
x_tolerance = 0.002 # Location uncertainty, m

def x_cp(x_data, location_tolerance):
 
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


def main():
    N = 1000
    x_cp_list = CNalpha_list = S_ref_list = aero_moment_list = np.array([])
    for i in range(N):
        [xbar, CNalpha, S_ref, aero_moment] = x_cp(x_data=data, location_tolerance=x_tolerance)
        x_cp_list = np.append(x_cp_list, xbar)
        CNalpha_list = np.append(CNalpha_list, CNalpha)
        S_ref_list = np.append(S_ref_list, S_ref)
        aero_moment_list = np.append(aero_moment_list, aero_moment)

    fig, ax = plt.subplots()
    ax.hist(x_cp_list, bins=50)
    plt.show()

if __name__ == "__main__":
    main()