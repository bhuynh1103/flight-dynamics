import numpy as np
import matplotlib.pyplot as plt

file = "rocket_data.csv"

data = np.loadtxt(file, delimiter=",")

location = data[:, 0] # m
mass = data[:, 1] # kg

# Define tolerances
x_tolerance = 0.002 # Location uncertainty, +/- 2mm
mass_tolerance = 0.1 # Mass uncertainity, +/- 0.1kg

def x_cp(x_data, location_tolerance):
    # Ref "The Theoretical Prediction of the Center of Pressure", James Barrowman
    # Assume small angles of attack (< 10 deg)

    # Geometry w/ uncertainty
    nose_L = 1.2 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    body_L = 6 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    body_D = 0.4 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    nozzle_L = 0.8 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    fin_cr = 0.8 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    fin_ct = 0.7 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    fin_station = 1.2 + np.random.uniform(low=-location_tolerance, high=location_tolerance) # Distance from the fins to end of body
    fin_D = (0.8 + np.random.uniform(low=-location_tolerance, high=location_tolerance)) * np.sqrt(2) # Fin-to-fin diameter
    fin_span = (fin_D-body_D) / 2
    rocket_length = nose_L + body_L + nozzle_L

    # N = np.size(x_data)
    # x_data = [1.2, 6, 0.8] # Nosecone, body, nozzle
    # x_with_tolerance = x_data + np.random.uniform(low=-location_tolerance, high=location_tolerance, size=N)
    # rocket_length = np.sum(x_with_tolerance)
    # rocket_length = x_with_tolerance[-1]

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
    # delta_xf = m*(a+2*b) / 3*(a+b) + 1/6 * (a+b - (a*b/(a+b)))
    xbar_fins = xf + delta_xf # CP location for fins

    # Full Rocket
    CNalpha = CNalpha_nose + CNalpha_body + CNalpha_fins
    xbar = ((CNalpha_nose*xbar_nose) + (CNalpha_body*xbar_body) + (CNalpha_fins*xbar_fins)) / CNalpha

    # More Parameters
    S_ref = np.pi * body_R**2
    cp_moment = (xbar_fins**2 * CNalpha_fins + + xbar_body**2 * CNalpha_body + xbar_nose**2 * CNalpha_nose) / CNalpha

    return (xbar, CNalpha, S_ref, cp_moment)

def main():
    print(x_cp(x_data=location, location_tolerance=x_tolerance))
        # N = 1000
        # x_cp_list = np.array([])

        # for i in range(N):
        #     get_xcp = x_cp(x_data=location, location_tolerance=x_tolerance)
        #     x_cp_list = np.append(x_cp_list, get_xcp[0])

        # fig, ax = plt.subplots()
        # ax.hist(x_cp_list, bins=50)
        # plt.show()

if __name__ == "__main__":
    main()