import numpy as np
import matplotlib.pyplot as plt

file = "rocket_data.csv"

data = np.loadtxt(file, delimiter=",")

location = data[:, 0] # m
mass = data[:, 1] # kg

# Define tolerances
location_tolerance = 0.002 # Location uncertainty, +/- 2mm
mass_tolerance = 0.1 # Mass uncertainity, +/- 0.1kg

def x_cp(x_data, location_tolerance):  # Ref "The Theoretical Prediction of the Center of Pressure", James Barrowman
    # Assume small angles of attack (< 10 deg)

    x_data = [1.2, 6] # Only the nosecone, body, and fins contribute to CP
    N = np.size(x_data)

    x_with_tolerance = x_data + np.random.uniform(low=-location_tolerance, high=location_tolerance, size=N)
    rocket_length = np.sum(x_with_tolerance)
    # rocket_length = x_with_tolerance[-1]

    # Nose
    CNalpha_nose = 2 # Coefficient of Normal Force W.R.T. AoA for an ogive nosecone
    nose_L = 1.2 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    xbar_nose = rocket_length - 0.466*nose_L # CP location for an ogive nosecone, as measured from the engine

    # Body
    CNalpha_body = 0 # Cylindrical body, axisymmetric
    xbar_body = 0
    body_L = 6 + np.random.uniform(low=-location_tolerance, high=location_tolerance)

    # Check number of fins on Halcyon

    # Fins
    a = 0.8 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    b = 0.7 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    m = (a-b)/2
    l = s = 0.2 + np.random.uniform(low=-location_tolerance, high=location_tolerance)
    n = 4 # Number of fins
    R = 0.2 + np.random.uniform(low=-location_tolerance, high=location_tolerance) # Radius of rocket body
    d = 2*R
    noseTip_to_finRoot = body_L + nose_L # Distance from nose tip to front edge of fin root
    xf = rocket_length - noseTip_to_finRoot # As measured from the engine

    CNalpha_fins_noInterference = (4*n*(s/d)**2) / (1 + np.sqrt(1 + (2*l/(a+b))**2))
    fin_interference = 1 + R/(s+R)
    CNalpha_fins = CNalpha_fins_noInterference * fin_interference
    delta_xf = m*(a+2*b) / 3*(a+b) + 1/6 * (a+b - (a*b/(a+b)))
    xbar_fins = xf + delta_xf # CP location for fins

    # Full Rocket
    CNalpha = CNalpha_nose + CNalpha_body + CNalpha_fins
    xbar = ((CNalpha_nose*xbar_nose) + (CNalpha_body*xbar_body) + (CNalpha_fins*xbar_fins)) / CNalpha

    return xbar

def main():
        N = 1000
        x_cp_list = np.array([])

        for i in range(N):
            x_cp_list = np.append(x_cp_list, x_cp(location, location_tolerance))

        fig, ax = plt.subplots()
        ax.hist(x_cp_list, bins=50)
        plt.show()

if __name__ == "__main__":
    main()