import numpy as np
import matplotlib.pyplot as plt

path = "mass_placement.csv"

data = np.loadtxt(path, delimiter=",")

x_data    = data[:, 0] # [m]
mass_data = data[:, 1] # [kg]

# Tolerances are +/- unless otherwise defined
x_tolerance    = 2/1000 # [m]
mass_tolerance = 0.1    # [kg]

# debugging
# x_tolerance    = 0
# mass_tolerance = 0

def x_cg(masses, m_tolerance, x_distances, x_tolerance):
    N = np.size(masses)

    m = masses + np.random.uniform(low=-m_tolerance, high=m_tolerance, size=N)
    x = x_distances + np.random.uniform(low=-x_tolerance, high=x_tolerance, size=N)
    x_cg = m @ x / np.sum(m)

    return x_cg

def dryden_turbulence(t, dt, V, h, gust, gust_intensity):
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
    noise = np.normal(0, 1, size=[3, 1])

    gust_next = (1 - V_fps * dt / L) * gust + sigma * np.sqrt(2 * V_fps * dt / L) * noise # [m/s]

    gust_rotated = np.array([gust_next[2], gust_next[1], gust_next[0]]);
    return gust_rotated

def main():
    N = 1000
    x_cg_list = np.array([])

    for i in range(N):
        x_cg_list = np.append(x_cg_list, x_cg(mass_data, mass_tolerance, x_data, x_tolerance))

    fig, ax = plt.subplots()
    ax.hist(x_cg_list, bins=50)
    plt.show()

    t = np.linspace(0, 60, 1000)
    

if __name__ == "__main__":
    main()