import numpy as np
import matplotlib.pyplot as plt

def integration_sim(self, dt, num_iterations):
    for i in range(num_iterations):
        if i == 0:
            state = np.zeros(9)
            state[4] = np.deg2rad(90)
            state[3] = 0.001
            state[1] = 0.01
            state_matrix = np.zeros((num_iterations, 9))

        dx = self.state_dot(state=state, dt=0.01, gust_intensity=4.5)
        state = dx*dt + state
        state_matrix[num_iterations-1, :] = state
        # print(f"x [m] = {state_matrix[:, 0]}, z [m] = {state_matrix[:, 1]}, v_x [m/s] = {state_matrix[:, 2]}, v_z [m/s] = {state_matrix[:, 3]}, theta [rad] = {state_matrix[:, 4]}, q [rad/s] = {state_matrix[:, 5]}")
        # print()

    # print(state_matrix)

    start_time = 0
    end_time = num_iterations * dt
    time = np.arange(start_time, end_time, dt)
    plt.plot(time, state_matrix[:, 1])
    plt.xlabel("Time [s]")
    plt.ylabel("z [m]")
    plt.show()