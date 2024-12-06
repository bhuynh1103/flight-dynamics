import numpy as np
import matplotlib.pyplot as plt

def integration_sim(self, dt, num_iterations):
        for i in range(num_iterations):
            if i == 0:
                state = np.zeros(9)
                state[4] = np.deg2rad(90)
                state[3] = 0.001
                state[1] = 0.01
                state_matrix = np.zeros((num_iterations, 6))

            dx = self.state_dot(state=state, dt=0.01, gust_intensity=4.5)
            state[0:6] = dx[0:6]*dt + state[0:6]
            state[6:9] = dx[6:9]
            state_matrix[i, :] = state[0:6]
        
        # start_time = 0
        # end_time = num_iterations * dt
        # time = np.arange(start_time, end_time, dt)

        # fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1)

        # ax1.plot(time, state_matrix[:, 0])

        # ax2.plot(time, state_matrix[:, 1])

        # ax3.plot(time, state_matrix[:, 2])

        # ax4.plot(time, state_matrix[:, 3])

        # ax5.plot(time, state_matrix[:, 4])

        # ax6.plot(time, state_matrix[:, 5])

        # plt.show()

        return state_matrix

        # plt.plot(time, state_matrix[:, 3])
        # plt.xlabel("Time [s]")
        # plt.ylabel("v_z [m/s]")

        # print(f"x = {state[0]}")
        # print(f"z = {state[1]}")
        # print(f"vx = {state[2]}")
        # print(f"vz = {state[3]}")
        # print(f"theta = {state[4]}")
        # print(f"q = {state[5]}")
        # print(f"gust_u = {state[6]}")
        # print(f"gust_v = {state[7]}")
        # print(f"gust_w = {state[8]}")
        # print()