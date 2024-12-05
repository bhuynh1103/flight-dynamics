import numpy as np

def integration_sim(self, dt, num_iterations):
    for i in range(num_iterations):
        if i == 0:
            state = np.zeros(9)
            state[4] = np.deg2rad(90)
            state[3] = 0.001

        dx = self.state_dot(state=state, dt=0.01, gust_intensity=4.5) * dt
        state = dx*dt + state
        return state