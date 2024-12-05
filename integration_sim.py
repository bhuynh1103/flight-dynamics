def integration_sim(state, dx, dt, num_iterations)

    # from state: vx, vz, q
    # constants: mass, thrust
    # calculate: lift, drag, static margin
    # misc.: aero_moment

    for i in range(num_iterations):
         state = dx*dt + state