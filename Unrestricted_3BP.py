import numpy as np
import math
import matplotlib.pyplot as plt
from random import random
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')
from tools import *

def main():
    # Set Initial Conditions
    r_0 = np.array([[-15*10**9, 7*10**8, -100*10**9], [147.1*10**9, -15*10**8, 105*10**8], [-147.1*10**9 + 384.4*10**6, 205*10**7, 84*10**6]])   # [m]
    v_0 = np.array([[-2*10, 10**3, 45], [-2.5*10**3, 16*10, 52], [-1.7-10**3, 1.2*10, 4*10]])        # [m/s]
    masses = np.array([5.972*10**28, 5.972*10**28, 5.972*10**28])                               # [kg]
    G = 6.6743*10**(-11)                                                                        # [m^3/kg*s^2]
    print(r_0[0,:])
    
    r_02 = np.zeros_like(r_0)
    v_02 = np.zeros_like(v_0)
    masses_2 = np.zeros_like(masses)
    r_03 = np.zeros_like(r_0)
    v_03 = np.zeros_like(v_0)
    masses_3 = np.zeros_like(masses)
    for i in range(0,3):
        masses_2[i] = masses[i]
        masses_3[i] = masses[i]
        for j in range(0,3):
            r_02[i,j] = r_0[i,j] + random()*(r_0[i,j]/(10**6))
            v_02[i,j] = v_0[i,j] + random()*(v_0[i,j]/(10**6))
            r_03[i,j] = r_0[i,j] + random()*(r_0[i,j]/(10**6))
            v_03[i,j] = v_0[i,j] + random()*(v_0[i,j]/(10**6))

    print(r_0 - r_02)
    print(v_0 - v_02)
    # Time and Step Size
    t_0 = 0             # [s]
    t_f = 0  # [s] Seconds/day * days/year * years
    h = 5000           # [s]

    # Integrators

    r_1,v_1,r_2,v_2,r_3,v_3 = PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h)
    r_12,v_12,r_22,v_22,r_32,v_32 = PEFRL_3BP(r_02, v_02, masses_2, t_0, t_f, h)
    r_13,v_13,r_23,v_23,r_33,v_33 = PEFRL_3BP(r_03, v_03, masses_3, t_0, t_f, h)


    # Calculate Energy
    #E = total_energy(r_1, v_1, r_2, v_2, r_3, v_3, masses, G)
    #E2 = total_energy(r_12, v_12, r_22, v_22, r_32, v_32, masses_2, G)
    #E3 = total_energy(r_13, v_13, r_23, v_23, r_33, v_33, masses_3, G)
    #plot_energy(E)
    #plot_energies(E,E2,E3,t_0,t_f)

    # Errors
    #re = rel_error(E)
    #ae = absolute_error(E)

    r1_ani, r2_ani, r3_ani = prep_animation(r_1, r_2, r_3, 500)
    r12_ani, r22_ani, r32_ani = prep_animation(r_12, r_22, r_32, 500)
    r13_ani, r23_ani, r33_ani = prep_animation(r_13, r_23, r_33, 500)

    animate_3BP(r1_ani, r2_ani, r3_ani,"unpert.gif")
    animate_3BP(r12_ani, r22_ani, r32_ani,"pert1.gif")
    animate_3BP(r13_ani, r23_ani, r33_ani,"pert2.gif")
    animate_3BP(r1_ani, r12_ani, r13_ani, "body1.gif")
    animate_3BP(r2_ani, r22_ani, r23_ani, "body2.gif")
    animate_3BP(r3_ani, r32_ani, r33_ani, "body3.gif")
    # Plots
    plot_3BP(r_1, r_2, r_3)
    plot_3BP(r_12, r_22, r_32)
    plot_3BP(r_13, r_23, r_33)
    plot_3BP(r_1, r_12, r_13)
    plot_3BP(r_2, r_22, r_23)
    plot_3BP(r_3, r_32, r_33)
    # plot_difference(r_1, r_12,t_0,t_f)
    return

main()