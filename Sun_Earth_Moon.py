import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')
from tools import *

def main():
    # Set Initial Conditions
    r_0 = np.array([[0.0, 0.0, 0.0], [147.1*10**9, 0, 0], [147.1*10**9 + 384.4*10**6, 0, 0]])   # [m]
    v_0 = np.array([[0.0, 0.0, 0.0], [0, 29.783*10**3, 0], [0, 29.783*10**3 + 1022, 0]])        # [m/s]
    masses = np.array([1.98847*10**30, 5.972*10**24, 7.34767309*10**22])                        # [kg]
    G = 6.6743*10**(-11)                                                                        # [m^3/kg*s^2]

    # Time and Step Size
    t_0 = 0             # [s]
    t_f = 86400*365.24*2  # [s] Seconds/day * days/year * years
    h = 86400           # [s]

    # Integrators
    #r_1,v_1,r_2,v_2,r_3,v_3 = Naive3BP(r_0, v_0, masses, t_0, t_f, h)
    #r_1,v_1,r_2,v_2,r_3,v_3 = Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h)
    r_1,v_1,r_2,v_2,r_3,v_3 = PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h)

    h2 = 10000
    #r_12,v_12,r_22,v_22,r_32,v_32 = Naive3BP(r_0, v_0, masses, t_0, t_f, h2)
    #r_12,v_12,r_22,v_22,r_32,v_32 = Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h2)
    r_12,v_12,r_22,v_22,r_32,v_32 = PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h2)

    h3 = 100
    #r_13,v_13,r_23,v_23,r_33,v_33 = Naive3BP(r_0, v_0, masses, t_0, t_f, h3)
    #r_13,v_13,r_23,v_23,r_33,v_33 = Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h3)
    r_13,v_13,r_23,v_23,r_33,v_33 = PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h3)

    # Calculate Energy
    #E = total_energy(r_1, v_1, r_2, v_2, r_3, v_3, masses, G)
    #E2 = total_energy(r_12, v_12, r_22, v_22, r_32, v_32, masses, G)
    #E3 = total_energy(r_13, v_13, r_23, v_23, r_33, v_33, masses, G)

    # Errors
    #re = rel_error(E)
    #mre = max_rel_error(E)
    #re2 = rel_error(E2)
    #mre2 = max_rel_error(E2)
    #re3 = rel_error(E3)
    #mre3 = max_rel_error(E3)

    #ae = absolute_error(E)
    #mae = max_abs_error(E)
    #ae2 = absolute_error(E2)
    #mae2 = max_abs_error(E2)
    #ae3 = absolute_error(E3)
    #mae3 = max_abs_error(E3)

    #plot_abserrors(ae, ae2, ae3, t_0, t_f)
    #plot_relerrors(re, re2, re3, t_0, t_f)
    #plot_max_abserror(mae, mae2, mae3, h, h2, h3)
    #plot_max_relerror(mre, mre2, mre3, h, h2, h3)
    
    # Plots
    #plot_3BP(r_1, r_2, r_3)
    #plot_3BP(r_12, r_22, r_32)
    #plot_3BP(r_13, r_23, r_33)

    # Animation
    r1_ani,r2_ani,r3_ani = prep_animation(r_1, r_2, r_3, 10)
    r12_ani,r22_ani,r32_ani = prep_animation(r_12, r_22, r_32, 10)
    r13_ani,r23_ani,r33_ani = prep_animation(r_13, r_23, r_33, 1000)
    animate_3BP(r1_ani,r2_ani,r3_ani, "P_day.gif")
    animate_3BP(r12_ani,r22_ani,r32_ani, "P_10000.gif")
    animate_3BP(r13_ani,r23_ani,r33_ani, "P_100.gif")

    return

main()