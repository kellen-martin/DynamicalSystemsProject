import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')
from tools import *

def main():
    # Relative Satelite Initial Conditions
    Earth_SOI = 9.24*10**8
    v_mag_1 = 1000      # [m/s]
    v_mag_2 = 2000
    v_mag_3 = 3000
    eta = 0             # [radians]

    deltav_i_1 = v_mag_1*np.cos(eta)
    deltav_j_1 = v_mag_1*np.sin(eta)
    deltav_i_2 = v_mag_2*np.cos(eta)
    deltav_j_2 = v_mag_2*np.sin(eta)
    deltav_i_3 = v_mag_3*np.cos(eta)
    deltav_j_3 = v_mag_3*np.sin(eta)
  

    # Set Initial Conditions
    r_0 = np.array([[0.0, 0.0, 0.0], [147.1*10**9, 0, 0], [147.1*10**9, Earth_SOI, 0]])                          # [m]
    v_0_1 = np.array([[0.0, 0.0, 0.0], [0, 29.783*10**3, 0], [deltav_i_1, 29.783*10**3 + deltav_j_1, 0]])        # [m/s]
    v_0_2 = np.array([[0.0, 0.0, 0.0], [0, 29.783*10**3, 0], [deltav_i_2, 29.783*10**3 + deltav_j_2, 0]])        # [m/s]
    v_0_3 = np.array([[0.0, 0.0, 0.0], [0, 29.783*10**3, 0], [deltav_i_3, 29.783*10**3 + deltav_j_3, 0]])        # [m/s]
    masses = np.array([1.98847*10**30, 5.972*10**24, 0])                                                         # [kg]
    G = 6.6743*10**(-11)                                                                                         # [m^3/kg*s^2]

    # Time and Step Size
    t_0 = 0                  # [s]
    t_f = 86400*365.24*1     # [s] Seconds/day * days/year * years
    h = 1000                 # [s]

    # Integrators
    r_1,v_1,r_2_1,v_2,r_3_1,v_3 = PEFRL_3BP(r_0, v_0_1, masses, t_0, t_f, h)
    r_1,v_1,r_2_2,v_2,r_3_2,v_3 = PEFRL_3BP(r_0, v_0_2, masses, t_0, t_f, h)
    r_1,v_1,r_2_3,v_2,r_3_3,v_3 = PEFRL_3BP(r_0, v_0_3, masses, t_0, t_f, h)

    # Find Relative Position
    delta_r_1 = r_3_1 - r_2_1
    delta_r_2 = r_3_2 - r_2_2
    delta_r_3 = r_3_3 - r_2_3

    # Circle of Radius Earth SOI
    theta = np.linspace( 0 , 2 * np.pi , 150 )
    radius = Earth_SOI
    a = radius * np.cos( theta )
    b = radius * np.sin( theta )

    # Plot s/c relative position
    plt.figure(figsize=(12,12))
    plt.plot(delta_r_1[:,0],delta_r_1[:,1], label = "1000 m/s")
    plt.plot(delta_r_2[:,0],delta_r_2[:,1], label = "2000 m/s")
    plt.plot(delta_r_3[:,0],delta_r_3[:,1], label = "3000 m/s")
    plt.plot(a,b, label = "Earth SOI")
    plt.plot(c,d, label = ".1 AU")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.legend()
    plt.title("1.5 Year s/c Earth Relative Position [eta = 0]")
    plt.axis('equal')
    plt.show()
    
    # Plots
    plot_3BP(r_1, r_2_1, r_3_1)
    plot_3BP(r_1, r_2_2, r_3_2)
    plot_3BP(r_1, r_2_3, r_3_3)


    # Animation
    r1_ani_1,r2_ani_1,r3_ani_1 = prep_animation(r_1, r_2_1, r_3_1, 200)
    r1_ani_2,r2_ani_2,r3_ani_2 = prep_animation(r_1, r_2_2, r_3_2, 200)
    r1_ani_3,r2_ani_3,r3_ani_3 = prep_animation(r_1, r_2_3, r_3_3, 200)
    animate_3BP(r1_ani_1,r2_ani_1,r3_ani_1, "Earth_Flyby_1.gif")
    animate_3BP(r1_ani_2,r2_ani_2,r3_ani_2, "Earth_Flyby_2.gif")
    animate_3BP(r1_ani_3,r2_ani_3,r3_ani_3, "Earth_Flyby_3.gif")

    return

main()