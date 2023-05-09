import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')
from tools import *

def main():
    G = 6.6743*10**(-11)# [m^3/kg*s^2]
    # Sun-Earth-Moon IC
    r_0 = np.array([[0.0, 0.0, 0.0], [147.1*10**9, 0, 0], [147.1*10**9 + 384.4*10**6, 0, 0]])
    v_0 = np.array([[0.0, 0.0, 0.0], [0, 29.783*10**3, 0], [0, 29.783*10**3 + 1022, 0]])
    masses = np.array([1.98847*10**30, 5.972*10**24, 7.34767309*10**22])

    plot_3BP(r_0[0,:])
    # Random ICs
    #r_0 = np.array([[0.0, 0.0, 0.0], [147.1*10**9, 0, 0], [302*10**8, 0, 0]])
    #v_0 = np.array([[7*10**3, -6*10**3, 0.0], [0, 9.783*10**3, 0], [0, 50*10**3 + 1022, 0]])
    #masses = np.array([1.989*10**30, 5.972*10**24, 7.342*10**22])
    t_0 = 0
    t_f = 86400*365*2
    h = 86400

    #r_1,v_1,r_2,v_2,r_3,v_3 = Naive3BP(r_0, v_0, masses, t_0, t_f, h)
    r_1,v_1,r_2,v_2,r_3,v_3 = Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h)
    #r_1,v_1,r_2,v_2,r_3,v_3 = PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h)

    E = total_energy(r_1, v_1, r_2, v_2, r_3, v_3, masses, G)
    #print(r_3)
    r1_ani, r2_ani, r3_ani = prep_animation(r_1, r_2, r_3, 5)
    animate_3BP(r1_ani,r2_ani,r3_ani, "ting.gif")
    #plot_3BP(r_1, r_2, r_3)

    #rela_error = rel_error(E)
    #plot_error(rela_error, t_0, t_f)
    return

main()