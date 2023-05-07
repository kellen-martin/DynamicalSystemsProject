import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')
from tools import *

def main():
    # Sun-Earth-Moon IC
    r_0 = np.array([[0.0, 0.0, 0.0], [147.1*10**9, 0, 0], [147.1*10**9 + 384.4*10**6, 0, 0]])
    v_0 = np.array([[0.0, 0.0, 0.0], [0, 29.783*10**3, 0], [0, 29.783*10**3 + 1022, 0]])
    masses = np.array([1.98847*10**30, 5.972*10**24, 7.34767309*10**22])

    # Random ICs
    #r_0 = np.array([[0.0, 0.0, 0.0], [147.1*10**9, 0, 0], [302*10**8, 0, 0]])
    #v_0 = np.array([[7*10**3, -6*10**3, 0.0], [0, 9.783*10**3, 0], [0, 50*10**3 + 1022, 0]])
    #masses = np.array([1.989*10**30, 5.972*10**24, 7.342*10**22])
    t_0 = 0
    t_f = 86400*365
    h = 100

    #r_1,v_1,r_2,v_2,r_3,v_3 = Naive3BP(r_0, v_0, masses, t_0, t_f, h)
    #r_1,v_1,r_2,v_2,r_3,v_3 = Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h)
    r_1,v_1,r_2,v_2,r_3,v_3 = PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h)


    #print(r_3)
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection='3d')
    plt.gca().patch.set_facecolor('black')
    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 1.0)), ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 1.0)), ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
    # plot the trajectory of objects
    ax.plot(r_1[:,0],r_1[:,1],r_1[:,2],'y-')
    ax.plot(r_2[:,0],r_2[:,1],r_2[:,2], 'b-')
    ax.plot(r_3[:,0],r_3[:,1],r_3[:,2], 'r--')

    # plot starting points
    ax.plot(r_1[0,0],r_1[0,1],r_1[0,2],'y*')
    ax.plot(r_2[0,0],r_2[0,1],r_2[0,2], 'b*')
    ax.plot(r_3[0,0],r_3[0,1],r_3[0,2], 'r*')


    # Mark Starting Point

    plt.show()
    return

main()