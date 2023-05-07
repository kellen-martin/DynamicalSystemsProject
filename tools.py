# Usefull functions for simulating 2 and 3 body problem
import numpy as np
import math
import matplotlib.pyplot as plt

def Naive3BP(r_0, v_0, masses, t_0, t_f, h):

    # Gravitational Constant
    G = 6.6743*10**(-11)# [m^3/kg*s^2]


    # Find number of steps
    steps = int((t_f + t_0)/h)

    # Initialize Vector
    r_1 = np.zeros((steps,3))
    r_2 = np.zeros((steps,3))
    r_3 = np.zeros((steps,3))
    v_1 = np.zeros((steps,3))
    v_2 = np.zeros((steps,3))
    v_3 = np.zeros((steps,3))

    # Enter initial conditions
    r_1[0,0:3] = r_0[0,:]
    r_2[0,0:3] = r_0[1,:]
    r_3[0,0:3] = r_0[2,:]
    v_1[0,0:3] = v_0[0,:]
    v_2[0,0:3] = v_0[1,:]
    v_3[0,0:3] = v_0[2,:]

    m_1 = masses[0]
    m_2 = masses[1]
    m_3 = masses[2]


    for i in range(1,steps):
        # Relative Positions
        r_12 = r_1[i-1,:] - r_2[i-1,:]
        r_13 = r_1[i-1,:] - r_3[i-1,:]
        r_23 = r_2[i-1,:] - r_3[i-1,:]
        r_21 = -r_12
        r_31 = -r_13
        r_32 = -r_23

        # Calculate Accelerations
        a_1 = -G*m_2*(r_12)/(distance(r_12)**3) - G*m_3*(r_13)/(distance(r_13)**3)
        a_2 = -G*m_1*(r_21)/(distance(r_21)**3) - G*m_3*(r_23)/(distance(r_23)**3)
        a_3 = -G*m_1*(r_31)/(distance(r_31)**3) - G*m_2*(r_32)/(distance(r_32)**3)
        if i == 1:
            print(a_3)

        # Update Vectors
        v_1[i,:] = v_1[i-1,:] + h*a_1
        v_2[i,:] = v_2[i-1,:] + h*a_2
        v_3[i,:] = v_3[i-1,:] + h*a_3

        r_1[i,:] = r_1[i-1,:] + h*v_1[i-1,:]
        r_2[i,:] = r_2[i-1,:] + h*v_2[i-1,:]
        r_3[i,:] = r_3[i-1,:] + h*v_3[i-1,:]

    return r_1, v_1, r_2, v_2, r_3, v_3

def Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h):
    # Gravitational Constant
    G = 6.6743*10**(-11)# [m^3/kg*s^2]

    # Find number of steps
    steps = int((t_f + t_0)/h)

    # Initialize Vector
    r_1 = np.zeros((steps,3))
    r_2 = np.zeros((steps,3))
    r_3 = np.zeros((steps,3))
    v_1 = np.zeros((steps,3))
    v_2 = np.zeros((steps,3))
    v_3 = np.zeros((steps,3))

    # Enter initial conditions
    r_1[0,0:3] = r_0[0,:]
    r_2[0,0:3] = r_0[1,:]
    r_3[0,0:3] = r_0[2,:]
    v_1[0,0:3] = v_0[0,:]
    v_2[0,0:3] = v_0[1,:]
    v_3[0,0:3] = v_0[2,:]

    m_1 = masses[0]
    m_2 = masses[1]
    m_3 = masses[2]

    # Leapfrog Method

    
    return r_1, v_1, r_2, v_2, r_3, v_3
def distance(r_ij):
    distance = np.sqrt( (r_ij[0])**2 + (r_ij[1])**2 + (r_ij[2])**2)
    return distance