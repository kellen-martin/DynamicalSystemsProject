# Usefull functions for simulating 2 and 3 body problem
import numpy as np
import math
import matplotlib.pyplot as plt

# Euler Time Stepping 3-body problem
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
        a_1, a_2, a_3 = get_acceleration(r_1[i-1,:], r_2[i-1,:], r_3[i-1,:], m_1, m_2, m_3, G)

        # Update Vectors
        v_1[i,:] = v_1[i-1,:] + h*a_1
        v_2[i,:] = v_2[i-1,:] + h*a_2
        v_3[i,:] = v_3[i-1,:] + h*a_3

        r_1[i,:] = r_1[i-1,:] + h*v_1[i-1,:]
        r_2[i,:] = r_2[i-1,:] + h*v_2[i-1,:]
        r_3[i,:] = r_3[i-1,:] + h*v_3[i-1,:]

    return r_1, v_1, r_2, v_2, r_3, v_3

# Symplectic Leapfrog 3-body problem
def Symplectic_Leapfrog(r_0, v_0, masses, t_0, t_f, h):
    # Gravitational Constant
    G = 6.6743*10**(-11)# [m^3/kg*s^2]

    # Find number of steps
    steps = int((t_f + t_0)/h) + 1

    # Initialize Vector
    r_1 = np.zeros((steps,3))
    r_2 = np.zeros((steps,3))
    r_3 = np.zeros((steps,3))
    v_1 = np.zeros((steps,3))
    v_2 = np.zeros((steps,3))
    v_3 = np.zeros((steps,3))

    # Enter initial conditions
    r_1[0,:] = r_0[0,:]
    r_2[0,:] = r_0[1,:]
    r_3[0,:] = r_0[2,:]
    v_1[0,:] = v_0[0,:]
    v_2[0,:] = v_0[1,:]
    v_3[0,:] = v_0[2,:]

    m_1 = masses[0]
    m_2 = masses[1]
    m_3 = masses[2]

    # Leapfrog Method
    for i in range(0,steps-1):

        # Velocity half-step
        a1, a2, a3 = get_acceleration(r_1[i,:], r_2[i,:], r_3[i,:], m_1, m_2, m_3, G)
        v1_half = v_1[i,:] + .5*h*a1
        v2_half = v_2[i,:] + .5*h*a2
        v3_half = v_3[i,:] + .5*h*a3

        # Position full-step
        r_1[i+1,:] = r_1[i,:] + h*v1_half
        r_2[i+1,:] = r_2[i,:] + h*v2_half
        r_3[i+1,:] = r_3[i,:] + h*v3_half

        # Velocity full-step
        a1, a2, a3 = get_acceleration(r_1[i+1,:], r_2[i+1,:], r_3[i+1,:], m_1, m_2, m_3, G)
        v_1[i+1,:] = v1_half + .5*h*a1
        v_2[i+1,:] = v2_half + .5*h*a2
        v_3[i+1,:] = v3_half + .5*h*a3

    print('Simulation End')
    return r_1, v_1, r_2, v_2, r_3, v_3

def PEFRL_3BP(r_0, v_0, masses, t_0, t_f, h):
      # Gravitational Constant
    G = 6.6743*10**(-11)# [m^3/kg*s^2]

    # Find number of steps
    steps = int((t_f + t_0)/h) + 1

    # Initialize Vector
    r_1 = np.zeros((steps,3))
    r_2 = np.zeros((steps,3))
    r_3 = np.zeros((steps,3))
    v_1 = np.zeros((steps,3))
    v_2 = np.zeros((steps,3))
    v_3 = np.zeros((steps,3))

    # Enter initial conditions
    r_1[0,:] = r_0[0,:]
    r_2[0,:] = r_0[1,:]
    r_3[0,:] = r_0[2,:]
    v_1[0,:] = v_0[0,:]
    v_2[0,:] = v_0[1,:]
    v_3[0,:] = v_0[2,:]

    m_1 = masses[0]
    m_2 = masses[1]
    m_3 = masses[2]

    # PEFRL variables
    zeta = 0.1786178958448091
    lambda_ = -0.2123418310626054
    chi = -0.6626458266981849*10**(-1)

    # PEFRL time
    for i in range(0,steps-1):
         # Step 1
         r1_temp = r_1[i,:] + zeta*h*v_1[i,:]
         r2_temp = r_2[i,:] + zeta*h*v_2[i,:]
         r3_temp = r_3[i,:] + zeta*h*v_3[i,:]

         # Step 2
         a1, a2, a3 = get_acceleration(r1_temp, r2_temp, r3_temp, m_1, m_2, m_3, G)
         v1_temp = v_1[i,:] + (1 - 2*lambda_)*.5*h*a1
         v2_temp = v_2[i,:] + (1 - 2*lambda_)*.5*h*a2
         v3_temp = v_3[i,:] + (1 - 2*lambda_)*.5*h*a3

         # Step 3 
         r1_temp = r1_temp + chi*h*v1_temp
         r2_temp = r2_temp + chi*h*v2_temp
         r3_temp = r3_temp + chi*h*v3_temp

         # Step 4 
         a1, a2, a3 = get_acceleration(r1_temp, r2_temp, r3_temp, m_1, m_2, m_3, G)
         v1_temp = v1_temp + lambda_*h*a1
         v2_temp = v3_temp + lambda_*h*a2
         v3_temp = v3_temp + lambda_*h*a3

         # Step 5
         r1_temp = r1_temp + (1 - 2*(chi + zeta))*h*v1_temp
         r2_temp = r2_temp + (1 - 2*(chi + zeta))*h*v2_temp
         r3_temp = r3_temp + (1 - 2*(chi + zeta))*h*v3_temp

         # Step 6 
         a1, a2, a3 = get_acceleration(r1_temp, r2_temp, r3_temp, m_1, m_2, m_3, G)
         v1_temp = v1_temp + lambda_*h*a1
         v2_temp = v2_temp + lambda_*h*a2
         v3_temp = v3_temp + lambda_*h*a3

         # Step 7 
         r1_temp = r1_temp + chi*h*v1_temp
         r2_temp = r2_temp + chi*h*v2_temp
         r3_temp = r3_temp + chi*h*v3_temp

         # Step 8
         a1, a2, a3 = get_acceleration(r1_temp, r2_temp, r3_temp, m_1, m_2, m_3, G)
         v1_temp = v1_temp + (1 - 2*lambda_)*h*.5*a1
         v2_temp = v2_temp + (1 - 2*lambda_)*h*.5*a2
         v3_temp = v3_temp + (1 - 2*lambda_)*h*.5*a3

         # Step 9 
         r1_temp = r1_temp + zeta*h*v1_temp
         r2_temp = r2_temp + zeta*h*v2_temp
         r3_temp = r3_temp + zeta*h*v3_temp

         r_1[i+1,:] = r1_temp
         v_1[i+1,:] = v1_temp
         r_2[i+1,:] = r2_temp
         v_2[i+1,:] = v2_temp
         r_3[i+1,:] = r3_temp
         v_3[i+1,:] = v3_temp

         

    return r_1, v_1, r_2, v_2, r_3, v_3

def distance(r_ij):
    distance = np.sqrt( (r_ij[0])**2 + (r_ij[1])**2 + (r_ij[2])**2)
    return distance

# Calculates acceleration for each body 
def get_acceleration(r_1,r_2,r_3,m_1,m_2,m_3,G):
# Relative Positions
        r_12 = r_1 - r_2
        r_13 = r_1 - r_3
        r_23 = r_2 - r_3
        r_21 = -r_12
        r_31 = -r_13
        r_32 = -r_23

        # Calculate Accelerations
        a_1 = -G*m_2*(r_12)/(distance(r_12)**3) - G*m_3*(r_13)/(distance(r_13)**3)
        a_2 = -G*m_1*(r_21)/(distance(r_21)**3) - G*m_3*(r_23)/(distance(r_23)**3)
        a_3 = -G*m_1*(r_31)/(distance(r_31)**3) - G*m_2*(r_32)/(distance(r_32)**3)
    
        return a_1, a_2, a_3

# calculates center of mass
def get_com(r_1,r_2,r_3,m_1,m_2,m_3):
     com = (r_1*m_1 + r_2*m_2 + r_3*m_3)/(m_1 + m_2 + m_3)
     return com

