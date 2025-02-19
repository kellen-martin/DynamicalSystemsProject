# Usefull functions for simulating 2 and 3 body problem
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

####################### Integrators ##################################
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
         v2_temp = v2_temp + lambda_*h*a2
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

         
    print('Simulation End')
    return r_1, v_1, r_2, v_2, r_3, v_3
######################################################################

###################### Tools #########################################
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

# Calculates the total kinetic energy of the 3BP
def kinetic_energy(v_1,v_2,v_3, m_1, m_2, m_3):
     # KE = sum(1/2mv^2)
     KE = .5*np.linalg.norm(v_1)**2*m_1 + .5*np.linalg.norm(v_2)**2*m_2 + .5*np.linalg.norm(v_3)**2*m_3
     return KE

def potnetial_energy(r_1, r_2, r_3, m_1, m_2, m_3, G):
        r_12 = r_1 - r_2
        r_13 = r_1 - r_3
        r_23 = r_2 - r_3
     
        U_12 = -G*(m_1*m_2)/distance(r_12);
        U_13 = -G*(m_1*m_3)/distance(r_13);
        U_23 = -G*(m_2*m_3)/distance(r_23);

        U = U_12 + U_13 + U_23;

        return U

def total_energy(r_1, v_1, r_2, v_2, r_3, v_3, masses, G):
        E = np.zeros((len(r_1)))
        m_1 = masses[0]
        m_2 = masses[1]
        m_3 = masses[2]

        for i in range(0,len(r_1)):
             KE = kinetic_energy(v_1[i,:],v_2[i,:],v_3[i,:], m_1, m_2, m_3)
             U = potnetial_energy(r_1[i,:], r_2[i,:], r_3[i,:], m_1, m_2, m_3, G)
             E[i] = KE + U

        
        return E

def rel_error(E):
     rel_error = abs(E[0] - E)/abs(E[0])
     return rel_error

def max_rel_error(E):
     mre = max(abs(E[0] - E)/abs(E[0]))
     return mre 

def absolute_error(E):
     abs_error = abs(E[0] - E)
     return abs_error

def max_abs_error(E):
     mae = max(abs(E[0] - E))
     return mae 
#####################################################################

########################### Plots ###################################
def plot_3BP(r_1,r_2,r_3):
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(r_1[:,0],r_1[:,1],r_1[:,2], label = 'Body 1')
    ax.plot(r_2[:,0],r_2[:,1],r_2[:,2], label = 'Body 2')
    ax.plot(r_3[:,0],r_3[:,1],r_3[:,2], label = 'Body 3')

    # plot starting points
    ax.plot(r_1[0,0],r_1[0,1],r_1[0,2],'*')
    ax.plot(r_2[0,0],r_2[0,1],r_2[0,2],'*')
    ax.plot(r_3[0,0],r_3[0,1],r_3[0,2],'*')

    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.title('Earth Fly-by')
    plt.legend()
    plt.show()
    return

def plot_difference(r_1, r_2, t_0, t_f):
     diff = np.zeros((len(r_1)))
     for i in range(0,len(diff)):
        diff[i] = np.linalg.norm(abs(r_1[i,:] - r_2[i,:]))

     times = np.linspace(t_0, t_f, len(diff))
     plt.plot(times, diff)
     plt.show()
     return

def plot_energy(E, t_0, t_f):
     times = np.linspace(t_0, t_f, len(E))

     plt.plot(times, E)
     plt.ylabel("Energy [J]")
     plt.xlabel("Time [s]")
     plt.title("Total Energy vs Time")

     plt.show()

def plot_energies(E1, E2, E3, t_0, t_f):
     times1 = np.linspace(t_0, t_f, len(E1))
     times2 = np.linspace(t_0, t_f, len(E2))
     times3 = np.linspace(t_0, t_f, len(E3))

     plt.plot(times1, E1, label = 'unperturbed')
     plt.plot(times2, E2, label = 'perturbation 1')
     plt.plot(times3, E3, label = 'perturbation 2')

     plt.ylabel("Energy [J]")
     plt.xlabel("Time [s]")
     plt.title("Total Energy vs Time")

     plt.legend()
     plt.show()
     return

def plot_relerrors(re1, re2, re3, t_0, t_f):
     times1 = np.linspace(t_0, t_f, len(re1))
     times2 = np.linspace(t_0, t_f, len(re2))
     times3 = np.linspace(t_0, t_f, len(re3))

     plt.plot(times1, re1, label = 'h = 1 day')
     plt.plot(times2, re2, label = 'h = 10000 seconds')
     plt.plot(times3, re3, label = 'h = 100 seconds')

     plt.ylabel("Relative Error")
     plt.xlabel("Time [s]")
     plt.title("Relative Error vs Time")

     plt.legend()
     plt.show()
     return

def plot_abserrors(ae1, ae2, ae3, t_0, t_f):
     times1 = np.linspace(t_0, t_f, len(ae1))
     times2 = np.linspace(t_0, t_f, len(ae2))
     times3 = np.linspace(t_0, t_f, len(ae3))

     plt.plot(times1, ae1, label = 'h = 1 day')
     plt.plot(times2, ae2, label = 'h = 10000 seconds')
     plt.plot(times3, ae3, label = 'h = 100 seconds')

     plt.ylabel("Absolute Error")
     plt.xlabel("Time [s]")
     plt.title("Absolute Error vs Time")

     plt.legend()
     plt.show()
     return

def plot_max_abserror(ae1, ae2, ae3, h1, h2, h3):
     aes = np.array([ae3, ae2, ae1])
     hs = np.array([h3, h2, h1])

     plt.plot(hs, aes)
     plt.xlabel('Step Size [s]')
     plt.ylabel('Absolute Error')
     plt.title('Absolute Error vs Step Size')
     plt.show()
     return
def plot_max_relerror(re1, re2, re3, h1, h2, h3):
     res = np.array([re3, re2, re1])
     hs = np.array([h3, h2, h1])

     plt.plot(hs, res)
     plt.xlabel('Step Size [s]')
     plt.ylabel('Relative Error')
     plt.title('Relative Error vs Step Size')
     plt.show()
     return

def plot_error(rel_error, t_0, t_f):
     times = np.linspace(t_0, t_f, len(rel_error))

     plt.plot(times, rel_error)
     plt.ylabel("Relative Error")
     plt.xlabel("Time [s]")
     plt.title("Relative Error vs Time")

     plt.show()
     return

def prep_animation(r_1, r_2, r_3, inte):
     r1_ani = r_1[::inte]
     r_2ani = r_2[::inte]
     r_3ani = r_3[::inte]
     return r1_ani, r_2ani, r_3ani

def animate_3BP(r_1, r_2, r_3,file_name):
    num_points = len(r_1)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_axis_off()

    def Animate_Func(num):
         ax.clear()
         ax.set_axis_off()
         ax.scatter(r_1[num,0], r_1[num,1], r_1[num,2])
         ax.scatter(r_2[num,0], r_2[num,1], r_2[num,2])
         ax.scatter(r_3[num,0], r_3[num,1], r_3[num,2])
         ax.plot3D(r_1[0,0], r_1[0,1],r_1[0,2])

         ax.set_xlim3d([min(r_3[:,0])-1, max(r_3[:,0])+1])
         ax.set_ylim3d([min(r_3[:,1])-1, max(r_3[:,1])+1])
         ax.set_zlim3d([min(r_3[:,2])-1, max(r_3[:,2])+1])


    line_ani = animation.FuncAnimation(fig, Animate_Func, interval=100, frames=num_points)
    plt.show()

    f = r"C:/Users/kmartin6/Desktop/Orbital Mechanics"
    f = f + file_name
    writergif = animation.PillowWriter(fps=num_points/6)
    line_ani.save(f, writer=writergif)
         
    
    return