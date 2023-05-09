import numpy as np
from scipy.optimize import newton

# given initial conditions/force fxn, return 2 ellipses
# a = major axes
# bigC = centers
# b = minor axes
# e = eccentric axes
def exact_orbit(r,v,m,mu):
    c = np.cross(r,v)
    e = np.empty((2,3))
    for i in range(2):
       e[i] = np.cross(v[i],c[i])/mu[i]-r[i]/np.linalg.norm(r[i])

    a = np.empty(2)
    bigC = np.empty(2)
    b = np.empty(2)

    for i in range(2):
        a[i] = (-1)**i*np.linalg.norm(c[i])**2/(mu[i]*abs(np.linalg.norm(e[i])**2-1))
        bigC[i] = np.linalg.norm(e[i])*a[i]
        b[i] = (-1)**i*np.sqrt(a[i]**2-bigC[i]**2)

    return[e,a,b]

# Estimate angle u (eccentric anomaly) for a given t
# Using scipy newton, but I can find my old code as well..
def t_to_u(t,n,T,e):
    f = lambda u: n*(t-T)-u+np.linalg.norm(e[0])*np.sin(u)
    fp = lambda u: np.linalg.norm(e[0])*np.cos(u)-1
    x0 = n*(t-T)
    tol = 1e-8
    Nmax = 1000
    u = newton(f,x0,fp,rtol=tol,maxiter=Nmax)

    return(u)

def u_to_t(u,n,T,e):
    t = (u-np.linalg.norm(e[0])*np.sin(u))/n+T

    return(t)

#def exact_pos()

# RK4
def RK4_step(r,v,h,f,mu):
    [k1r,k1v] = f(r,v,mu)
    [k2r,k2v] = f(r+h/2*k1r,v+h/2*k1v,mu)
    [k3r,k3v] = f(r+h/2*k2r,v+h/2*k2v,mu)
    [k4r,k4v] = f(r+h*k3r,v+h*k3v,mu)

    rnew = r+(h/6)*(k1r+2*k2r+2*k3r+k4r)
    vnew = v+(h/6)*(k1v+2*k2v+2*k3v+k4v)

    return[rnew,vnew]



# Euler
def euler_step(r,v,h,f,mu):
    [fr,fv] = f(r,v,mu)
    rnew = r+fr*h
    vnew = v+fv*h
    #print(rnew)
    return [rnew,vnew]

# Leapfrog
def leapfrog(r,v,h,mu):


    for i in range(1,N+1):
        fi = f(t[i-1],w[i-1])  # i = i+1
        vhalf = w[i-1,2:] +.5*h*fi[2:]
        w[i,0:2] = w[i-1,0:2] + h*vhalf   # update r_i+1
        fi = f(t[i],w[i])
        w[i,2:] =vhalf+.5*h*fi[2:]   # update v_i+1

    return[t,w]

def leap_step(r,v,h,f,mu):
    k1v = h/2*f(r, v, mu)[1]
    vhalf = v + k1v

    # update position by full step
    k1r = h*vhalf
    rnew = r + k1r

    # update velocity by half step
    k2v = h/2*f(rnew, vhalf, mu)[1]
    vnew = vhalf + k2v

    return[rnew,vnew]
