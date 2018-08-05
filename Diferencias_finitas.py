# wave.py
# Numerical solver for 1D wave equation
# Finite difference method
# u_tt = (c^2)*u_xx, 0<x<L
# Boundary conditions: u(0,t) = 0, u(L,t) = 0
# Initial conditions: u(x,0)=u0(x), u_t(x,0)=u_t0(x)
# 

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import animation

# Wave parameters
c = 10.0    # wave speed
time = 10.0 # length of time domain
L = 4.0     # length of space domain

# Finite Difference parameters
nx = 201            # discrete space intervals
sigma = 0.9        # Courant number to satisfy CFL condition
dx = L / nx         # space step size
dt = sigma*dx/c     # time step size
nt = int(time/dt)   # discrete time intervals

# Initial condition parameters
m = 1   # number of nodes in initial conditions
A = 1.0 # amplitude of initial condition

# returns initial value for u
def u0Init(nx,A,L,m):
    x = np.linspace(0,L,nx)
    stand = A*np.sin(m*np.pi*x/L) # standing wave
    gauss = A*np.exp(-10*(x-0.25*L)**2) # Gaussian
    #u0 = stand + gauss
    #u0 = stand + 0.7*A*np.sin((m+2)*np.pi*x/L) + 0.3*A*np.sin((m+4)*np.pi*x/L)
    u0 = A*np.exp(-10*(x-0.2*L)**2) - A*np.exp(-10*(x-0.6*L)**2)

    return u0

# returns initial value for u_t
def u_t0Init(nx,A,L,m):
    x = np.linspace(0,L,nx)
    u_t0 = np.zeros(nx) # no initial velocity
    #u_t0 = A*np.cos(m*np.pi*x/L)
    #u_t0 = 5*A*np.exp(-10*(x-0.3*L)**2)

    return u_t0


# Animations
def animate(data):
    x = np.linspace(0,L,nx)
    y = data
    line.set_data(x,y)
    return line,

# Finite difference solver
def wave1D(nt,dt,dx,c,L,A,m):
    un = np.zeros((nt,nx))      # initialize solution 
    u0 = u0Init(nx,A,L,m)       # get initial position
    u_t0 = u_t0Init(nx,A,L,m)   # get initial velocity
    rSq = (c*dt/dx)**2          # combine step sizes

    # first time step, to account for initial velocity
    un[0,:] = u0
    un[1,1:-2] = un[0,1:-2] - dt*u_t0[1:-2] + 0.5*rSq*(un[0,2:-1]-2*un[0,1:-2]+un[0,0:-3])
    un[1,0] = 0; un[1,-1] = 0; # Boundary conditions
    
    # rest of time steps
    for i in range(1,nt-1):
        un[i+1,1:-2] = 2*un[i,1:-2]-un[i-1,1:-2] + \
                rSq*(un[i,2:-1]-2*un[i,1:-2]+un[i,0:-3])

        un[i+1,0] = 0 # Boundary conditions
        un[i+1,-1] = 0

    return un


fig = plt.figure()
ax = plt.axes(xlim=(0,L),ylim=(-1.5*A,1.5*A))
line, = ax.plot([],[],lw=2)

un = wave1D(nt,dt,dx,c,L,A,m)

anim = animation.FuncAnimation(fig,animate,frames=un,interval=1)
plt.show()













