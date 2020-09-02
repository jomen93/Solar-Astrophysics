import numpy as np
from fluid_solver import fluid_solver_2D
from fluid_diff import Fluid
import matplotlib.pyplot as plt 



Nx = 256
Ny = 256
Nz = 0

xmin = 0.; xmax = 1.0
ymin = 0.; ymax = 1.0
zmin = 0.; zmax = 0.0

dx = (xmax - xmin)/(Nx-1)
dy = (ymax - ymin)/(Ny-1)
dz = (zmax - zmin)/(Nz-1)
 
fluid = fluid_solver_2D(Nx,Ny,Nz,dx,dy,dz)

fluid.start = 1
fluid.S = 1e-1
fluid.iterations = int(1e2)
fluid.rho = 25./(36*np.pi)
fluid.mu = 1.
fluid.nu = 0.
fluid.eta = 0.
fluid.dt = 1e-3
fluid.nt = int(1.05/fluid.dt)
switch = False

x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
y,x = np.meshgrid(y,x)

xx,yy = np.mgrid[1:Nx-1,1:Ny-1]

u = v = p = src = Bx = By = A = np.zeros((Nx,Ny))

# Intital conditions

u = np.sin((y-ymax/2.))  #-np.exp(-((x-xmax/2.)/1.)**2)
v = -np.sin((x-xmax/2.)) #np.exp(-((y-xmax/2.)/1.)**2)
Bo = 1.#/np.sqrt(4.*np.pi)
A =  Bo*np.cos(4*np.pi*x)/(4*np.pi) + Bo*np.cos(4*np.pi*y)/(4*np.pi)


#Boundary conditions 

u[0,:] = u[-2, :]
u[-1, :] = u[1,:]
u[:,0] = u[:,-2]
u[:,-1] = u[:,1]

v[0,:] = v[-2, :]
v[-1, :] = v[1,:]
v[:,0] = v[:,-2]
v[:,-1] = v[:,1]

A[0,:] = A[-2, :]
A[-1, :] = A[1,:]
A[:,0] = A[:,-2]
A[:,-1] = A[:,1]


Bx[1:-1,1:-1] = (A[1:-1,2:]-A[1:-1,:-2])/dy
By[1:-1,1:-1] = (A[2:,1:-1]-A[:-2,1:-1])/dx

Bx[0,:] = Bx[-2, :]
Bx[-1, :] = Bx[1,:]
Bx[:,0] = Bx[:,-2]
Bx[:,-1] = Bx[:,1]

By[0,:] = By[-2, :]
By[-1, :] = By[1,:]
By[:,0] = By[:,-2]
By[:,-1] = By[:,1]

M = np.hypot(Bx,By)
plt.quiver(x, y, Bx,By, M , cmap=plt.cm.jet, width=0.01)
plt.title("Campo Magnetico inicial")
plt.colorbar()
plt.xlabel("$u_{x}$")
plt.ylabel("$u_{y}$")
plt.savefig("Campo_Mag_Inicial")
plt.show()


for i in range(fluid.start,fluid.nt+1):

	u_old = u.copy()
	u[1:-1,1:-1] = (fluid.nonlinear_advect_implicit_periodic_2d(u,u,v,xx,yy)- fluid.nonlinear_advect_implicit_periodic_2d(Bx,Bx,By,xx,yy)+ Bx[1:-1,1:-1] - fluid.apply_pressure_2dX(Bx**2+By**2,0.5))
	v[1:-1,1:-1] = (fluid.nonlinear_advect_implicit_periodic_2d(v,u_old,v,xx,yy)-fluid.nonlinear_advect_implicit_periodic_2d(By,Bx,By,xx,yy)+ By[1:-1,1:-1]-fluid.apply_pressure_2dY(Bx**2+By**2,0.5))

	u[0,:] = u[-2, :]; u[-1,:]=u[1, :]
	u[ :,0] = u[ :,-2]; u[ :,-1] = u[:, 1]

	v[0,:] = v[-2, :]; v[-1,:]=v[1, :]
	v[ :,0] = v[ :,-2]; v[ :,-1] = v[:, 1]

	src[1:-1,1:-1] = fluid.calc_source_2d(u,v)
	p = fluid.transform_pressure_poisson_2d(p,src)

	p[0,:] = p[-2, :]; p[-1,:]=p[1, :]
	p[ :,0] = p[ :,-2]; p[ :,-1] = p[:, 1]

	u[1:-1,1:-1] -= fluid.apply_pressure_2dX(p,1./fluid.rho)
	v[1:-1,1:-1] -= fluid.apply_pressure_2dY(p,1./fluid.rho)

	
	u[0, :] = u[-2, : ]; u[-1, :] = u[1, : ]
	u[:,0] = u[: , -2]; u[ : , -1] = u [:,1]

	v[0, :] = v[-2, : ]; v[-1, :] = v[1, : ]
	v[:,0] = v[: , -2]; v[ : , -1] = v [:,1]

	A[1:-1,1:-1] = fluid.nonlinear_advect_implicit_periodic_2d(A,u,v,xx,yy)
	A[0, :] = A[-2, : ]; A[-1, :] = A[1, : ]
	A[:,0] = A[: , -2]; A[ : , -1] = A [:,1]

	Bx[1:-1,1:-1] = ((A[1:-1, 2:]- A[1:-1, :-2])/fluid.dy)
	By[1:-1,1:-1] = ((A[2: , 1:-1]- A[:-2, 1:-1])/fluid.dx)
	

	Bx[0, :] = Bx[-2, : ]; Bx[-1, :] = Bx[1, : ]
	Bx[:,0] = Bx[: , -2]; Bx[ : , -1] = Bx [:,1]

	By[0, :] = By[-2, : ]; By[-1, :] = By[1, : ]
	By[:,0] = By[: , -2]; By[ : , -1] = By [:,1]

	if (i % 10 == 5 ) and switch == True:
		M = np.hypot(u,v)
		#plt.streamplot(x,y,u,v,color = "k",linewidth=0.8,density=0.4, arrowstyle='->', arrowsize=1.5)
		plt.quiver(x, y, u,v, M , cmap=plt.cm.jet, width=0.01)#,scale=10)
		plt.title("Campo de velocidades ")
		#plt.savefig("campoFinal"+str(i))
		plt.xlabel("$u_{x}$")
		plt.ylabel("$u_{y}$")
		plt.colorbar()
		plt.show()

	print(i)




M = np.hypot(Bx,By)
#plt.streamplot(x,y,u,v,color = "k",linewidth=0.8,density=0.4, arrowstyle='->', arrowsize=1.5)
plt.quiver(x, y, Bx,By, M , cmap=plt.cm.jet, width=0.01)#,scale=10)
plt.title("Campo Magnetico final")
plt.savefig("Campo_Mag_Final")
plt.colorbar()
plt.xlabel("$u_{x}$")
plt.ylabel("$u_{y}$")
plt.show()













