import sys
from fluid_diff import Fluid

class solver_2D(Fluid):

	def bd_1st_2d(self, f,c1,c2): # Back difference 
		g = c1*(f[1:-1,1:-1]-f[:-2,1:-1])/self.dx
		g += c2*(f[1:-1,1:-1]-f[1:-1,:-2])/self.dy
		return g

	def cdi_1st_2d(self,f,c1,c2): # Central difference implicit
		g = c1*(f[2: ,1:-1]+f[:-2,1:-1])/(2*self.dx)
		g += c2*(f[1:-1,2:]+f[1:-1, :2])/(2*self.dy)
		return g

	def cd_1st_2dX(self,f,c):# Central differences x
		g = c*(f[2: ,1:-1]-f[ :-2,1:-1])/(2*self.dx)
		return g

	def cd_1st_2dY(self,f,c): # Central differences y
		g = c*(f[1:-1,2:]-f[1:-1, :-2])/(2*self.dy)
		return g

	def cd_2nd_2d(self,f,c1,c2):
		g =(c1/self.dx**2 * (f[2:, 1:-1]-2*f[1:-1,1:-1] + f[ :-2,1:-1])+c2/self.dy**2 * (f[1:-1,2:]-2*f[1:-1,1:-1] + f[ 1:-1,:-2]))
		return g

