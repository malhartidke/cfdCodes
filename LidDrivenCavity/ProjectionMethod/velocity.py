import numpy as np

class uVelocity():
	"""docstring for uVelocity"""

	def __init__(self, mesh, solverParameters):
		super(uVelocity, self).__init__()
		self.mesh = mesh
		self.solverParameters = solverParameters
		self.values = np.zeros((self.mesh.Ny,self.mesh.Nx))

	def __str__(self):
		return f"Horizontal Velocity at edges in Staggered Grid"
		
	def dudx(self):
		return (self.values[1:-1,2:] - self.values[1:-1,:-2])/(2*self.mesh.delX)

	def dudy(self):
		return (self.values[2:,1:-1] - self.values[:-2,1:-1])/(2*self.mesh.delY)

	def duudxx(self):
		return (self.values[1:-1,2:] - (2*self.values[1:-1,1:-1]) + self.values[1:-1,:-2])/(self.mesh.delX**2)

	def duudyy(self):
		return (self.values[2:,1:-1] - (2*self.values[1:-1,1:-1]) + self.values[:-2,1:-1])/(self.mesh.delY**2)

	def xMomentumTerm(self, vValues):
		return - (self.values[1:-1,1:-1]*self.dudx()) - (((vValues.values[1:-1,2:]+vValues.values[1:-1,1:-1]+vValues.values[:-2,2:]+vValues.values[:-2,1:-1])/4)*self.dudy()) + (self.solverParameters.nu*(self.duudxx()+self.duudyy())) + self.solverParameters.gx

	def ugetAvgNodalValues(self):
		uNodes = (self.values[:-1,:] + self.values[1:,:]) / 2
		return (uNodes[:,:-1] + uNodes[:,1:]) / 2

class vVelocity():
	"""docstring for vVelocity"""
	def __init__(self, mesh, solverParameters):
		super(vVelocity, self).__init__()
		self.mesh = mesh
		self.solverParameters = solverParameters
		self.values = np.zeros((self.mesh.Ny,self.mesh.Nx))
	
	def __str__(self):
		return f"Vertical Velocity at edges in Staggered Grid"
		
	def dvdx(self):
		return (self.values[1:-1,2:] - self.values[1:-1,:-2])/(2*self.mesh.delX)

	def dvdy(self):
		return (self.values[2:,1:-1] - self.values[:-2,1:-1])/(2*self.mesh.delY)

	def dvvdxx(self):
		return (self.values[1:-1,2:]-(2*self.values[1:-1,1:-1])+self.values[1:-1,:-2])/(self.mesh.delX**2)

	def dvvdyy(self):
		return (self.values[2:,1:-1]-(2*self.values[1:-1,1:-1])+self.values[:-2,1:-1])/(self.mesh.delY**2)

	def yMomentumTerm(self, uValues):
		return - (((uValues.values[2:,1:-1]+uValues.values[1:-1,1:-1]+uValues.values[2:,:-2]+uValues.values[1:-1,:-2])/4)*self.dvdx()) - (self.values[1:-1,1:-1]*self.dvdy()) + (self.solverParameters.nu*(self.dvvdxx()+self.dvvdyy())) + self.solverParameters.gy

	def vgetAvgNodalValues(self):
		vNodes = (self.values[:-1,:] + self.values[1:,:]) / 2
		return (vNodes[:,:-1] + vNodes[:,1:]) / 2
