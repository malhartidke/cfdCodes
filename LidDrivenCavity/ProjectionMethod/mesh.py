import numpy as np

class Mesh(object):
	"""docstring for Mesh"""
	def __init__(self, Nx=30, Ny=30):
		super(Mesh, self).__init__()
		self.Nx = Nx
		self.Ny = Ny
		self.delX = 1/(Nx-1)
		self.delY = 1/(Ny-1)
		self.delXY = (self.delX**2) + (self.delY**2)
		self.a = (self.delY**2)/(2*self.delXY)
		self.b = (self.delX**2)/(2*self.delXY)

	def getMeshGrid(self):
		X, Y = np.meshgrid(np.linspace(0,1,self.Nx),np.linspace(0,1,self.Ny))
		return X, Y