import numpy as np

class solverParameters(object):
	"""docstring for solverParameters"""
	def __init__(self, mesh, delT=0.0001, rho=1000, nu=1, minL2Norm=1e-5, gx=0, gy=-9.81):
		super(solverParameters, self).__init__()
		self.mesh = mesh
		self.delT = delT
		self.rho = rho
		self.minL2Norm = minL2Norm
		self.gx = gx
		self.gy = gy
		self.nu = nu
		self.k = ((self.rho/self.delT)*(self.mesh.delX*self.mesh.delY))/(2*self.mesh.delXY)