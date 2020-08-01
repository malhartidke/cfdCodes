import numpy as np

class Pressure(object):
	"""docstring for Pressure"""
	def __init__(self, mesh, solverParameters, minError=1e-5, solveMethod='Jacobi'):
		super(Pressure, self).__init__()
		self.mesh = mesh
		self.solverParameters = solverParameters
		self.values = np.ones((self.mesh.Ny,self.mesh.Nx))
		self.solveMethod = solveMethod
		self.minError = minError

	def __str__():
		return f"Pressure at nodes in a Staggered Grid"

	def dpdx(self):
		return (self.values[1:-1,2:] - self.values[1:-1,:-2])/self.mesh.delX

	def dpdy(self):
		return (self.values[2:,1:-1] - self.values[:-2,1:-1])/self.mesh.delY

	def boundary(self):
		pass
		# return self.values

	def solve(self, u, v):

		error = 1
		itr = 0

		while error > self.minError:

			itr += 1
			pOld = np.copy(self.values)

			if self.solveMethod == 'Jacobi':
				self.values[1:-1,1:-1] = (self.mesh.a*(self.values[1:-1,2:] + self.values[1:-1,:-2])) + (self.mesh.b*(self.values[:-2,1:-1,]+self.values[2:,1:-1,])) - (self.solverParameters.k*((self.mesh.delX*(u.values[1:-1,1:-1]-u.values[1:-1,:-2]))+(self.mesh.delY*(v.values[1:-1,1:-1]-v.values[:-2,1:-1]))))
			elif self.solveMethod == 'GaussSeidel':
				for i in range(1,self.mesh.Ny-1):
					for j in range(1,self.mesh.Nx-1):
						self.values[i,j] = (self.mesh.a*(self.values[i+1,j] + self.values[i-1,j])) + (self.mesh.b*(self.values[i,j+1]+self.values[i,j-1])) - (self.solverParameters.k*((self.mesh.delX*(u.values[i,j]-u.values[i-1,j]))+(self.mesh.delY*(v.values[i,j]-v.values[i,j-1]))))
			else:
				raise Exception("Invalid Method to Solve for Pressure")

			self.values[:, -1] = self.values[:, -2] 
			self.values[0, :] = self.values[1, :]   
			self.values[:, 0] = self.values[:, 1]
			self.values[-1, :] = 0
			# self.values = self.boundary()
			error = np.linalg.norm(self.values - pOld, 2)

		return self.values, itr