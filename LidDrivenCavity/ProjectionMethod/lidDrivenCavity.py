import numpy as np
from timeit import default_timer as timer
from mesh import Mesh
from parameters import solverParameters
from velocity import uVelocity, vVelocity
from pressure import Pressure
from resultPlot import results

np.set_printoptions(precision=4)
Nx = 30
Ny = 30
l2Norm = 1
itrOuter = 0
minL2Norm = 1e-5

mesh = Mesh(Nx,Ny)
sP = solverParameters(mesh=mesh, minL2Norm=minL2Norm)

u = uVelocity(mesh=mesh, solverParameters=sP)
uOld = uVelocity(mesh=mesh, solverParameters=sP)

v = vVelocity(mesh=mesh, solverParameters=sP)
vOld = vVelocity(mesh=mesh, solverParameters=sP)

p = Pressure(mesh=mesh, solverParameters=sP)

errorTracker = []
startTimeTotal = timer()
while (l2Norm > sP.minL2Norm):

	startTimeOuter = timer()
	itrOuter += 1
	print('\nIteration No.: '+str(itrOuter))
	u.values[0,:] = 0.0 
	v.values[0,:] = 0.0      
	u.values[-1,:] = 1
	v.values[-1,:] = 0.0
	u.values[:,0] = 0.0
	v.values[:,0] = 0.0
	u.values[:,-1] = 0.0
	v.values[:,-1] = 0.0

	uStar = uVelocity(mesh=mesh, solverParameters=sP)
	vStar = vVelocity(mesh=mesh, solverParameters=sP)

	print('Predicting Velocities')
	uStar.values[1:-1,1:-1] = u.values[1:-1,1:-1] + ((sP.delT)*(u.xMomentumTerm(v)))
	vStar.values[1:-1,1:-1] = v.values[1:-1,1:-1] + ((sP.delT)*(v.yMomentumTerm(u)))

	print('Solving for Pressure')
	startTimePressure = timer()
	p.values, itr = p.solve(uStar,vStar)
	totalTimePressure = (timer() - startTimePressure)
	print('Completed '+str(itr)+' iterations for solving pressure in '+str(totalTimePressure)+' s')

	print('Correcting Velocities')
	u.values[1:-1,1:-1] = uStar.values[1:-1,1:-1] - ((sP.delT/sP.rho)*(p.dpdx()))
	v.values[1:-1,1:-1] = vStar.values[1:-1,1:-1] - ((sP.delT/sP.rho)*(p.dpdy()))

	l2Norm = np.sqrt(np.sum(np.square(u.values -  uOld.values)) + np.sum(np.square(v.values - vOld.values)))
	print('Error in Velocties is: ',l2Norm)
	errorTracker.append(l2Norm)
	uOld.values = np.copy(u.values)
	vOld .values= np.copy(v.values)
	print('\nSimulation Time: '+str(itrOuter*sP.delT)+' s')
	totalTimeOuter = (timer() - startTimeOuter)
	print('Iteration Execution Time: '+str(totalTimeOuter)+' s')

np.savetxt('uVelocity.csv',u.values,delimiter=",")
np.savetxt('vVelocity.csv',v.values,delimiter=",")
np.savetxt('pressure.csv',p.values,delimiter=",")
np.savetxt('error.csv',np.asarray(errorTracker),delimiter=",")

totalTime = (timer() - startTimeTotal)
print('\nTotal Simulation Time: '+str(itrOuter*sP.delT)+' s')
print('Total Execution Time: '+str(totalTime)+' s\n')

results(mesh, itrOuter, u ,v, p, errorTracker)

# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# import matplotlib.cbook as cbook
# import pandas as pd

# file_location = 'uVelocity.csv'
# df = pd.read_csv(file_location)
# u.values = np.array(df)

# file_location = 'vVelocity.csv'
# df = pd.read_csv(file_location)
# v.values = np.array(df)

# file_location = 'error.csv'
# df = pd.read_csv(file_location)
# errorTracker = np.array(df)
# itrOuter = 1411
# X, Y = mesh.getMeshGrid()

# uNodes = u.ugetAvgNodalValues()
# vNodes = v.vgetAvgNodalValues()
# vel = np.sqrt(np.square(uNodes)+np.square(vNodes))

# fig, ax = plt.subplots(1, 1)
# pcm = ax.pcolor(X, Y, vel,norm=colors.LogNorm(vmin=vel.min(), vmax=vel.max()),cmap='RdBu_r')
# fig.colorbar(pcm, ax=ax, extend='max')
# plt.quiver(X, Y, u.values, v.values)

# plt.figure(2)
# plt.plot(np.linspace(1,itrOuter,itrOuter),np.asarray(errorTracker))
# plt.show()
