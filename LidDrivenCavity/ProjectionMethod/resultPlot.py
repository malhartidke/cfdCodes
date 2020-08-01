import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook

def results(mesh, itrOuter, u ,v, p, errorTracker):
	
	X, Y = mesh.getMeshGrid()
	uNodes = u.ugetAvgNodalValues()
	vNodes = v.vgetAvgNodalValues()
	vel = np.sqrt(np.square(uNodes)+np.square(vNodes))

	fig, ax = plt.subplots(1, 1)
	pcm = ax.pcolor(X, Y, vel,norm=colors.LogNorm(vmin=vel.min(), vmax=vel.max()),cmap='RdBu_r')
	fig.colorbar(pcm, ax=ax, extend='max')
	plt.quiver(X, Y, u.values, v.values)

	plt.figure(2)
	plt.plot(np.linspace(1,itrOuter,itrOuter),np.asarray(errorTracker))
	plt.show()