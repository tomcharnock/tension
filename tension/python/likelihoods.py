import numpy as np
import matplotlib.pyplot as plt
import sys

plot = False

CMB = np.loadtxt('data/CMB/chains/CMB_1.txt',usecols=(2,3,4,5,6))
LSS = np.loadtxt('data/Strong/chains/Strong_L_1.txt',usecols=(2,3,4,5,6))

dim = 5
num = 10

ranges = np.zeros([dim,2])
for i in xrange(dim):
	ranges[i,0] = min(np.min(CMB[:,i]),np.min(LSS[:,i]))
	ranges[i,1] = max(np.max(CMB[:,i]),np.max(LSS[:,i]))

CMBhist,bins = np.histogramdd(CMB[:,:dim],bins=num,range=ranges,normed=True)
LSShist,bins = np.histogramdd(LSS[:,:dim],bins=num,range=ranges,normed=True)

widths = np.zeros([dim])
binvolume = 1
for i in xrange(dim):
	widths[i] = bins[i][1]-bins[i][0]
	binvolume = binvolume*widths[i]

print np.sum(CMBhist*LSShist)*binvolume

if plot and dim == 2:
	X, Y = np.meshgrid(bins[0],bins[1])

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)
	ax1.pcolormesh(X, Y, CMBhist)
	ax1.pcolormesh(X, Y, LSShist)

	hist = CMBhist*LSShist

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	ax2.pcolormesh(X, Y, hist)

	plt.show()

