from datetime import datetime
import os
import psutil
process = psutil.Process(os.getpid())
startTime = datetime.now()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interpn
from scipy.special import erfi

def readshit(chain, file_num, columns):
	data = []
	for num in xrange(file_num):
		filename = "data/" + chain[0] + "/chains/" + chain[1] + "_" + str(num+1) + ".txt"
		with open(filename) as f:
			for read_line in f:
				line = filter(None, read_line.replace('\n','').split(' '))
				repeat = int(float(line[0]))
				for i in xrange(repeat):
					data.append([float(line[j]) for j in columns])
	return np.array(data)

print "Reading CMB"
CMB = readshit(['CMB', 'CMB'], 6, [2,3,4,6,7])
print "Number of CMB samples = ",len(CMB)
print "Reading LSS"
LSS = readshit(['Strong', 'Strong_L'], 6, [2,3,4,5,6])
print "Number of LSS samples = ",len(LSS)

print "Calculating difference vector"
diff = np.zeros([len(CMB)*len(LSS),5])
for i in xrange(len(LSS)):
	diff[i*len(CMB):(i+1)*len(CMB)] = CMB-LSS[i]

print "Histogramming"
bin_num = 2
hist, edges = np.histogramdd(diff,bins=bin_num)

ranges = np.array([[edges[j][i]+(edges[j][i+1]-edges[j][i])/2 for i in xrange(bin_num)] for j in xrange(5)])
norm = (len(CMB)*len(LSS))

print "Finding number of samples at 0"
val = interpn(ranges, hist, [0, 0, 0, 0, 0])
print "Number of samples at 0 = ", val
print "Tension = ", np.sqrt(2)*erfi(np.sum(hist[hist>=val])/norm)

print datetime.now() - startTime, ' seconds.'
total_mem = process.memory_info().rss >> 30
print total_mem, 'GB'

fig = plt.figure(figsize=(7, 6))
gs = GridSpec(4, 4)
gs.update(left=0.1, right=0.95, wspace=0, hspace=0)
labels = ['$\Omega_{\\rm b}h^2$', '$\Omega_{\\rm c}h^2$', '$\Theta_{\\rm MC}$', '$A_{\\rm s}$', '$n_{\\rm s}$'] 
for i in xrange(4):
	for j in xrange(i, 4):
		axes = tuple([k for k in xrange(5) if ((k != i) and (k != j+1))])
		ax = plt.subplot(gs[i,j])
		if i != j:
			ax.set_xticks([])
			ax.set_yticks([])
		else:
			ax.set_xlabel(labels[j+1],labelpad=0)
			ax.set_ylabel(labels[i],labelpad=0)
		
		plt.contourf(ranges[j+1],ranges[i],np.sum(hist, axis=axes))
		plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

plt.show()

