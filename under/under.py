from datetime import datetime
import os
import psutil
process = psutil.Process(os.getpid())
startTime = datetime.now()

import argparse
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interpn
from scipy.special import erfinv

def readshit(chain, file_num, columns):
	data = []
	for num in xrange(file_num):
		filename = "data/" + chain[0] + "/chains/" + chain[1] + "_" + str(num+1) + ".txt"
		num_lines = 0
		with open(filename) as f:
			for read_line in f:
				num_lines += 1
		skip_lines = int(num_lines/3.)
		num_lines = 0
		with open(filename) as f:
			for read_line in f:
				if num_lines > skip_lines:
					line = filter(None, read_line.replace('\n','').split(' '))
					repeat = int(float(line[0]))
					for i in xrange(repeat):
						data.append([float(line[j]) for j in columns])
				num_lines += 1
					
	return np.array(data)

def histogram(i):
	diff = CMB-LSS[i]
	res, edges = np.histogramdd(diff, bins=args.bins, range=bins)
	return res

def sigma(P, l):
	P_s = np.sort(P.ravel())[::-1]
	sigma = 0
	for i in P_s:
		sigma += i
		if sigma >= l:
			return i
	return 1	

parser = argparse.ArgumentParser()
parser.add_argument("bins", type=int)
parser.add_argument("CMB_folder", type=str)
parser.add_argument("CMB_data", type=str)
parser.add_argument("LSS_folder", type=str)
parser.add_argument("LSS_data", type=str)
args = parser.parse_args()

print "Reading CMB"
CMB = readshit([args.CMB_folder, args.CMB_data], 6, [2, 3, 4, 6, 7])
print "Number of CMB samples = ",len(CMB)
print "Reading LSS"
LSS = readshit([args.LSS_folder, args.LSS_data], 6, [2, 3, 4, 5, 6])
print "Number of LSS samples = ",len(LSS)

bins = np.array([[min(np.min(CMB[:, i]),np.min(LSS[:, i])), max(np.max(CMB[:, i]),np.max(LSS[:, i]))] for i in xrange(5)])

CMB_hist, CMB_edges = np.histogramdd(CMB, bins = args.bins)
CMB_ranges = np.array([[CMB_edges[j][i]+(CMB_edges[j][i+1]-CMB_edges[j][i])/2 for i in xrange(args.bins)] for j in xrange(5)])

LSS_hist, LSS_edges = np.histogramdd(LSS, bins = args.bins)
LSS_ranges = np.array([[LSS_edges[j][i]+(LSS_edges[j][i+1]-LSS_edges[j][i])/2 for i in xrange(args.bins)] for j in xrange(5)])

CMB_domainsize = 1
LSS_domainsize = 1
for i in xrange(5):
	CMB_domainsize = CMB_domainsize * (CMB_ranges[i, 1] - CMB_ranges[i, 0])
	LSS_domainsize = LSS_domainsize * (LSS_ranges[i, 1] - LSS_ranges[i, 0])

CMB_hist = CMB_hist/(len(CMB)*CMB_domainsize)
LSS_hist = LSS_hist/(len(LSS)*LSS_domainsize)

print datetime.now() - startTime, ' seconds.'
total_mem = process.memory_info().rss >> 20
print total_mem, 'MB'

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'legend.fontsize': 8, 'font.size': 8}
plt.rcParams.update(params)

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
			ax.set_xlabel(labels[j+1], labelpad = 0)
			ax.set_ylabel(labels[i], labelpad = 0)

		CMB_domains = (CMB_ranges[j+1, 1] - CMB_ranges[j+1, 0])*(CMB_ranges[i, 1] - CMB_ranges[i, 0])
		LSS_domains = (LSS_ranges[j+1, 1] - LSS_ranges[j+1, 0])*(LSS_ranges[i, 1] - LSS_ranges[i, 0])
		CMB_P = np.sum(CMB_hist*CMB_domainsize/CMB_domains, axis = axes)
		LSS_P = np.sum(LSS_hist*LSS_domainsize/LSS_domains, axis = axes)
		gridsize = 100
		x = np.linspace(bins[j+1, 0], bins[j+1, 1], gridsize)
		y = np.linspace(bins[i, 0], bins[i, 1], gridsize)
		interp_domains = (x[1] - x[0]) * (y[1] - y[0])
		CMB_plot = np.zeros([gridsize, gridsize])
		LSS_plot = np.zeros([gridsize, gridsize])
		for k in xrange(gridsize):
			for l in xrange(gridsize):
				CMB_plot[k, l] = interpn([CMB_ranges[j+1], CMB_ranges[i]], CMB_P, [x[k], y[l]], bounds_error = False, fill_value = 0)
				LSS_plot[k, l] = interpn([LSS_ranges[j+1], LSS_ranges[i]], LSS_P, [x[k], y[l]], bounds_error = False, fill_value = 0)
		CMB_levels = [sigma(CMB_plot, 0.997/interp_domains), sigma(CMB_plot, 0.95/interp_domains), sigma(CMB_plot, 0.68/interp_domains)]
		LSS_levels = [sigma(LSS_plot, 0.997/interp_domains), sigma(LSS_plot, 0.95/interp_domains), sigma(LSS_plot, 0.68/interp_domains)]
		ax.contour(x, y, CMB_plot, colors = 'blue', levels = CMB_levels)
		ax.contour(x, y, LSS_plot, colors = 'red', levels = LSS_levels)
		ax.locator_params(nbins=5)


#plt.savefig('difference_vector/plots/contours_' + str(args.bins) + '.pdf', bbox_inches = 'tight')
plt.show()

