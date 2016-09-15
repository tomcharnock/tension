from datetime import datetime
import os
import psutil
process = psutil.Process(os.getpid())
startTime = datetime.now()

import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d

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

def extrapolate(interpolator):
	xs = interpolator.x
	ys = interpolator.y
	
	def pointwise(x):
        	if x < xs[0]:
        		return 0
        	elif x > xs[-1]:
            		return 0
        	else:
        		return interpolator(x)

	def ufunclike(xs):
        	return np.array(map(pointwise, np.array(xs)))
	
	return ufunclike

def domainsize(ranges):
	domains = []
	for i in xrange(5):
		domains.append(ranges[i][1]-ranges[i][0])
	return domains

def sigma(data,domains):
	sorted_val = np.sort(data.ravel())[::-1]
	domain = 1
	for i in xrange(len(domainsize)):
		domain = domain*domains[i]
	sigma = 0
	levels = []
	one_sigma = True
	two_sigma = True
	three_sigma = True
	for i in sorted_val:
		sigma += i
		if (sigma > 0.68/domain) and (one_sigma):
			levels.append(i)
			one_sigma = False
		if (sigma > 0.95/domain) and (two_sigma):
			levels.append(i)
			two_sigma = False
		if (sigma > 0.997/domain) and (three_sigma):
			levels.append(i)
			three_sigma = False
		i += 1
	return levels[::-1]

def sigma_values(max_prob, min_prob, ranges):
	max_prob = max_prob*(ranges[1]-ranges[0])
	min_prob = min_prob*(ranges[1]-ranges[0])
	sorted_val = np.sort(min_prob.ravel())[::-1]
	sigma = 0
	for i in sorted_val:
		if np.sum(max_prob[min_prob>i]) > 0.68:
			return sigma
		sigma += i

	return 1

def likelihoods(CMB_hist, CMB_edges, LSS_hist, LSS_edges):

	CMB_domains = domainsize(CMB_ranges)
	LSS_domains = domainsize(LSS_ranges)

	params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','legend.fontsize':8,'font.size':8}
	plt.rcParams.update(params)
	fig = plt.figure(figsize=(8,2))
	gs = GridSpec(1,5)
	gs.update(left=0.05, right=0.95, bottom=0.25, wspace=0, hspace=0)
	labels = ['$\Omega_{\\rm b}h^2$', '$\Omega_{\\rm c}h^2$', '$\Theta_{\\rm MC}$', '$A_{\\rm s}$', '$n_{\\rm s}$'] 
	for i in xrange(5):
		axes = [k for k in xrange(5) if (k != i)]
		
		CMB_domain = 1
		LSS_domain = 1
		for j in axes:
			CMB_domain = CMB_domain*CMB_domains[j]
			LSS_domain = LSS_domain*LSS_domains[j]

		ax = plt.subplot(gs[i])
		ax.set_yticks([])
		ax.set_xlabel(labels[i],labelpad=0)

		ranges = np.linspace(min(CMB_ranges[i][0],LSS_ranges[i][0]), max(CMB_ranges[i][-1],LSS_ranges[i][-1]),100)

		CMB_plot = np.sum(CMB_hist*CMB_domain, axis=tuple(axes))
		LSS_plot = np.sum(LSS_hist*LSS_domain, axis=tuple(axes))
		if np.max(CMB_plot) >= np.max(LSS_plot):
			max_prob = CMB_plot
		        max_ranges = CMB_ranges[i]
			min_prob = LSS_plot
			min_ranges = LSS_ranges[i]
		else:
			max_prob = LSS_plot
			max_ranges = LSS_ranges[i]
			min_prob = CMB_plot
			min_ranges = CMB_ranges[i]

		max_spline = extrapolate(interp1d(max_ranges, max_prob))(ranges)
		ax.plot(ranges,max_spline, color='blue')

		min_spline = extrapolate(interp1d(min_ranges, min_prob))(ranges)
		ax.plot(ranges, min_spline, color='red')
		print np.max(np.sum(CMB_hist,axis=tuple(axes))),np.max(max_spline)
		print np.max(np.sum(LSS_hist,axis=tuple(axes))),np.max(min_spline)
		'''


		sigma_val = sigma_values(max_prob, min_prob, ranges)	
		sigma_levels = sigma(max_prob, [ranges[1]-ranges[0]])

		ax.fill_between(ranges, 0, max_prob, where = max_prob > sigma_levels[-1],color='blue',alpha=0.5)
		
		ax.text(0.05, 0.95, '{0:.3f}'.format(sigma_val), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes,color='black')
		'''
		ax.set_xticks(ax.get_xticks()[:-1:2])

	plt.savefig('overlap/plots/likelihoods_'+str(args.bins)+'.pdf')
	plt.show()

def contours():
	fig = plt.figure(figsize=(7, 6))
	gs = GridSpec(4, 4)
	gs.update(left=0.1, right=0.95, wspace=0, hspace=0)
	labels = ['$\Omega_{\\rm b}h^2$', '$\Omega_{\\rm c}h^2$', '$\Theta_{\\rm MC}$', '$A_{\\rm s}$', '$n_{\\rm s}$'] 
	for i in xrange(4):
		for j in xrange(i, 4):
			axes = [k for k in xrange(5) if ((k != i) and (k != j+1))]
			CMB_domain = 1
			LSS_domain = 1
			for k in axes:
				CMB_domain = CMB_domain*CMB_domainsize[k]
				LSS_domain = LSS_domain*LSS_domainsize[k]
				
			ax = plt.subplot(gs[i,j])
			if i != j:
				ax.set_xticks([])
				ax.set_yticks([])
			else:
				ax.set_xlabel(labels[j+1],labelpad=0)
				ax.set_ylabel(labels[i],labelpad=0)
			CMB_plot = np.sum(CMB_hist*CMB_domain, axis=tuple(axes))
			CMB_levels = sigma(CMB_plot, CMB_domainsize[i], CMB_domainsize[j+1])
			plt.contour(CMB_ranges[j+1],CMB_ranges[i],CMB_plot,levels = CMB_levels)
			LSS_plot = np.sum(LSS_hist*LSS_domain, axis=tuple(axes))
			LSS_levels = sigma(LSS_plot, LSS_domainsize[i], LSS_domainsize[j+1])
			plt.contour(LSS_ranges[j+1],LSS_ranges[i],LSS_plot,levels=LSS_levels)
			plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	plt.savefig('overlap/plots/contours_' + str(args.bins) + '.pdf')
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

CMB_hist, CMB_edges =  np.histogramdd(CMB, bins=args.bins)
CMB_ranges = np.array([[CMB_edges[j][i]+(CMB_edges[j][i+1]-CMB_edges[j][i])/2 for i in xrange(args.bins)] for j in xrange(5)])
CMB_hist = CMB_hist/np.sum(CMB_hist)
del CMB

LSS_hist, LSS_edges =  np.histogramdd(LSS, bins=args.bins)
LSS_ranges = np.array([[LSS_edges[j][i]+(LSS_edges[j][i+1]-LSS_edges[j][i])/2 for i in xrange(args.bins)] for j in xrange(5)])
LSS_hist = LSS_hist/(np.sum(LSS_hist))
del LSS

likelihoods(CMB_hist, CMB_edges, LSS_hist, LSS_hist)

print datetime.now() - startTime, ' seconds.'
total_mem = process.memory_info().rss >> 20
print total_mem, 'MB'

#contours()

