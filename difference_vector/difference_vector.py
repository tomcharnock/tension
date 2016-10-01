from datetime import datetime
import os
import psutil
process = psutil.Process(os.getpid())
startTime = datetime.now()

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interpn
from scipy.special import erfinv

def get_columns(root, chain, params):
	columns = []
	filename = root + chain + ".paramnames"
	with open(filename) as f:
		i = 0
		for read_line in f:
			line = filter(None, read_line.replace('\n','').split('\t'))
			if line[0] in params:
				columns.append(i+2)
			i += 1
	return columns

def get_labels(root, chain, params):
	label_dict = {'omegabh2': '$\Delta\Omega_{\\rm b}h^2$', 'omegach2': '$\Delta\Omega_{\\rm c}h^2$', 'theta': '$\Delta\Theta_{\\rm MC}$', 'logA': '$\Delta\log A_{\\rm s}$', 'ns': '$\Delta n_{\\rm s}$', 'mnu': '$\Delta\sum m_\\nu$', 'meffsterile': '$\Delta m_{\\rm eff}^{\\rm sterile}$', 'nnu': '$\Delta N_{\\rm eff}$'}
	labels = []
	filename = root + chain + ".paramnames"
	with open(filename) as f:
		i = 0
		for read_line in f:
			line = filter(None, read_line.replace('\n','').split('\t'))
			if line[0] in params:
				labels.append(label_dict[line[0]])
			i += 1
	return labels

def readshit(root, chain, file_num, columns):
	data = []
	for num in xrange(file_num):
		filename = root + chain + "_" + str(num+1) + ".txt"
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

def sci_lims(x, y):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot([x[0], x[-1]], [0, 0])
	ax.plot([0, 0], [y[0], y[-1]])
	ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	plt.tight_layout()
	xoffset = ax.get_xaxis().get_offset_text().get_text()
	xoffset = xoffset.replace('1e', '')
	xoffset = xoffset.replace(u'\u2212', '-')
	if xoffset == '':
		xoffset = '0'
	xoffset = float(xoffset)
	yoffset = ax.get_yaxis().get_offset_text().get_text()
	yoffset = yoffset.replace('1e', '')
	yoffset = yoffset.replace(u'\u2212', '-')
	if yoffset == '':
		yoffset = '0'
	yoffset = float(yoffset)
	plt.close()
	return [xoffset, yoffset]
	

parser = argparse.ArgumentParser()
parser.add_argument("-bins", type = int)
parser.add_argument("-params", nargs = '+')
parser.add_argument("-root", type = str)
parser.add_argument("-CMB", type = str)
parser.add_argument("-LSS", type = str)
args = parser.parse_args()

CMB_columns = get_columns(args.root, args.CMB, args.params)
LSS_columns = get_columns(args.root, args.LSS, args.params)
num_params = len(args.params)

print "Reading CMB"
CMB = readshit(args.root, args.CMB, 6, CMB_columns)
CMB = CMB[:100]
print "Number of CMB samples = ",len(CMB)
print "Reading LSS"
LSS = readshit(args.root, args.LSS, 6, LSS_columns)
LSS = LSS[:100]
print "Number of LSS samples = ",len(LSS)

bins = [[np.min(CMB[:, i])-np.max(LSS[:, i]), np.max(CMB[:, i])-np.min(LSS[:, i])] for i in xrange(num_params)]
a, edges = np.histogramdd([[0] for i in xrange(num_params)], bins=args.bins, range=bins)
ranges = np.array([[edges[j][i]+(edges[j][i+1]-edges[j][i])/2 for i in xrange(args.bins)] for j in xrange(num_params)])

sci = []
for i in xrange(num_params-1):
	for j in xrange(i, num_params-1):
		sci.append(sci_lims(ranges[j+1], ranges[i]))

domainsize = 1
for i in xrange(num_params):
	domainsize = domainsize * (ranges[i, 1] - ranges[i, 0])
	
norm = (len(CMB)*len(LSS))*domainsize

hist = np.zeros([args.bins, args.bins, args.bins, args.bins, args.bins])

print "Using ", multiprocessing.cpu_count(), " processors"
pool = Pool()
for i in pool.imap_unordered(histogram, xrange(len(LSS))):
	hist += i

hist = hist/norm
print "Finding number of samples at 0"
val = interpn(ranges, hist, [0 for i in xrange(num_params)], bounds_error = False, fill_value = 0)
C = np.sum(hist[hist >= val]*domainsize)
print "Number of samples at 0 = ", val
print "C = ",C
tension = np.sqrt(2)*erfinv(C)
print "Tension = ", tension

print datetime.now() - startTime, ' seconds.'
total_mem = process.memory_info().rss >> 20
print total_mem, 'MB'

plot_params = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'legend.fontsize': 8, 'font.size': 8}
plt.rcParams.update(plot_params)

fig = plt.figure(figsize=(1.4 * num_params, 1.2 * num_params))
gs = GridSpec(num_params-1, num_params-1)
gs.update(left=0.1, right=0.95, wspace=0, hspace=0)

labels = get_labels(args.root, args.CMB, args.params)
a = 0
for i in xrange(4):
	for j in xrange(i, 4):
		axes = tuple([k for k in xrange(5) if ((k != i) and (k != j+1))])
		ax = plt.subplot(gs[i,j])
		if i != j:
			ax.set_xticks([])
			ax.set_yticks([])
		else:
			ax.set_xlabel(labels[j+1] + ' [$\\times10^{'+ '{0:0d}'.format(int(sci[a][0])) + '}$]',labelpad=0)
			ax.set_ylabel(labels[i] + ' [$\\times10^{'+ '{0:0d}'.format(int(sci[a][0])) + '}$]',labelpad=0)

		domains = (ranges[j+1, 1] - ranges[j+1, 0])*(ranges[i, 1] - ranges[i, 0])
		P = np.sum(hist*domainsize/domains, axis = axes)
		interpolate = False
		if interpolate:
			gridsize = 100
			x = np.linspace(ranges[j+1, 0], ranges[j+1, -1], gridsize)
			y = np.linspace(ranges[i, 0], ranges[i, -1], gridsize)
			interp_domains = (x[1] - x[0]) * (y[1] - y[0])
			P_plot = np.zeros([gridsize, gridsize])
			for k in xrange(gridsize):
				for l in xrange(gridsize):
					P_plot[k, l] = 	interpn([ranges[j+1], ranges[i]], P, [x[k], y[l]], bounds_error = False, fill_value = 0)
			P_levels = [sigma(P_plot, 0.997/interp_domains), sigma(P_plot, 0.95/interp_domains), sigma(P_plot, 0.68/interp_domains)]
			ax.contour(x/(10**sci[a][0]), y/(10**sci[a][1]), P_plot, colors = 'purple', levels = P_levels)
			ax.plot([x[0]/(10**sci[a][0]), x[-1]/(10**sci[a][0])], [0, 0], color = 'black', linestyle = ':')
			ax.plot([0, 0], [y[0]/(10**sci[a][1]), y[-1]/(10**sci[a][1])], color = 'black', linestyle = ':')
		else:
			P_levels = [sigma(P, 0.997/domains), sigma(P, 0.95/domains), sigma(P, 0.68/domains)]
			ax.contour(ranges[j+1]/(10**sci[a][0]), ranges[i]/(10**sci[a][1]), P, colors = 'purple', levels = P_levels)
			ax.plot([ranges[j+1, 0]/(10**sci[a][0]), ranges[j+1, -1]/(10**sci[a][0])], [0, 0], color = 'black', linestyle = ':')
			ax.plot([0, 0], [ranges[i, 0]/(10**sci[a][1]), ranges[i, -1]/(10**sci[a][1])], color = 'black', linestyle = ':')
		ax.locator_params(nbins=5)

		a += 1

plt.savefig('difference_vector/plots/' + args.LSS + '_' + str(args.bins) + '_' + str(tension) +'.pdf', bbox_inches = 'tight')
