import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_setup(rows, columns):
	from matplotlib.gridspec import GridSpec

	plot_params = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'legend.fontsize': 8, 'font.size': 8}
	plt.rcParams.update(plot_params)

	fig = plt.figure(figsize=((3.5/3)*columns, 1.5*rows))
	gs = GridSpec(rows, columns)
	gs.update(wspace = 0, hspace=0, bottom = 0.2)

	if ((rows == 1) and (columns == 1)):
		ax = [plt.subplot(gs[0, 0])]
	else:
		ax = []
		for i in xrange(rows):
			for j in xrange(i, columns):
				ax.append(plt.subplot(gs[i, j]))
				if i != j:
					ax[-1].set_xticks([])
					ax[-1].set_yticks([])
	return ax

def plot_save(filename, show = True):
	plt.savefig(filename + '.pdf', bbox_inches='tight')
	if show:
		plt.show()
	plt.close()

def sigma(P, l):
	import numpy as np
	P_s = np.sort(P.ravel())[::-1]
	sigma = 0
	for i in P_s:
		sigma += i
		if sigma >= l:
			return i
	return 1	

def sci_lims(x, y):
	import matplotlib.pyplot as plt
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

def sci_ranges(num_params, ranges):
	sci = []
	for i in xrange(num_params-1):
		for j in xrange(i, num_params-1):
			sci.append(sci_lims(ranges[j+1], ranges[i]))
	return sci

def interpolate(bins, ranges, sci, P, i, j, a):
	import numpy as np
	from scipy.interpolate import interpn
	x = np.linspace(ranges[j+1, 0], ranges[j+1, -1], bins)
	y = np.linspace(ranges[i, 0], ranges[i, -1], bins)
	interp_domains = (x[1] - x[0]) * (y[1] - y[0])
	P_plot = np.zeros([bins, bins])
	for k in xrange(bins):
		for l in xrange(bins):
			P_plot[k, l] = 	interpn([ranges[j+1], ranges[i]], P, [x[k], y[l]], bounds_error = False, fill_value = 0)
	interp_x = x/(10**sci[a][0])
	interp_y = y/(10**sci[a][1])
	return P_plot, interp_domains, interp_x, interp_y

def crosshairs(ax, x, y):
	ax.plot([x[0], x[-1]], [0, 0], color = 'black', linestyle = ':')
	ax.plot([0, 0], [y[0], y[-1]], color = 'black', linestyle = ':')
	return ax

def plot(args, ranges, domainsize, hist, filename, plot_crosshairs = False):
	from utils import get_labels

	ax = plot_setup(args.num_params-1, args.num_params-1)

	labels = get_labels(args.chain_dir, args.CMB, args.params, args.method)
	sci = sci_ranges(args.num_params, ranges)
	a = 0
	for i in xrange(args.num_params-1):
		for j in xrange(i, args.num_params-1):
			axes = tuple([k for k in xrange(args.num_params) if ((k != i) and (k != j+1))])
			if i == j:
				ax[a].set_xlabel(labels[j+1] + ' [$\\times10^{'+ '{0:0d}'.format(int(sci[a][0])) + '}$]',labelpad=0)
				ax[a].set_ylabel(labels[i] + ' [$\\times10^{'+ '{0:0d}'.format(int(sci[a][0])) + '}$]',labelpad=0)

			domains = (ranges[j+1, 1] - ranges[j+1, 0])*(ranges[i, 1] - ranges[i, 0])
			colours = ['blue', 'red']
			linestyles = ['-', '--']
			for k in xrange(len(hist)):
				P = np.sum(hist[k]*domainsize/domains, axis = axes)
				if args.plot_interpolate != None:
					P_plot, dom, x, y = interpolate(args.plot_interpolate, ranges, sci, P, i, j, a)
				else:
					P_plot, dom, x, y = P, domains, ranges[j+1]/(10**sci[a][0]), ranges[i]/(10**sci[a][1])
				P_levels = [sigma(P_plot, 0.997/dom), sigma(P_plot, 0.95/dom), sigma(P_plot, 0.68/dom)]	
				if not all(q<r for q, r in zip(P_levels, P_levels[1:])):
					if len(hist) == 1: 
						ax[a].contourf(x, y, P_plot)
						if plot_crosshairs:
							ax[a] = crosshairs(ax[a], x, y)
					else:
						ax[a].contourf(x, y, P_plot)
				else:
					if len(hist) == 1: 
						ax[a].contour(x, y, P_plot, colors = 'purple', levels = P_levels)
						if plot_crosshairs:
							ax[a] = crosshairs(ax[a], x, y)
					else:
						ax[a].contour(x, y, P_plot, colors = colours[k], linestyles = linestyles[k], levels = P_levels)

				ax[a].locator_params(nbins=5)

			a += 1
	plot_save(filename, show = False)

def plot_surprise(D1_int, S1, D2_int, S2, filename):
	ax = plot_setup(1, 1)[0]
	ax.plot([0, 0], [0, 3], color = 'black')
	ax.barh(2, D1_int, 0.8, color = 'blue')
	ax.barh(2.1, S1, 0.6,color = 'dodgerblue')
	ax.barh(1, D2_int, 0.8, color = 'firebrick')	
	ax.barh(1.1, S2, 0.6, color = 'red')	
	xmin = min(0, np.min(D1_int))
	xmin = min(xmin, np.min(D2_int))
	xmin = min(xmin, np.min(S1))
	xmin = min(xmin, np.min(S2))
	xmax = max(0, np.max(D1_int))
	xmax = max(xmax, np.max(D2_int))
	xmax = max(xmax, np.max(S1))
	xmax = max(xmax, np.max(S2))
	lim = []
	if xmin < 0:
		lim.append(xmin)
		xmin = xmin + 0.1 * xmin
	if xmax > 0:
		lim.append(xmax)
		xmax = xmax + 0.1 * xmax
	if (xmin is 0) and (xmax is 0):
		xmin = -1
		xmax = 1
		lim.append(0)
	ax.set_xticks(lim)
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([0.8, 3])
	ax.set_yticks([])
	ax.set_xlabel('${\\rm Bits}\\times\log2$', labelpad = 0)
	plot_save(filename, show = False)
