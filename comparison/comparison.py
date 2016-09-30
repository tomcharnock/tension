import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d, interpn

def plot_setup(rows, columns):
	params = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'legend.fontsize': 8, 'font.size': 8}
	plt.rcParams.update(params)
	fig = plt.figure(figsize=((3.5/3)*columns, 1.5*rows))
	gs = GridSpec(rows, columns)
	gs.update(wspace = 0, hspace=0)
	ax = []
	for i in xrange(rows*columns):
		ax.append(plt.subplot(gs[i])) 
		ax[-1].set_xticks([])
		ax[-1].set_yticks([])
	return ax

def plot_save(filename, show = True):
	plt.savefig('comparison/plots/' + filename + '.pdf', bbox_inches='tight')
	if show:
		plt.show()
	plt.close()

def dimension(data):
	try:
		value = data.shape[0]
	except IndexError:
    		value = None
	return value

def distributions(dim, ranges, mean, std, skew = None):
	pos = np.zeros((dim,) + ranges[0].shape)
	for i in xrange(dim):
		pos[i] = ranges[i]
	pos = np.swapaxes(pos, 0, -1)
	P1 = ss.multivariate_normal(mean[0], std[0]).pdf(pos)
	P2 = ss.multivariate_normal(mean[1], std[1]).pdf(pos)
	if skew is not None:
		P1_skew = ss.multivariate_normal(skew[0][0], skew[1][0]).pdf(pos)
		P1 = (P1+P1_skew)/2.
		P2_skew = ss.multivariate_normal(skew[0][1], skew[1][1]).pdf(pos)
		P2 = (P2+P2_skew)/2.
	return P1, P2

def distribution_3(dim, ranges, mean, std, skew = None):
	pos = np.zeros((dim,) + ranges[0].shape)
	for i in xrange(dim):
		pos[i] = ranges[i]
	pos = np.swapaxes(pos, 0, -1)
	P3 = ss.multivariate_normal(mean[0], std[1]).pdf(pos)
	if skew is not None:
		P3_skew = ss.multivariate_normal(skew[0][0], skew[1][1]).pdf(pos)
		P3 = (P3+P3_skew)/2.
	return P3

def sigma(P,l):
	P_s = np.sort(P.ravel())[::-1]
	sigma = 0
	for i in P_s:
		sigma += i
		if sigma >= l:
			return i
	return 1	

def Bhattacharrya(ax, gridsize, ranges, mean = None, std = None, skew = None, P12 = None, dim = 1):
	def one_dimension(ax, ranges, P1, P2, P, BC):
		ax.plot(ranges, P1, color='blue')
		ax.plot(ranges, P2, color='red', linestyle = '--')
	        ax.plot(ranges, P, color = 'black', linestyle = ':')
        	ax.fill_between(ranges, 0, P, color = 'black', alpha = 0.5)
	        ax.text(0.05, 0.95, '$B = ' + '{0:.3f}'.format(BC) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	
	
	def two_dimensions(ax, ranges, P1, P2, P, BC, domainsize):
		P1_levels = [sigma(P1, 0.997/domainsize), sigma(P1, 0.95/domainsize), sigma(P1, 0.68/domainsize)]
		P2_levels = [sigma(P2, 0.997/domainsize), sigma(P2, 0.95/domainsize), sigma(P2, 0.68/domainsize)]
	        ax.contourf(ranges[0], ranges[1], P, alpha = 0.5, levels = np.linspace(0, max(np.max(P1), np.max(P2)), 100))
		ax.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
		ax.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
	        ax.text(0.05, 0.95, '$B = ' + '{0:.3f}'.format(BC) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	

	if P12 is None:
		P1, P2 = distributions(dim, ranges, mean, std, skew = skew)
	else:
		P1, P2 = P12

	P = np.sqrt(P1*P2)
	domainsize = 1
	for i in xrange(dim):
		domainsize = domainsize * gridsize
		
        BC = np.sum(P*domainsize)
	if dim is 1:
		one_dimension(ax, ranges[0], P1, P2, P, BC)
	else:
		two_dimensions(ax, ranges, P1, P2, P, BC, domainsize)


def overlap(ax, gridsize, ranges, mean = None, std = None, skew = None, P12 = None, dim = 1):
	def one_dimension(ax, ranges, P1, P2, M, M_val):
		ax.plot(ranges, P1, color='blue')
		ax.plot(ranges, P2, color='red', linestyle = '--')
		ax.fill_between(ranges, 0, M, color = 'black', alpha = 0.5)
		ax.text(0.05, 0.95, '$O = ' + '{0:.3f}'.format(M_val) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

	def two_dimensions(ax, ranges, P1, P2, M, M_val, domainsize):
		P1_levels = [sigma(P1, 0.997/domainsize), sigma(P1, 0.95/domainsize), sigma(P1, 0.68/domainsize)]
		P2_levels = [sigma(P2, 0.997/domainsize), sigma(P2, 0.95/domainsize), sigma(P2, 0.68/domainsize)]
	        ax.contourf(ranges[0], ranges[1], M, alpha = 0.5, levels = np.linspace(0, max(np.max(P1), np.max(P2)), 100))
		ax.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
		ax.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
	        ax.text(0.05, 0.95, '$O = ' + '{0:.3f}'.format(M_val) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	

	if P12 is None:
		P1, P2 = distributions(dim, ranges, mean, std, skew = skew)
	else:
		P1, P2 = P12

	M = np.minimum(P1,P2)
	domainsize = 1
	for i in xrange(dim):
		domainsize = domainsize * gridsize

	M_val = np.sum(M*domainsize)
	if dim is 1:
		one_dimension(ax, ranges[0], P1, P2, M, M_val)
	else:
		two_dimensions(ax, ranges, P1, P2, M, M_val, domainsize)

def under(ax, ranges, mean = None, std = None, skew = None, P12 = None, level = 0.997300, dim = 1):
	def one_dimension(ax, ranges, P1, P2, under_P1, under_P2, domainsize):
		ax1 = ax[0]
		ax2 = ax[1]
		ax1.plot(ranges, P1, color = 'blue')
		ax1.plot(ranges, P2, color = 'red', linestyle = '--')
		ax2.plot(ranges, P1, color = 'blue')
		ax2.plot(ranges, P2, color = 'red', linestyle = '--')
		ax1.fill_between(ranges, 0, P1, where = P2 >= sigma(P2, level/domainsize), color = 'black', alpha = 0.5)
		ax1.text(0.05, 0.95, '$I_1 = ' + '{0:.3f}'.format(np.sum(under_P1 * domainsize)) + '$', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes,color='black')
		ax2.fill_between(ranges, 0, P2, where = P1 >= sigma(P1, level/domainsize), color = 'black', alpha = 0.5)
		ax2.text(0.05, 0.95, '$I_2 = ' + '{0:.3f}'.format(np.sum(under_P2 * domainsize)) + '$', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes,color='black')

	def two_dimensions(ax, ranges, P1, P2, under_P1, under_P2, domainsize):
		ax1 = ax[0]
		ax2 = ax[1]
		P1_levels = [sigma(P1, 0.997/domainsize), sigma(P1, 0.95/domainsize), sigma(P1, 0.68/domainsize)]
		P2_levels = [sigma(P2, 0.997/domainsize), sigma(P2, 0.95/domainsize), sigma(P2, 0.68/domainsize)]
	        ax1.contourf(ranges[0], ranges[1], under_P1, alpha = 0.5, levels = np.linspace(np.min(under_P1[under_P1>0]), np.max(under_P1), 100))
	        ax2.contourf(ranges[0], ranges[1], under_P2, alpha = 0.5, levels = np.linspace(np.min(under_P2[under_P2>0]), np.max(under_P2), 100))
		ax1.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
		ax2.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
		ax1.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
		ax2.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
	        ax1.text(0.05, 0.95, '$I_1 = ' + '{0:.3f}'.format(np.sum(under_P1 * domainsize)) + '$', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)	
	        ax2.text(0.05, 0.95, '$I_2 = ' + '{0:.3f}'.format(np.sum(under_P2 * domainsize)) + '$', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)	

	if P12 is None:
		P1, P2 = distributions(dim, ranges, mean, std, skew = skew)
	else:
		P1, P2 = P12

	domainsize = 1
	for i in xrange(dim):
		domainsize = domainsize * gridsize
	P1_norm = P1 * domainsize
	P2_norm = P2 * domainsize
	under_P1 = np.copy(P1)
	under_P1[P2_norm < sigma(P2_norm, level)] = 0
	under_P2 = np.copy(P2)
	under_P2[P1_norm < sigma(P1_norm, level)] = 0
	if dim is 1:
		one_dimension(ax, ranges[0], P1, P2, under_P1, under_P2, domainsize)
	else:
		two_dimensions(ax, ranges, P1, P2, under_P1, under_P2, domainsize)

def Charnock(ax, mean, std, skew = None, num = 1e8, bins = 10, dim = 1, interp_bins = 100):
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

	def difference_vector(dim, num, mean, std, skew = None):
		P1 = np.random.multivariate_normal(mean[0], std[0], num)
		P2 = np.random.multivariate_normal(mean[1], std[1], num)
		if skew is not None:
			P1_skew = np.random.multivariate_normal(skew[0][0], skew[1][0], num/2)
			P2_skew = np.random.multivariate_normal(skew[0][1], skew[1][1], num/2)
			P1_out = np.zeros([num, dim])
			P2_out = np.zeros([num, dim])
			for i in xrange(dim):
				P1_out[:,i] = np.append(P1[:num/2,i],P1_skew[:,i])
				P2_out[:,i] = np.append(P2[:num/2,i],P2_skew[:,i])
			P1 = P1_out
			P2 = P2_out
		P = P1 - P2
		return P

	def one_dimension(ax, ranges, interp_ranges, P_hist, interp_domainsize):
		P_spline = extrapolate(interp1d(ranges, P_hist))
		P_plot = P_spline(interp_ranges)
		P_0 = P_spline([0])
		C = np.sum(P_plot[P_plot >= P_0] * interp_domainsize)
		ax.plot(interp_ranges, P_plot, color='purple')
		ax.fill_between(interp_ranges, 0, P_plot, where = P_plot >= P_0, color = 'black', alpha = 0.5)
		ax.text(0.05, 0.95, '$C = ' + '{0:.3f}'.format(C) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
		ax.set_xticks([0])	
	
	def two_dimension(ax, ranges, interp_bins, interp_ranges, P_hist, interp_domainsize, mean):
		P_plot = np.zeros([interp_bins, interp_bins])
		for i in xrange(interp_bins):
			for j in xrange(interp_bins):
				grid = np.array([interp_ranges[0, i], interp_ranges[1, j]])
				P_plot[i, j] = interpn(tuple([ranges[k] for k in xrange(dimension(mean))]), P_hist, grid, bounds_error = False, fill_value = 0)
		P_0 = interpn(ranges, P_hist, [0, 0], bounds_error = False, fill_value = 0)
		P_C = np.copy(P_plot)
		P_C[P_C <= P_0] = 0
		C = np.sum(P_C * interp_domainsize)
		P_levels = [sigma(P_plot, 0.997/interp_domainsize), sigma(P_plot, 0.95/interp_domainsize), sigma(P_plot, 0.68/interp_domainsize)]
		if P_0 != 0:
			ax.contourf(interp_ranges[0], interp_ranges[1], P_C, alpha = 0.5, levels = np.linspace(np.min(P_plot[P_plot > 0]), np.max(P_plot), 100))
		else:
			ax.contourf(interp_ranges[0], interp_ranges[1], P_C, alpha = 0.5, levels = np.linspace(0, np.max(P_plot), 100))
		ax.contour(interp_ranges[0], interp_ranges[1], P_plot, colors = 'purple', levels = P_levels)
		ax.plot([np.min(interp_ranges[0]), np.max(interp_ranges[0])], [0, 0], linestyle = ':', color = 'black')
		ax.plot([0, 0], [np.min(interp_ranges[1]), np.max(interp_ranges[1])], linestyle = ':', color = 'black')
		ax.text(0.05, 0.95, '$C = ' + '{0:.3f}'.format(C) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

	P = difference_vector(dim, num, mean, std, skew = skew)
	P_hist, edges = np.histogramdd(P, bins = bins)
	ranges = np.zeros([dim, bins])
	domainsize = 1
	interp_ranges = np.zeros([dim, interp_bins])
	interp_domainsize = 1
	for i in xrange(dim):
		ranges[i] = np.array([edges[i][j]+(edges[i][j+1]-edges[i][j])/2 for j in xrange(int(bins))])
		domainsize = domainsize * (ranges[i][1] - ranges[i][0])
		interp_ranges[i] = np.linspace(min(-0.1, np.min(ranges[i])), max(-0.1, np.max(ranges[i])), interp_bins)
		interp_domainsize = interp_domainsize * (interp_ranges[i][1] - interp_ranges[i][0])
	P_hist = P_hist/np.sum(P_hist*domainsize)

	if dim is 1:
		one_dimension(ax, ranges[0], interp_ranges[0], P_hist, interp_domainsize)
	else:
		two_dimension(ax, ranges, interp_bins, interp_ranges, P_hist, interp_domainsize, mean)

def surprise(ax, gridsize, ranges, mean = None, std = None, skew = None, P12 = None, dim = 1):
	def one_dimension(ax, ranges, P1, P2, D1_int, D2_int, D1_av, D2_av, S1, S2):
		if len(ax) == 2:
			ax1 = ax[0]
			ax2 = ax[1]
			ax1.plot(ranges, P1, color = 'blue')
			ax1.plot(ranges, P2, color = 'red', linestyle = '--')
		else:
			ax2 = ax[0]
		ax2.plot([0, 0], [0, 3], color = 'black')
		ax2.barh(2, D1_int, 0.8, color = 'blue')
		ax2.barh(2.1, S1, 0.6,color = 'dodgerblue')
		ax2.barh(1, D2_int, 0.8, color = 'firebrick')	
		ax2.barh(1.1, S2, 0.6, color = 'red')	
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
		ax2.set_xticks(lim)
		ax2.set_xlim([xmin, xmax])
		ax2.set_ylim([0.8, 3])
	
	def two_dimensions(ax, ranges, P1, P2, D1_int, D2_int, D1_av, D2_av, S1, S2, domainsize):
		if len(ax) == 2:
			ax1 = ax[0]
			ax2 = ax[1]
			P1_levels = [sigma(P1, 0.997/domainsize), sigma(P1, 0.95/domainsize), sigma(P1, 0.68/domainsize)]
			P2_levels = [sigma(P2, 0.997/domainsize), sigma(P2, 0.95/domainsize), sigma(P2, 0.68/domainsize)]
			ax1.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
			ax1.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
		else:
			ax2 = ax[0]
		ax2.plot([0, 0], [0, 3], color = 'black')
		ax2.barh(2, D1_int, 0.8, color = 'blue')
		ax2.barh(2.1, S1, 0.6,color = 'dodgerblue')
		ax2.barh(1, D2_int, 0.8, color = 'firebrick')	
		ax2.barh(1.1, S2, 0.6, color = 'red')	
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
			lim.append(xmin/3)
			xmin = xmin + 0.1 * xmin
		if xmax > 0:
			lim.append(xmax)
			xmax = xmax + 0.1 * xmax
		if (xmin is 0) and (xmax is 0):
			xmin = -1
			xmax = 1
			lim.append(0)
		ax2.set_xticks(lim)
		ax2.set_xlim([xmin, xmax])
		ax2.set_ylim([0.8, 3])

	if P12 is None:
		P1, P2 = distributions(dim, ranges, mean, std, skew = skew)
	else:
		P1, P2 = P12

	domainsize = 1
	for i in xrange(dim):
		domainsize = domainsize * gridsize

	D1 = P2*np.log(P2/P1)
	D1_int = np.sum(D1*domainsize)
	D2 = P1*np.log(P1/P2)
	D2_int = np.sum(D2*domainsize)
	P = P1*P2
	D1P = D1*P
	D2P = D2*P
	D1_av = np.sum(D1P)
	D2_av = np.sum(D2P)
	S1 = D1_int - np.abs(D1_av)
	S2 = D2_int - np.abs(D2_av)
	print 'D1 = ', D1_int, ' and D2 = ', D2_int
	print '<D1> = ', D1_av, ' and <D2> = ', D2_av
	print 'S1 = ', S1, ' and S2 = ', S2	
	if dim is 1:
		one_dimension(ax, ranges[0], P1, P2, D1_int, D2_int, D1_av, D2_av, S1, S2)
	else:
		two_dimensions(ax, ranges, P1, P2, D1_int, D2_int, D1_av, D2_av, S1, S2, domainsize)

def marshall(ax, gridsize, ranges, mean = None, std = None, skew = None, P12 = None, dim = 1):
	def one_dimension(ax, ranges, P1, P2, logR):
		ax.plot(ranges, P1, color = 'blue')
		ax.plot(ranges, P2, color = 'red', linestyle = '--')
	        ax.text(0.05, 0.95, '$\log R = ' + '{0:.3f}'.format(logR) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	
	
	def two_dimensions(ax, ranges, P1, P2, logR, domainsize):
		P1_levels = [sigma(P1, 0.997/domainsize), sigma(P1, 0.95/domainsize), sigma(P1, 0.68/domainsize)]
		P2_levels = [sigma(P2, 0.997/domainsize), sigma(P2, 0.95/domainsize), sigma(P2, 0.68/domainsize)]
		ax.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
		ax.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
	        ax.text(0.05, 0.95, '$\log R = ' + '{0:.3f}'.format(logR) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	

	if P12 is None:
		P1, P2 = distributions(dim, ranges, mean, std, skew = skew)
	else:
		P1, P2 = P12

	domainsize = 1
	for i in xrange(dim):
		domainsize = domainsize * gridsize
	
	prior = np.ones(P1.shape)
	prior = prior/np.sum(prior*domainsize)
	P = np.sum(P1*P2*domainsize)
	P2_int = np.sum(P2*prior*domainsize)
	logR = np.log(P/P2_int)

	if dim is 1:
		one_dimension(ax, ranges[0], P1, P2, logR)
	else:
		two_dimensions(ax, ranges, P1, P2, logR, domainsize)

def verde(ax, gridsize, ranges, P3, mean = None, std = None, skew = None, P12 = None, dim = 1):
	def one_dimension(ax, ranges, P1, P2, logT):
		ax.plot(ranges, P1, color = 'blue')
		ax.plot(ranges, P2, color = 'red', linestyle = '--')
	        ax.text(0.05, 0.95, '$\log R = ' + '{0:.3f}'.format(logT) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	
	
	def two_dimensions(ax, ranges, P1, P2, logT, domainsize):
		P1_levels = [sigma(P1, 0.997/domainsize), sigma(P1, 0.95/domainsize), sigma(P1, 0.68/domainsize)]
		P2_levels = [sigma(P2, 0.997/domainsize), sigma(P2, 0.95/domainsize), sigma(P2, 0.68/domainsize)]
		ax.contour(ranges[0], ranges[1], P1, colors = 'blue', levels = P1_levels)
		ax.contour(ranges[0], ranges[1], P2, colors = 'red', linestyles = 'dashed', levels = P2_levels)
	        ax.text(0.05, 0.95, '$\log R = ' + '{0:.3f}'.format(logT) + '$', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)	

	if P12 is None:
		P1, P2 = distributions(dim, ranges, mean, std, skew = skew)
	else:
		P1, P2 = P12

	domainsize = 1
	for i in xrange(dim):
		domainsize = domainsize * gridsize
	
	P_shift = np.sum(P1*P3*gridsize)
	P = np.sum(P1*P2*gridsize)
	logT = np.log(P_shift/P)

	if dim is 1:
		one_dimension(ax, ranges[0], P1, P2, logT)
	else:
		two_dimensions(ax, ranges, P1, P2, logT, domainsize)

gridsize = 0.02
x, y = np.mgrid[-10:10:gridsize, -10:10:gridsize]

#Paper Plots
mean = np.array([[[0], [0]], [[0], [0]], [[5], [-5]], [[0], [1.427]], [[0], [1.427]], [[0, 0], [0, 0]], [[0, 0], [0, 0]], [[5, 5], [-5, -5]], [[0, 0], [1.427, 1.427]], [[0, 0], [1.427, 1.427]]])
std = np.array([[[[1]], [[1]]], [[[1]], [[3**2]]], [[[1]], [[1]]], [[[1]], [[1]]], [[[1]], [[1]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]], [[[1, 0], [0, 1]], [[3**2, 0], [0, 3**2]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]]])
skew1D = [np.array([[-2], [4]]), np.array([[[1]], [[2**2]]])]
skew2D = [np.array([[-2, -2], [4, 4]]), np.array([[[1, 0], [0, 1]], [[2**2, 0], [0, 2**2]]])]
'''
#Bhattacharrya
ax = plot_setup(2, 5)
for i in xrange(10):
	if i < 5:
		if i != 4:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i])
		else:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i], skew = skew1D)
		Bhattacharrya(ax[i], gridsize, [x[:, 0]], P12 = np.array([P1, P2]), dim = 1)
	else:
		if i != 9:
			P1, P2 = distributions(2, [x, y], mean[i], std[i])
		else:
			P1, P2 = distributions(2, [x, y], mean[i], std[i], skew = skew2D)
		Bhattacharrya(ax[i], gridsize, [x, y], P12 = np.array([P1, P2]), dim = 2)
plot_save('Bhattacharrya')
'''

'''
#Overlap
ax = plot_setup(2, 5)
for i in xrange(10):
	if i < 5:
		if i != 4:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i])
		else:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i], skew = skew1D)
		overlap(ax[i], gridsize, [x[:, 0]], P12 = np.array([P1, P2]), dim = 1) 
	else:
		if i != 9:
			P1, P2 = distributions(2, [x, y], mean[i], std[i])
		else:
			P1, P2 = distributions(2, [x, y], mean[i], std[i], skew = skew2D)
		overlap(ax[i], gridsize, [x, y], P12 = np.array([P1, P2]), dim = 2) 
plot_save('Overlap')
'''

'''
#Under
ax = plot_setup(4, 5)
for i in xrange(10):
	if i < 5:
		if i != 4:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i])
		else:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i], skew = skew1D)
		under([ax[i], ax[i+5]], [x[:, 0]], P12 = np.array([P1, P2]), dim = 1)
	else:
		if i != 9:
			P1, P2 = distributions(2, [x, y], mean[i], std[i])
		else:
			P1, P2 = distributions(2, [x, y], mean[i], std[i], skew = skew2D)
		under([ax[5+i], ax[10+i]], [x, y], P12 = np.array([P1, P2]), dim = 2)
plot_save('Under')
'''


#Surprise
surprise_mean_1D = mean[:5] # np.array([[[0], [0]], [[0], [0]], [[5], [-5]], [[0], [1.427]], [[0], [1.75]], [[0], [1.427]]])
surprise_std_1D = std[:5] #np.array([[[[1]], [[1]]], [[[1]], [[3**2]]], [[[1]], [[1]]], [[[1]], [[1]]], [[[1]], [[1]]], [[[1]], [[1]]]])
#ax = plot_setup(2, 6)
#for i in xrange(6):
#	if i != 5:
ax = plot_setup(1, 5)
for i in xrange(5):
	if i != 4:
		P1, P2 = distributions(1, [x[:, 0]], surprise_mean_1D[i], surprise_std_1D[i])
	else:
		P1, P2 = distributions(1, [x[:, 0]], surprise_mean_1D[i], surprise_std_1D[i], skew = skew1D)
	#surprise([ax[i], ax[i+6]], gridsize, [x[:, 0]], P12 = np.array([P1, P2]), dim = 1)
	surprise([ax[i]], gridsize, [x[:, 0]], P12 = np.array([P1, P2]), dim = 1)

plot_save('Surprise1D')

surprise_mean_2D = np.array([[[0, 0], [0, 0]], [[0, 0], [0, 0]], [[5, 5], [-5, -5]], [[0, 0], [1.8, 1.8]], [[0, 0], [1.427, 1.427]]])
surprise_std_2D = np.array([[[[1, 0], [0, 1]], [[1, 0], [0, 1]]], [[[1, 0], [0, 1]], [[3**2, 0], [0, 3**2]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]], [[[1, 0], [0, 1]], [[1, 0], [0, 1]]]])
#ax = plot_setup(2, 6)
#for i in xrange(6):
#	if i != 5:
ax = plot_setup(1, 5)
for i in xrange(5):
	if i != 4:
		P1, P2 = distributions(2, [x, y], surprise_mean_2D[i], surprise_std_2D[i])
	else:
		P1, P2 = distributions(2, [x, y], surprise_mean_2D[i], surprise_std_2D[i], skew = skew2D)
	#surprise([ax[i], ax[i+6]], gridsize, [x, y], P12 = np.array([P1, P2]), dim = 2)
	surprise([ax[i]], gridsize, [x, y], P12 = np.array([P1, P2]), dim = 2)

plot_save('Surprise2D')

'''
#Marshall
ax = plot_setup(2, 5)
for i in xrange(10):
	if i < 5:
		if i != 4:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i])
		else:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i], skew = skew1D)
		marshall(ax[i], gridsize, [x[:, 0]], P12 = np.array([P1, P2]), dim = 1)
	else:
		if i != 9:
			P1, P2 = distributions(2, [x, y], mean[i], std[i])
		else:
			P1, P2 = distributions(2, [x, y], mean[i], std[i], skew = skew2D)
		marshall(ax[i], gridsize, [x, y], P12 = np.array([P1, P2]), dim = 2)
plot_save('Marshall')
'''

'''
#Verde
ax = plot_setup(2, 5)
for i in xrange(10):
	if i < 5:
		if i != 4:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i])
		else:
			P1, P2 = distributions(1, [x[:, 0]], mean[i], std[i], skew = skew1D)
		P3 = distribution_3(1, [x[:, 0]], mean[i], std[i])
		verde(ax[i], gridsize, [x[:, 0]], P3, P12 = np.array([P1, P2]), dim = 1)
	else:
		if i != 9:
			P1, P2 = distributions(2, [x, y], mean[i], std[i])
		else:
			P1, P2 = distributions(2, [x, y], mean[i], std[i], skew = skew2D)
		P3 = distribution_3(2, [x, y], mean[i], std[i])
		verde(ax[i], gridsize, [x, y], P3, P12 = np.array([P1, P2]), dim = 2)
plot_save('Verde')
'''

'''
#Charnock
ax = plot_setup(2, 5)
for i in xrange(10):
	if i < 5:
		if i != 4:
			Charnock(ax[i], mean[i], std[i], num = 10000000, dim = 1 , bins = 100)
		else:
			Charnock(ax[i], mean[i], std[i], skew = skew1D, num = 10000000, dim = 1 , bins = 100)
	else:
		if i != 9:
			Charnock(ax[i], mean[i], std[i], num = 10000000, dim = 2, bins = 20)
		else:
			Charnock(ax[i], mean[i], std[i], skew = skew2D, num = 10000000, dim = 2, bins = 20)
plot_save('Charnock')
'''
