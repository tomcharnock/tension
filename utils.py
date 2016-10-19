import numpy as np

class parameters():
	def __init__(self, arguments):
		import sys
		self.method = arguments['method']
		self.chain_dir = arguments['chain_dir']
		self.CMB = arguments['CMB']
		self.CMB_chains = arguments['CMB_chains']
		self.LSS = arguments['LSS']
		self.LSS_chains = arguments['LSS_chains']
		self.params = arguments['params']
		self.num_params = len(self.params)
		self.bins = arguments['bins']
		self.plot_dir = arguments['plot_dir']
		self.save_dir = arguments['save_dir']
		self.load = arguments['load']
		self.priors = arguments['priors']
		self.smoothing = arguments['smoothing']
		self.interpolate_bins = arguments['interpolate_bins']
		self.integration_bounds = arguments['integration_bounds']
		self.plot_interpolate = arguments['plot_interpolate'] 
                self.sampling_method = arguments['sampling_method']
		self.sampling_parameter = arguments['sampling_parameter']
                self.sampling_constraints = arguments['sampling_constraints']
		if not 'difference_vector' in self.method:
			if self.interpolate_bins != None:
				if self.interpolate_bins < 1:
					print 'Number of interpolate bins too low'
					sys.exit()
			if (('marshall' in self.method) or ('verde' in self.method)):
				if self.load != None:
					if self.priors == None:
						print 'Need to set the directory to load priors from'
						sys.exit()
				if 'verde' in self.method:
					if self.load != None:
						if len(self.load) < 3:
							print 'Need to set the directory to load the shifted distribution from'
			if 'ibi' in self.method:
				if self.integration_bounds == None:
					print 'Need to set the integration bounds'
					sys.exit()
				if self.integration_bounds > 1:
					print 'Integration bounds are too high'
					sys.exit()
				if self.integration_bounds < 0:
					print 'Integration bounds are too low'
					sys.exit()

def apply_smoothing(hist, amount):
	from scipy.ndimage.filters import gaussian_filter
	smoothing = [amount for i in xrange(len(hist.shape))]
	hist = gaussian_filter(hist, smoothing)
	return hist

def combine(args):
	if args.load != None:
		CMB_hist, ranges, domainsize = load(args.load[0])
		LSS_hist, ranges, domainsize = load(args.load[1])
		if 'marshall' in args.method:
			shifted_hist = None
			priors, ranges, domainsize = load(args.priors)
		elif 'verde' in args.method:
			shifted_hist, ranges, domainsize = load(args.load[2])
			priors, ranges, domainsize = load(args.priors)
		else:
			priors = None
			shifted_hist = None
	else:
		from itertools import product
		from scipy.interpolate import interpn
		CMB_samples = get_samples(args, args.CMB, args.CMB_chains)
		CMB_hist, CMB_ranges, CMB_domainsize = get_histogram(args, CMB_samples)
		LSS_samples = get_samples(args, args.LSS, args.LSS_chains)
		LSS_hist, LSS_ranges, LSS_domainsize = get_histogram(args, LSS_samples)
		ranges, bins = interpolation_bins(CMB_ranges, LSS_ranges, args.num_params, num_bins = args.interpolate_bins)
		domainsize = get_domainsize(ranges)
		interp_ranges = list(product(*ranges))
		CMB_hist_interp = interpn(CMB_ranges, CMB_hist, interp_ranges, bounds_error = False, fill_value = 0)
		CMB_hist_interp[CMB_hist_interp < 0] = 0
		CMB_hist_interp = np.reshape(CMB_hist_interp, tuple([bins[i] for i in xrange(args.num_params)]))
		CMB_hist_interp = CMB_hist_interp/np.sum(CMB_hist_interp)/domainsize
		if args.save_dir != None:
			CMB_filename = get_filename(args.save_dir, args.CMB, 'interpolated')
			save(CMB_filename, CMB_hist_interp, ranges, domainsize)
			del CMB_hist, CMB_ranges, CMB_domainsize
		LSS_hist_interp = interpn(LSS_ranges, LSS_hist, interp_ranges, bounds_error = False, fill_value = 0)
		LSS_hist_interp[LSS_hist_interp < 0] = 0
		LSS_hist_interp = np.reshape(LSS_hist_interp, tuple([bins[i] for i in xrange(args.num_params)]))
		LSS_hist_interp = LSS_hist_interp/np.sum(LSS_hist_interp)/domainsize
		if args.save_dir != None:
			LSS_filename = get_filename(args.save_dir, args.LSS, 'interpolated')
			save(LSS_filename, LSS_hist_interp, ranges, domainsize)
			del LSS_hist, LSS_ranges, LSS_domainsize
		if 'verde' in args.method:
			shifted_hist, shifted_ranges = get_shifted(args, CMB_samples, LSS_samples)
			shifted_hist_interp = interpn(shifted_ranges, shifted_hist, interp_ranges, bounds_error = False, fill_value = 0)
			shifted_hist_interp[shifted_hist_interp < 0] = 0
			shifted_hist_interp = np.reshape(shifted_hist_interp, tuple([bins[i] for i in xrange(args.num_params)]))
			shifted_hist_interp = shifted_hist_interp/np.sum(shifted_hist_interp)/domainsize
			if args.save_dir != None:
				shifted_filename = get_filename(args.save_dir, args.LSS, 'shifted_interpolated')
				save(shifted_filename, shifted_hist_interp, ranges, domainsize)
				del shifted_hist, shifted_ranges
			shifted_hist = shifted_hist_interp
		else:
			shifted_hist = None
		if (('marshall' in args.method) or ('verde' in args.method)):
			if args.priors != None:
				priors, priors_ranges, priors_domainsize = load(args.priors)
			else:
				priors = get_priors(args, CMB_hist_interp, domainsize)
				if args.save_dir != None:
					priors_filename = get_filename(args.save_dir, args.LSS, 'priors_interpolated')
					save(priors_filename, priors, ranges, domainsize)
		else:
			priors = None
		CMB_hist = CMB_hist_interp
		LSS_hist = LSS_hist_interp
	return CMB_hist, LSS_hist, shifted_hist, priors, ranges, domainsize

def difference_vector(args):
	if args.load != None:
		hist, ranges, domainsize = load(args.load[0])
	else:
		import multiprocessing
		from multiprocessing import Pool
		from functools import partial
		CMB = get_samples(args, args.CMB, args.CMB_chains)
		LSS = get_samples(args, args.LSS, args.LSS_chains)
		ranges, bins, domainsize = get_ranges(args.num_params, args.bins, CMB = CMB, LSS = LSS)
		hist = np.zeros([args.bins, args.bins, args.bins, args.bins, args.bins])
		print "Using ", multiprocessing.cpu_count(), " processors"
		pool = Pool()
		histogram_pass = partial(difference_histogram, CMB, args.bins, bins)
		for i in pool.imap_unordered(histogram_pass, LSS):
			hist += i
		if args.smoothing != 0:
			hist = apply_smoothing(hist, args.smoothing)
		hist = hist/(np.sum(hist)*domainsize)
		if args.save_dir != None:
			filename = get_filename(args.save_dir, args.LSS, 'difference_vector')
			save(filename, hist, ranges, domainsize)
	return hist, ranges, domainsize

def difference_histogram(CMB, num_bins, bins, LSS):
	diff = CMB-LSS
	res, edges = np.histogramdd(diff, bins = num_bins, range = bins)
	return res

def define_hist(CMB, LSS, num_params, num_bins):
	bins = [[np.min(CMB[:, i])-np.max(LSS[:, i]), np.max(CMB[:, i])-np.min(LSS[:, i])] for i in xrange(num_params)]
	a, edges = np.histogramdd([[0] for i in xrange(num_params)], bins = num_bins, range = bins)
	return bins, edges

def get_columns(args, chain, params):
	columns = [0 for i in xrange(len(params))]
	filename = args.chain_dir + chain + ".paramnames"
	with open(filename) as f:
		i = 0
		for read_line in f:
			line = filter(None, read_line.replace('\n','').split('\t'))
			if line[0] in params:
				columns[params.index(line[0])] = int(i+2)
			i += 1
	return columns

def get_D(P1, P2, domainsize):
	D = P2*np.log(P2/P1)
	D[np.isnan(D)] = 0.
	D_int = np.sum(D*domainsize)
	return D, D_int

def get_domainsize(ranges):
	num_params = len(ranges)
	domainsize = 1
	for i in xrange(num_params):
		domainsize = domainsize * (ranges[i][1] - ranges[i][0])
	return domainsize

def get_filename(output, data, method, bins = None, extra = None):
	filename = output + data + '_' + method
	if bins != None:
		filename = filename +'_' + str(bins)
	if extra != None:
		filename = filename + '_' + str(extra)
	return filename

def get_histogram(args, samples):
	hist, edges = np.histogramdd(samples, bins = args.bins)
	if args.smoothing != 0:
		hist = apply_smoothing(hist, args.smoothing)
	ranges, domainsize = get_ranges(args.num_params, args.bins, edges = edges)
	hist = hist/np.sum(hist)/domainsize
	return hist, ranges, domainsize

def get_ibi(P1, P2, domainsize, integration_bounds):
	from plot import sigma
	under_P2 = np.copy(P1)
	under_P2[P2 < sigma(P2, integration_bounds/domainsize)] = 0
	I = np.sum(under_P2*domainsize)
	return I

def get_labels(args, chain):
	if args.method == 'difference_vector':
		label_dict = {'omegabh2': '$\Delta\Omega_{\\rm b}h^2$', 'omegach2': '$\Delta\Omega_{\\rm c}h^2$', 'theta': '$\Delta\Theta_{\\rm MC}$', 'logA': '$\Delta\log A_{\\rm s}$', 'ns': '$\Delta n_{\\rm s}$', 'mnu': '$\Delta\sum m_\\nu$', 'meffsterile': '$\Delta m_{\\rm eff}^{\\rm sterile}$', 'nnu': '$\Delta N_{\\rm eff}$'}
	else:
		label_dict = {'omegabh2': '$\Omega_{\\rm b}h^2$', 'omegach2': '$\Omega_{\\rm c}h^2$', 'theta': '$\Theta_{\\rm MC}$', 'logA': '$\log A_{\\rm s}$', 'ns': '$n_{\\rm s}$', 'mnu': '$\sum m_\\nu$', 'meffsterile': '$m_{\\rm eff}^{\\rm sterile}$', 'nnu': '$N_{\\rm eff}$'}
	labels = []
	filename = args.chain_dir + chain + ".paramnames"
	with open(filename) as f:
		i = 0
		for read_line in f:
			line = filter(None, read_line.replace('\n','').split('\t'))
			if line[0] in args.params:
				labels.append(label_dict[line[0]])
			i += 1
	return labels

def get_likelihood(args, dataset, dataset_num, column = 1):
	data = []
	for num in xrange(int(dataset_num)):
		filename = args.chain_dir + dataset + "_" + str(num+1) + ".txt"
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
						data.append(float(line[column]))
				num_lines += 1
	likelihood = np.array(data)
	return likelihood

def get_priors(args, hist, domainsize):
	priors = np.ones(hist.shape)
	priors = priors/np.sum(priors*domainsize)
	return priors
	
def get_ranges(num_params, num_bins, edges = None, CMB = None, LSS = None):
	bins = None
	if edges == None:
		bins, edges = define_hist(CMB, LSS, num_params, num_bins)
	ranges = np.array([[edges[j][i]+(edges[j][i+1]-edges[j][i])/2 for i in xrange(num_bins)] for j in xrange(num_params)])
	domainsize = get_domainsize(ranges)
	if bins == None:
		return ranges, domainsize
	else:
		return ranges, bins, domainsize

def get_S(P1, P2, D_int):
	D_av = np.mean(np.log(P1)-np.log(P2))
	S = D_int - np.abs(D_av)
	return S

def get_samples(args, dataset, dataset_num):
	samples = read_chains(args, dataset, dataset_num)
	print "Number of samples in " + dataset + " = ",len(samples)
	print "Shape of data = ", samples.shape
	return samples

def get_shifted(args, CMB_samples, LSS_samples):
	CMB_means = np.mean(CMB_samples, axis = 0)
	LSS_means = np.mean(LSS_samples, axis = 0)
	shifted_samples = LSS_samples - LSS_means + CMB_means
	shifted_hist, shifted_ranges, shifted_domainsize = get_histogram(args, shifted_samples)
	return shifted_hist, shifted_ranges

def get_surprise(P1, P2, domainsize):
	idx1 = np.where(P1 != 0)
	P1 = P1[idx1]
	P2 = P2[idx1]
	idx2 = np.where(P2 != 0)
	P1 = P1[idx2]
	P2 = P2[idx2]
	D, D_int = get_D(P1, P2, domainsize)
	S = get_S(P1, P2, D_int)
	return D_int, S

def get_tension(hist, ranges, domainsize, num_params):
	from scipy.interpolate import interpn
	from scipy.special import erfinv
	val = interpn(ranges, hist, [0 for i in xrange(num_params)], bounds_error = False, fill_value = 0)
	C = np.sum(hist[hist >= val]*domainsize)
	tension = np.sqrt(2)*erfinv(C)
	return C,tension

def importance_sampling(args, samples, dataset, dataset_num):
	likelihood = get_likelihood(args, dataset, dataset_num)
	new_samples = np.copy(samples)
	index = np.where(new_samples[:, 0] != np.nan)[0]
	for i in xrange(len(args.sampling_method)):
		sample = True
		if args.sampling_parameter[i] not in args.params:
			column = get_columns(args, dataset, [args.sampling_parameter[i]])[0]
			if column < 2:
				sample = False
			parameter_samples = get_likelihood(args, dataset, dataset_num, column = column)[index]
		else:
			column = np.where(np.array(args.params) == args.sampling_parameter[i])[0][0]
			parameter_samples = new_samples[index, column]
		if sample:
			if args.sampling_method[i] == 'uniform':
				prior = np.ones(len(parameter_samples))
				index = np.where(parameter_samples < args.sampling_constraints[i][0])[0]
				prior[index] = 0.
				index = np.where(parameter_samples > args.sampling_constraints[i][1])[0]
				prior[index] = 0.
				posterior = likelihood * prior
			if args.sampling_method[i] == 'gaussian':
				prior = np.exp(-0.5*((args.sampling_constraints[i][0]-parameter_samples)**2.)/(args.sampling_constraints[i][1]**2.))
				posterior = likelihood*prior
			diff = np.exp(np.log(posterior)-np.log(likelihood))
			diff[np.isnan(diff)] = 0.
			diff[diff > 1] = 1.
			rand = np.random.rand(len(diff))
			index = np.where(diff > rand)[0]
			new_samples = new_samples[index, :]
			likelihood = likelihood[index]
	return new_samples
	
def initialise(module):
	import sys
	sys.path.append('params/')
	arguments = {'method': ['combine'], 'chain_dir': 'chains/', 'CMB': 'CMB', 'CMB_chains': 6, 'LSS': 'Strong_L', 'LSS_chains': 6, 'params': ['omegabh2', 'omegach2', 'theta', 'logA', 'ns'], 'bins': 40, 'plot_dir':None, 'save_dir': None, 'load': None,  'priors': None, 'smoothing': 0, 'interpolate_bins': None, 'integration_bounds': None, 'plot_interpolate': None, 'sampling_method': None, 'sampling_parameter': None, 'sampling_constraints': None}
	params = __import__(module)
	for keys in params.parameters:
		arguments[keys] = params.parameters[keys]
	print arguments
	args = parameters(arguments)
	return args

def interpolation_bins(CMB_ranges, LSS_ranges, num_params, num_bins = None):
	if num_bins != None:
		bins = np.array([num_bins for i in xrange(num_params)])
	else:
		domains = np.array([min(CMB_ranges[i][1]-CMB_ranges[i][0], LSS_ranges[i][1]-LSS_ranges[i][0]) for i in xrange(num_params)])
		bins = np.array([int(max((CMB_ranges[i][-1]-CMB_ranges[i][0])/domains[i], (LSS_ranges[i][-1]-LSS_ranges[i][0])/domains[i])) for i in xrange(num_params)])

	min_range = np.array([min(CMB_ranges[i][0], LSS_ranges[i][0]) for i in xrange(num_params)])
	max_range = np.array([max(CMB_ranges[i][-1], LSS_ranges[i][-1]) for i in xrange(num_params)])
	ranges = np.array([np.linspace(min_range[i], max_range[i], bins[i]) for i in xrange(num_params)])
	return ranges, bins

def load(filename):
	hist = np.load(filename + '_hist.npy')
	ranges = np.load(filename + '_ranges.npy')
	domainsize = np.load(filename + '_domainsize.npy')
	return hist, ranges, domainsize

def read_chains(args, dataset, dataset_num):
	columns = get_columns(args, dataset, args.params)
	data = []
	for num in xrange(int(dataset_num)):
		filename = args.chain_dir + dataset + "_" + str(num+1) + ".txt"
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
	samples = np.array(data)
	if args.sampling_method != None:
		samples = importance_sampling(args, samples, dataset, dataset_num)
	return samples

def save(filename, hist, ranges, domainsize):
	np.save(filename + '_hist.npy', hist)
	np.save(filename + '_ranges.npy', ranges)
	np.save(filename + '_domainsize.npy', domainsize)
