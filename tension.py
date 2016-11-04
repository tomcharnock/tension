import numpy as np
import sys
import matplotlib
#matplotlib.use('Agg') # Uncomment this when using on a hpc with no X server when plots need to be made
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class tension():
	def __init__(self, parameters):
		# Initialises the parameters and calculates the desired tension measure (and plots)
		self.load_parameters(parameters)								#	 - Loads parameters
		if self.method != None:
			if 'difference_vector' in self.method:
				self.difference_vector()							#	 - Calculates the difference vector and can make plots
			else:
				self.combine()									#	 - Calculates CMB and LSS histograms and other stuff
				output_string = ''								# String - Creates a string to save results to
				if 'surprise' in self.method:
					self.D1_int, self.S1 = self.get_surprise(self.CMB_hist, self.LSS_hist)	# Float  - Gets the relative entropy when LSS updates CMB
														# Float	 - Gets the surprise when LSS updates CMB
					self.D2_int, self.S2 = self.get_surprise(self.LSS_hist, self.CMB_hist)	# Float  - Gets the relative entropy when CMB updates LSS
														# Float  - Gets the surprise when CMB updates LSS
					print 'D1 = ', self.D1_int, ' and D2 = ', self.D2_int	
					print 'S1 = ', self.S1, ' and S2 = ', self.S2	
					output_string += 'surprise_' + str(self.S1) + '_' + str(self.S2) + '_'	# String - Adds the surprise results to the string
					if self.plot_dir != None:
						filename = self.get_filename(self.LSS, string = output_string)  # String - Gets the filename for the surprise only plot
						self.plot_surprise(filename)					#	 - Plots the surprises
				if 'bhattacharyya' in self.method:
					P = np.sqrt(self.CMB_hist*self.LSS_hist)				# Float  - Square root of the joint distributions
			        	self.B = np.sum(P*self.domainsize)					# Float	 - The Bhattacharyya distance
					print 'B = ', self.B
					output_string += 'bhattacharrya_' + str(self.B) + '_'			# String - Adds the Bhattacharyya distance to the string
				if 'overlap_coefficient' in self.method:
					OVL = np.minimum(self.CMB_hist, self.LSS_hist)				# Float  - Finds the minimum distribution of the CMB and LSS
					self.O = np.sum(OVL*self.domainsize)					# Float  - The overlap coefficient
					print 'O = ', self.O
					output_string += 'overlap_coefficient_' + str(self.O) + '_'		# String - Adds the overlap coefficient to the string
				if 'ibi' in self.method:
					self.I1 = self.get_ibi(self.CMB_hist, self.LSS_hist)			# Float	 - The integral of CMB between the bounds of LSS
					self.I2 = self.get_ibi(self.LSS_hist, self.CMB_hist)			# Float  - The integral of LSS betwwen the bounds of CMB
					print 'I1 = ', self.I1
					print 'I2 = ', self.I2
					output_string += 'ibi_' + str(self.I1) + '_' + str(self.I2) + '_'	# String - Adds the ibi results to the string
				if ('marshall' in self.method):
					self.logR = np.log(np.sum(self.CMB_hist*self.LSS_hist*self.priors*self.domainsize)\
						    /(np.sum(self.CMB_hist*self.priors*self.domainsize)*np.sum(self.LSS_hist*self.priors*self.domainsize)))
														# Float  - Joint distribution over the product of each distribution
					print 'log R = ', self.logR
					output_string += 'logR_' + str(self.logR) + '_'				# String - Adds logR result to the string
				if ('verde' in self.method):
					self.logT = np.log(np.sum(self.CMB_hist*self.shifted_hist*self.priors*self.domainsize)\
						    /np.sum(self.CMB_hist*self.LSS_hist*self.priors*self.domainsize))
														# Float  - Joint (shifted) distribution over the joint distribution
					print 'log T = ', self.logT		
					output_string += 'logT_' + str(self.logT) + '_'				# String - Adds logT result to the string
				if self.plot_dir != None:
					output_string = output_string[:-1]					# String - Removes the trailing '_' from the result string
					filename = self.get_filename(self.LSS, string = output_string)		# String - Gets the filename for the final plot
					self.plot(filename)							#	 - Plots each of the distributions

	def apply_smoothing(self):
		# Applies gaussian filter to histogram
		# Used in difference_vector() and get_histogram()
		from scipy.ndimage.filters import gaussian_filter
		smoothing = [self.smoothing for i in xrange(self.num_params)] #List  [num_params]               - Defines smoothing std for each dimension
		self.hist = gaussian_filter(self.hist, smoothing)	      #Array [bins for i in num_params] - Smooths histogram using gaussian filtering
 
	def combine(self):
		# Finds the histograms from CMB samples and LSS samples (and creates priors and shifted distribution for 'marshall' and 'verde')
		# Used in __init__()
		if self.load_file != None:					      
			self.load(self.load_file[0], 'CMB')		               # Array [interpolation_bins for i in num_params]     - Loads CMB histogram (might crash if
										       #			                	      dimensions not num_params and bins) 
			self.load(self.load_file[1], 'LSS')			       # Array [interpolation_bins for i in num_params]     - Loads LSS histogram (might crash if 
										       #       						      dimensions not num_params and bins)
										       # Array [bins, num_params]               	    - Loads LSS ranges
										       # Float				                    - Loads LSS domainsize
			self.hist = [self.CMB_hist, self.LSS_hist]		       # List [2, [interpolation_bins for i in num_params]] - List of CMB and LSS histograms
										       #						      for plotting
			if (('marshall' in self.method) or ('verde' in self.method)): 
				self.load(self.load_priors, 'priors')		       # Array [interpolation_bins for i in num_params]     - Loads priors (might crash if dimensions
										       #						      not num_params and bins)
				if 'verde' in self.method:			       
					self.load(self.load_file[1], 'shifted')	       # Array [bins for i in num_params]                   - Loads shifted histogram if using 'verde'
				else:
					self.shifted_hist = None                       # None				                    - No shifted_hist if not using 'verde'
			else:
				priors = None					       # None				      		    - No priors if not using 'verde', 'marshall'
				shifted_hist = None				       # None				      		    - No shifted_hist if not using 'verde'
		else:
			from itertools import product
			CMB_samples = self.get_samples(self.CMB, self.CMB_chains) # Array [number of samples, num_params]              - CMB samples
			self.get_histogram(CMB_samples)				  # Array [bins for i in num_params]	               - CMB histogram (global)
					  				          # Array [bins, num_params]		               - CMB ranges (global)
										  # Float				               - CMB domainsize (global)
			CMB_hist = self.hist					  # Array [bins for i in num_params]                   - Makes CMB histogram local variable
			CMB_ranges = self.ranges			          # Array [bins, num_params]		               - Makes CMB ranges local variable
			CMB_domainsize = self.domainsize                          # Float 				               - Makes CMB domainsize local variable
			LSS_samples = self.get_samples(self.LSS, self.LSS_chains) # Array [number of samples, num_params]              - LSS samples
			self.get_histogram(LSS_samples)				  # Array [bins for i in num_params]                   - LSS histogram (global)
										  # Array [bins, num_params]		               - LSS ranges (global)
 										  # Float				               - LSS domainsize (global)
			LSS_hist = self.hist					  # Array [bins for i in num_params]	               - Makes LSS histogram local variable
			LSS_ranges = self.ranges				  # Array [bins, num_params]		               - Makes LSS ranges local variable
			LSS_domainsize = self.domainsize			  # Float				               - Makes LSS domainsize local variable
			self.interpolation_bins(CMB_ranges, LSS_ranges)		  # Array [interpolation_bins, num_params]             - Gets ranges over whole CMB and LSS range
										  # List [num_params]			               - Interpolation bin numbers for each parameter
			self.get_domainsize()					  # Float 				               - Gets global domainsize
			domainsize = self.domainsize				  # Float				               - Makes domainsize local variable
			ranges = self.ranges					  # Array [interpolation_bins, num_params]             - Makes ranges local variable
			interp_ranges = list(product(*ranges))			  # itertools product			               - List of all ranges for interpolation
			self.CMB_hist = self.interpolate(interp_ranges, CMB_hist, CMB_ranges) # Array [interpolation_bins for i in num_params] - Gets interpolated CMB histogram
			self.LSS_hist = self.interpolate(interp_ranges, LSS_hist, LSS_ranges) # Array [interpolation_bins for i in num_params] - Gets interpolated LSS histogram
			if (('marshall' in self.method) or ('verde' in self.method)):	    
				if self.load_priors != None:
					self.load(self.load_priors, 'priors')		      # Array [interpolation_bins for i in num_params] - Loads priors (might crash if dims 
											      #						 	 not num_params and bins)
				else:
					self.get_priors()				      # Array [interpolation_bins for i in num_params] - Uniform priors in all parameters
				if 'verde' in self.method:
					shifted_hist, shifted_ranges = self.get_shifted(CMB_samples, LSS_samples)	  # Array [bins for i num_params] - Shifted histogram
															  # Array [bins, num_params]      - Shifted ranges
					self.shifted_hist = self.interpolate(interp_ranges, shifted_hist, shifted_ranges) # Array [interpolation_bins for i in num_params]
															  #                         - Interpolated shifted histogram
															  # Array [interpolation_bins, num_params]
					self.ranges = ranges		   # Array [interpolation_bins, num_params] - Resets global ranges to interpolated ranges
					self.domainsize = domainsize	   # Float				    - Resets global domainsize to interpolated domainsize
				else:
					self.shifted_hist = None	   # None				    - No shifted histogram needed for 'marshall'
			else:
				self.priors = None			   # None				    - No priors need when not using 'marshall' or 'verde'
				self.shifted_hist = None	 	   # None				    - No shifted histogram needed when not 'marshall' or 'verde'
			if self.save_dir != None:
				CMB_filename = self.get_filename(self.CMB) # String				    - CMB save filename
				LSS_filename = self.get_filename(self.LSS) # String				    - LSS save filename
				self.save(CMB_filename, 'CMB')		   #					    - Saves CMB histogram
				self.save(LSS_filename, 'LSS')		   #					    - Saves LSS histogram
				del CMB_hist, CMB_ranges, CMB_domainsize   #					    - Deletes large CMB variables from memory
				del LSS_hist, LSS_ranges, LSS_domainsize   #					    - Deletes large LSS variables from memory
				if (('marshall' in self.method) or ('verde' in self.method)):
					if self.load_priors == None:
						priors_filename = self.get_filename(self.LSS, string = 'priors')   # String	- Priors save filename
						self.save(priors_filename, 'priors')				   #		- Save priors
					if 'verde' in self.method:
						shifted_filename = self.get_filename(self.LSS, string = 'shifted') # String	- Shifted_hist save filename
						self.save(shifted_filename, 'shifted')				   #		- Saves shifted histogram
						del shifted_hist, shifted_ranges				   #		- Deletes large shifted variables from memory
			self.hist = [self.CMB_hist, self.LSS_hist]         # List [2, [interpolation_bins for i in num_params]] - List of CMB and LSS histograms for plotting
 
	def crosshairs(self, a, x, y):
		# Plots lines crossing at the origin of the plot for the difference_vector plots
		# Used in plot()
		self.ax[a].plot([x[0], x[-1]], [0, 0], color = 'black', linestyle = ':') # List [number of subplots]	- Plots horizontal line
		self.ax[a].plot([0, 0], [y[0], y[-1]], color = 'black', linestyle = ':') # List [number of subplots]	- Plots vertical line

	def define_hist(self, CMB, LSS):
		# Gets edges used to pre-set up for the histogram for the difference_vector
		# Used in difference_vector()
		bins = [[np.min(CMB[:, i])-np.max(LSS[:, i]), np.max(CMB[:, i])-np.min(LSS[:, i])] for i in xrange(self.num_params)] # Array [2, num_params]
																     #			- Lower/upper ranges
		a, edges = np.histogramdd([[0] for i in xrange(self.num_params)], bins = self.bins, range = bins)		     # Array [bins for i in num_params]
																     #			- Histogram (not used)
																     # Array [bins + 1, num_params]
																     #			- Edge of the bins
		return bins, edges

	#def difference_histogram(self, CMB, bins, LSS):				# Commented out because multiprocessing won't work with functions defined in classes
	#	diff = CMB-LSS
	#	histogram, edges = np.histogramdd(diff, bins = self.bins, range = bins)
	#	return histogram

	def difference_vector(self):
		# Finds histogram of the difference vector of the CMB and LSS samples 
		# Used in __init__()
		if self.load_file != None:
			self.load(self.load_file[0], 'difference_vector') 				# Array [bins for i in num_params]      - Loads histogram (might crash if
		 										        #	                    	          dimensions not num_params and bins)
												        # Array [bins, num_params]		- Loads ranges (might crash if
													#					  dimensions not num_params and bins)
													# Float				        - Loads domainsize
		else:
			import multiprocessing
			from multiprocessing import Pool
			from functools import partial
			CMB = self.get_samples(self.CMB, self.CMB_chains)				# Array [number of samples, num_params] - CMB samples
			LSS = self.get_samples(self.LSS, self.LSS_chains)				# Array [number of samples, num_params] - LSS samples
			bins, edges = self.define_hist(CMB, LSS)					# Array [2, num_params 			- Upper and lower parameter ranges
													# Array [bins + 1, num_params]		- Edge of the bins in histogram
			self.ranges = np.array([[edges[j][i]+(edges[j][i+1]-edges[j][i])/2. \
				      for i in xrange(self.bins)] for j in xrange(self.num_params)])    # Array [bins, num_params]		- Ranges
			self.get_domainsize()								# Float 				- Domainsize
			self.hist = np.zeros([self.bins for i in xrange(self.num_params)])		# Array [bins for i in num_params]	- Pre-makes the histogram
			print "Using ", multiprocessing.cpu_count(), " processors"
			pool = Pool()									#					- Initialises multiprocessing workers
			#histogram_pass = partial(self, difference_histogram, CMB, bins)		# Commented out because multiprocessing won't work with functions in classes
			histogram_pass = partial(difference_histogram, CMB, bins, self.bins)		# functool function			- Allows extra parameters to be
													# 					  passed to multiprocessing function
			for i in pool.imap_unordered(histogram_pass, LSS):				#					- Calculates histogram for each LSS
													#					- Samples as soon as a worker is free
				self.hist += i								# Array [bins for i in num_params]	- Populates pre-made histogram
			if self.smoothing != 0:
				self.apply_smoothing()							# Array [bins for i in num_params]      - Use gaussian smoothing on histogram
			self.hist = self.hist/(np.sum(self.hist)*self.domainsize)			# Array [bins for i in num_params]	- Normalises the histogram
			if self.save_dir != None:
				filename = self.get_filename(self.LSS)					# String				- Histogram save filename
				self.save(filename, 'difference_vector')                                #					- Saves difference vector histogram
		self.get_tension()									# 					- Gets the values of C and tension
		print "C = ", self.C
		print "tension = ", self.tension,"sigma"
		if self.plot_dir != None:
			self.hist = [self.hist]								# List					- Histogram into a list for plotting
			filename = self.get_filename(self.LSS, string = 'difference_vector_C_' + str(self.C)) # String				- Filename to save plot
			self.plot(filename, plot_crosshairs = True)					#					- Plot the difference vector contours

	def get_columns(self, chain, params):
		# Gets the columns for the parameters in the chains
		# Used in get_samples() and importance_sampling()
		columns = [0 for i in xrange(len(params))]					# Array [num_params]	- Place holder for column number
		filename = self.chain_dir + '/' + chain + ".paramnames"				# String		- Chain file name
		with open(filename) as f:
			i = 0									# Int			- Integer to track line number
			for read_line in f:
				line = filter(None, read_line.replace('\n','').split('\t'))	# String		- Line with newline removed and elements split by tab
				if line[0] in params:
					columns[params.index(line[0])] = int(i+2)		# Array [num_params]	- For each parameter (in turn) save the column number
				i += 1								# Int			- Adds one to the line number
		return columns

	def get_D_int(self, P1, P2):
		# Gets the value of the relative entropy
		# Used in get_surprise()
		D = P2*np.log(P2/P1)			# Array	[shape of P1 and P2]	- Calulates the state of the update from P1 to P2 at every point
		D[np.isnan(D)] = 0.			# Array [shape of P1 and P2]	- Removes any NaN values (not that there should be as P1 and P2 are set up properly) from D
		D_int = np.sum(D*self.domainsize)	# Float				- Gets the relative entropy
		return D_int

	def get_domainsize(self):
		# Gets the domainsize (d theta) which is the size of the each of the bins for each parameter
		# Used in combine(), difference_vector() and get_histogram()
		self.domainsize = 1									# Float	- Sets up the initial value of the domainsize
		for i in xrange(self.num_params):		
			self.domainsize = self.domainsize * (self.ranges[i][1] - self.ranges[i][0])	# Float - Multiplies the domainsize by the bin size for each parameter

	def get_filename(self, data, numerical = None, string = None):
		# Gets the filename to use for saving histograms, ranges, domainsizes and plots
		# Used in __init__(), combine() and difference_vector()
		if self.filename != None:
			filename = self.filename + '_' + data		# String - Name of filename if filename parameter is set with the data type appended
		else:
			filename = data					# String - Name of the filename is the data type (normally LSS parameter name apart from for CMB histogram)
		if numerical != None:
			filename = filename + '_' + str(numerical)	# String - Name of the filename appended by a numerical value (usually the number of bins)
		if string != None:
			filename = filename + '_' + string		# String - Name of the filename appended by a string (usually a descriptor)
		return filename

	def get_histogram(self, samples):
		# Creates the histogram from the supplied samples and normalises.
		# Used in combine() and get_shifted()
		self.hist, edges = np.histogramdd(samples, bins = self.bins)			# Array [interpolation_bins for i in num_params] - Creates the histogram
		if self.smoothing != 0:
			self.apply_smoothing()							# Array [interpolation_bins for i in num_params] - Smooths the histogram	
		self.ranges = np.array([[edges[j][i]+(edges[j][i+1]-edges[j][i])/2. \
			      for i in xrange(self.bins)] for j in xrange(self.num_params)])	# Array [interpolation_bins, num_params]	 - Defines the ranges 
		self.get_domainsize()								# Float						 - Gets the domainsize
		self.hist = self.hist/np.sum(self.hist)/self.domainsize				# Array [interpolation_bins for i in num_params] - Normalises the histogram

	def get_ibi(self, P1, P2):
		# Gets the integral of P1 between the integration bounds of P2
		# Used in combine()
		under_P2 = np.copy(P1)								# Array [interpolation_bins for i in num_params] - Copies P1
		under_P2[P2 < self.sigma(P2, self.integration_bounds/self.domainsize)] = 0	# Array [interpolation_bins for i in num_params] - Sets all values of P1 to zero at
												#						   values of P2 which are within the 
												#						   integration bounds
												#						   
		I = np.sum(under_P2*self.domainsize)						# Float						 - Integral of P1 with values outside
												#						   P2 integration bounds set to zero
		return I

	def get_labels(self):
		# Gets labels for each of the subplots
		# Used in plots()
		if self.method == 'difference_vector':
			label_dict = {'omegabh2': '$\Delta\Omega_{\\rm b}h^2$', 'omegach2': '$\Delta\Omega_{\\rm c}h^2$', 'theta': '$\Delta\Theta_{\\rm MC}$', 'logA': '$\Delta\log A_{\\rm s}$', 'ns': '$\Delta n_{\\rm s}$', 'mnu': '$\Delta\sum m_\\nu$', 'meffsterile': '$\Delta m_{\\rm eff}^{\\rm sterile}$', 'nnu': '$\Delta N_{\\rm eff}$'}
			# Dict	-  Labels for omegabh2, omegach2, theta, logA, ns, mnu, meffsterile and nnu when using the difference vector
		else:
			label_dict = {'omegabh2': '$\Omega_{\\rm b}h^2$', 'omegach2': '$\Omega_{\\rm c}h^2$', 'theta': '$\Theta_{\\rm MC}$', 'logA': '$\log A_{\\rm s}$', 'ns': '$n_{\\rm s}$', 'mnu': '$\sum m_\\nu$', 'meffsterile': '$m_{\\rm eff}^{\\rm sterile}$', 'nnu': '$N_{\\rm eff}$'}
			# Dict	-  Labels for omegabh2, omegach2, theta, logA, ns, mnu, meffsterile and nnu when not using the difference vector
		self.labels = [0 for i in xrange(self.num_params)]					 # List [num_params] - Place holder for the labels
		filename = self.chain_dir + '/' + self.CMB + ".paramnames"			 # String	     - Chain file name
		with open(filename) as f:	
			for read_line in f:
				line = filter(None, read_line.replace('\n','').split('\t'))	 # String	     - Line with newline removed and elements split by tab
				if line[0] in self.params:
					self.labels[self.params.index(line[0])] = label_dict[line[0]] # List [num_params] - List containing the label name ordered by self.parmams
	
	def get_priors(self):
		# Creates uniform priors the same shape as the histograms
		# Used in combine()
		self.priors = np.ones(self.shape)				# Array [interpolation_bins for i in num_params]	- Uniform prior for each parameter
		self.priors = self.priors/np.sum(self.priors*self.domainsize)	# Array [interpolation_bins for i in num_params]	- Normalises the prior

	def get_S(self, P1, P2, D_int):
		# Gets the value of the surprise
		# Used in get_surprise()
		D_av = np.mean(np.log(P1)-np.log(P2))	# Float		- Calculates the expected relative entropy
		S = D_int - np.abs(D_av)		# Float		- Value of the suprised calulated from the difference between the relative entropy and the expected entropy
		return S

	def get_samples(self, dataset, dataset_num):
		# Gets the samples from the chains and performs importance sampling if necessary
		# Used in combine() and difference_vector()
		columns = self.get_columns(dataset, self.params)			  # List [num_params]				  - Gets the columns of the parameters 
		samples = self.read_chains(dataset, dataset_num, columns)		  # Array [number of samples, num_params]	  - Gets the samples from chains
		if self.sampling_method != None:
			if ((self.sample_CMB) and (dataset == self.CMB)):
				samples = self.importance_sampling(samples, dataset, dataset_num) # Array [number of samples, num_params] - Importance samples the chains
			if ((self.sample_LSS) and (dataset == self.LSS)):
				samples = self.importance_sampling(samples, dataset, dataset_num) # Array [number of samples, num_params] - Importance samples the chains
		print "Number of samples in " + dataset + " = ",len(samples)
		print "Shape of data = ", samples.shape
		return samples
	
	def get_shifted(self, CMB_samples, LSS_samples):
		# Gets the shifted distribution when using 'verde'
		# Used in combine()
		CMB_means = np.mean(CMB_samples, axis = 0)		# Array [num_params]			         - Means of each of the parameters from the CMB_samples
		LSS_means = np.mean(LSS_samples, axis = 0)		# Array [num_params]			         - Means of each of the parameters from the LSS_samples
		shifted_samples = LSS_samples - LSS_means + CMB_means   # Array [number of LSS samples, num_params]      - Shift the LSS samples so they are centred on the CMB means
		self.get_histogram(shifted_samples)			# Array [interpolation_bins for i in num_params] - Gets the histogram of the shifted samples
		shifted_hist, shifted_ranges = self.hist, self.ranges	# Array [interpolation_bins for i in num_params] - Makes the shifted histogram a local variable
									# Array [interpolation_bins, num_params]	 - Makes the range of the shifted histogram a local variable
		return shifted_hist, shifted_ranges

	def get_surprise(self, P1, P2):
		# Get the surprise by calculating the relative entropy and expected entropy, preprocesses P1 and P2 so they won't give NaN values when calculating the entropy
		# Used in __init__()
		idx1 = np.where(P1 != 0)	# Array - Indices of where P1 is zero
		P1 = P1[idx1]			# Array - Removes zero values of P1
		P2 = P2[idx1]			# Array - Removes regions of P2 where P1 is zero
		idx2 = np.where(P2 != 0)	# Array - Indices of where P2 is zero
		P1 = P1[idx2]			# Array - Removes zero values of P2
		P2 = P2[idx2]			# Array - Removes regions of P1 where P2 is zero
		D_int = self.get_D_int(P1, P2)	# Float - Gets the value of the relative entropy
		S = self.get_S(P1, P2, D_int)	# Float - Gets the value of the surprise
		return D_int, S

	def get_tension(self):
		# Gets tension by calculating the value of the probability distribution at the origin of the difference of the parameters and integrating above that value
		# Used in difference_vector()
		from scipy.interpolate import interpn
		from scipy.special import erfinv
		val = interpn(self.ranges, self.hist, [0 for i in xrange(self.num_params)], bounds_error = False, fill_value = 0) # Float - Value of histogram above params = 0
		self.C = np.sum(self.hist[self.hist >= val]*self.domainsize)							  # Float - Integrates the histogram above val
		self.tension = np.sqrt(2.)*erfinv(self.C)									  # Float - Finds tension implied by the integration

	def importance_sampling(self, samples, dataset, dataset_num):
		# Importance sample the chains
		# Used in get_samples()
		likelihood = self.read_chains(dataset, dataset_num, [1])   # Array [number of samples, 1]         	 - Gets the likelihood from the chains
		likelihood = np.reshape(likelihood, likelihood.shape[0])   # Array [number of samples]            	 - Reshapes the likelihood into an array
		new_samples = np.copy(samples)				   # Array [number of samples, num_params]	 - Copies the samples
		index = np.where(new_samples[:, 0] != np.nan)[0]	   # Array [number of samples]			 - Gets all indices
		for i in xrange(len(self.sampling_method)):
			sample = True					   # Boolean					 - Boolean to decide whether to sample
			if self.sampling_parameter[i] not in self.params:
				column = self.get_columns(dataset, [self.sampling_parameter[i]])[0] # Int		 - Gets column of the parameter to importance sample		
				if column < 2:
					sample = False			   # Boolean				 	 - Do not sample if parameter is not in chains
					print "Cannot importance sample this chain because parameter is not present"
				else:
					parameter_samples = self.read_chains(dataset, dataset_num, [column])[index]
									   # Array [number of samples, 1]		 - Gets the parameter samples to importance sample
					parameter_samples = np.reshape(parameter_samples, parameter_samples.shape[0])
									   # Array [number of samples]		 	 - Reshapes the samples into an array
			else:	
				column = np.where(np.array(self.params) == self.sampling_parameter[i])[0][0]
									   # Int 				         - Get the column of the sampling parameter
				parameter_samples = new_samples[:, column] # Array [number of samples, num_params] 	 - Puts sampling parameter samples into new array
			if sample:
				if self.sampling_method[i] == 'uniform':
					prior = np.ones(len(parameter_samples))	# Array [number of samples]		 - Sets uniform prior
					index = np.where(parameter_samples < self.sampling_constraints[i][0])[0]
									   # Array					 - Gets the indices below minimum parameter value
					prior[index] = 0.		   # Array [number of samples]			 - Sets prior below minimum value to zero
					index = np.where(parameter_samples > self.sampling_constraints[i][1])[0]
									   # Array					 - Gets the indices above maximum parameter value
					prior[index] = 0.		   # Array [number of samples]			 - Sets prior above minimum value to zero
					posterior = likelihood * prior	   # Array [number of samples]			 - Calculates the posterior distribution
				if self.sampling_method[i] == 'gaussian':
					prior = np.exp(-0.5*((self.sampling_constraints[i][0]-parameter_samples)**2.)/(self.sampling_constraints[i][1]**2.))
									  # Array [number of samples]			 - Calculates prior distribution as a gaussian
					posterior = likelihood*prior	  # Array [number of samples]			 - Calculates the posterior distribution
				diff = np.exp(np.log(posterior)-np.log(likelihood)) # Array [number of samples]		 - Difference between the posterior and likelihood
				diff[np.isnan(diff)] = 0.		  # Array [number of samples]			 - Sets NaN values to zero
				diff[diff > 1] = 1.	  		  # Array [number of samples]			 - Sets values above one to one
				rand = np.random.rand(len(diff))	  # Array [number of samples]			 - Random numbers between zero and one for each value of diff
				index = np.where(diff > rand)[0]	  # Array [number of samples]			 - Indices where diff is greater than the random number
				new_samples = new_samples[index, :]	  # Array [number of samples, num_params]	 - The samples reassigned by the importance sampling
				likelihood = likelihood[index]		  # Array [number of samples]			 - Removes likelihood values which are counted out
		return new_samples

	def interpolate(self, interp_ranges, hist, ranges):
		# Interpolates the CMB or LSS (or shifted) histogram onto a histogram covering the whole parameter space
		# Used in combine()
		from scipy.interpolate import interpn
		hist_interp = interpn(ranges, hist, interp_ranges, bounds_error = False, fill_value = 0)		    
										# Array [interpolation_bins**num_params]         - Interpolation of the histogram
		hist_interp[hist_interp < 0] = 0				# Array [interpolation_bins**num_params]	 - Sets values less than zero to zero
		hist_interp = np.reshape(hist_interp, tuple([self.interpolation_bins[i] for i in xrange(self.num_params)]))
										# Array [interpolation_bins for i in num_params] - Reshapes the histogram
		hist_interp = hist_interp/np.sum(hist_interp)/self.domainsize   # Array [interpolation_bins for i in num_params] - Normalises the histogram
		self.shape = hist_interp.shape					# Tuple						 - Shape of the histogram
		return hist_interp

	def interpolation_bins(self, CMB_ranges, LSS_ranges):
		# Sets the bins for the interpolated histograms. This can be preset in the params file or done automatically (although this can use up HUGE amounts of memory)
		# Used in combine()
		if self.interpolate_bins != None:
			bins = np.array([self.bins for i in xrange(self.num_params)])
			self.interpolation_bins = bins
		else:
			domains = np.array([min(CMB_ranges[i][1]-CMB_ranges[i][0], LSS_ranges[i][1]-LSS_ranges[i][0]) for i in xrange(self.num_params)])
			self.interpolation_bins = np.array([int(max((CMB_ranges[i][-1]-CMB_ranges[i][0])/domains[i], (LSS_ranges[i][-1]-LSS_ranges[i][0])/domains[i])) \
						  for i in xrange(self.num_params)])
		min_range = np.array([min(CMB_ranges[i][0], LSS_ranges[i][0]) for i in xrange(self.num_params)])
		max_range = np.array([max(CMB_ranges[i][-1], LSS_ranges[i][-1]) for i in xrange(self.num_params)])
		self.ranges = np.array([np.linspace(min_range[i], max_range[i], self.interpolation_bins[i]) for i in xrange(self.num_params)])

	def load(self, filename, choice):
		# Loads the histograms, ranges, domainsize (and priors), (dv = difference_vector)
		# Used in combine() and difference_params()
		if choice == 'difference_vector':
			self.hist = np.load(self.save_dir + '/' + filename + '_difference_vector_hist.npy')	       # Array [bins for i in num_params] - Loads histogram (dv)
			self.ranges = np.load(self.save_dir + '/' + filename + '_difference_vector_ranges.npy')	       # Array [bins, num_params]         - Loads ranges (dv)
			self.domainsize = np.load(self.save_dir + '/' + filename + '_difference_vector_domainsize.npy')# Float			          - Loads domainsize (dv)
		if choice == 'CMB':
			self.CMB_hist = np.load(self.save_dir + '/' + filename + '_hist.npy')		 # Array [interpolation_bins for i in num_params] - Loads CMB histogram
		if choice == 'LSS':
			self.LSS_hist = np.load(self.save_dir + '/' + filename + '_hist.npy')		 # Array [interpolation_bins for i in num_params] - Loads LSS histogram
			self.ranges = np.load(self.save_dir + '/' + filename + '_ranges.npy')		 # Array [interpolation_bins, num_params]         - Loads interpolation ranges
			self.domainsize = np.load(self.save_dir + '/' + filename + '_domainsize.npy')	 # Float					  - Interpolation domainsize
		if choice == 'priors':
			self.priors = np.load(self.save_dir + '/' + filename + '_priors.npy')		 # Array [interpolation_bins for i in num_params] - Loads priors
		if choice == 'shifted':
			self.shifted_hist = np.load(self.save_dir + '/' + filename + '_shifted.npy')     # Array [interpolation_bins for i in num_params] - Loads shifted histogram

	def load_parameters(self, parameters):
		# Loads the class parameters
		# Used in __init__()
		arguments = {								 #			- DEFAULT PARAMETERS
			     'method': None,						 # List [num_params]	- List of strings for each tension quantification method
			     'chain_dir': 'chains',					 # String		- Directory where the chains are
			     'CMB': 'CMB',						 # String		- CMB chain file name
			     'CMB_chains': 6,						 # Int			- Number of CMB chains
			     'LSS': 'Strong_L',						 # String		- LSS chain file name
			     'LSS_chains': 6,						 # Int			- Number of LSS chains
 			     'params': ['omegabh2', 'omegach2', 'theta', 'logA', 'ns'],	 # List	[num_params]	- List of strings for each parameter in the histogram (COSMOMC naming
											 #			  style) also requires 'CMB' file parameter file (*.paramnames)
			     'bins': 40,						 # Int			- Number of bins to use in the histogram
			     'plot_dir': None,						 # String		- Directory where the plots should be saved
			     'save_dir': None,						 # String		- Directory where the data should be saved
			     'load_file': None,						 # List [1 or 2]	- List of strings with the names of the difference vector filename
											 #			  or the CMB filename and LSS filename (contained in 'save_dir')
			     'load_priors': None,				 	 # String		- Name of the priors filename (contained in 'save_dir')
			     'smoothing': 0.5,						 # Float		- Amount of smoothing the gaussian filter applies to the histograms
			     'interpolate_bins': None,					 # Int			- Number of bins to interpolate the CMB and LSS histograms on to
			     'integration_bounds': None,				 # Float		- The fraction of samples to integrate to form the integration bounds
											 #			  in the integral between intervals method
			     'plot_show': False,					 # Boolean		- Show plots when they are made (True) or just save them (False)
			     'sampling_method': None,					 # List	[sample params]	- Strings with 'uniform' or 'gaussian' for each sampling parameter
			     'sampling_parameter': None,				 # List	[sample params]	- Strings containing each of the sampling parameters, the parameter 
											 #		          will be ignored if it's not in the MCMC chains
			     'sampling_constraints': None,				 # List [2 for i in sample params]
											 #			- List containing list of lower and upper bounds for 'uniform' or
											 #			- mean and standard deviation for 'gaussian'
			     'sample_CMB': False,					 # Boolean		- Whether to importance sample the CMB chains
			     'sample_LSS': False,					 # Boolean		- Whether to importance sample the LSS chains
			     'filename': None}						 # String		- Filename at the beginning of each saved files and plots
		for keys in parameters:
			arguments[keys] = parameters[keys]				 #			- Loads parameters from dictionary
		self.method = arguments['method']					 # List [num_params]	- List of strings for each tension quantification method
		self.chain_dir = arguments['chain_dir']                                  # String		- Directory where the chains are
		self.CMB = arguments['CMB']                                              # String		- CMB chain file name
		self.CMB_chains = arguments['CMB_chains']                                # Int			- Number of CMB chains
		self.LSS = arguments['LSS']                                              # String		- LSS chain file name
		self.LSS_chains = arguments['LSS_chains']                                # Int			- Number of LSS chains
		self.params = arguments['params']                                        # List	[num_params]	- List of strings for each parameter in the histogram (COSMOMC naming
                                                                                         #			  style) also requires 'CMB' file parameter file (*.paramnames)
		self.num_params = len(self.params)                                       # Int			- Number of parameters 
		self.bins = arguments['bins']                                            # Int			- Number of bins to use in the histogram
		self.plot_dir = arguments['plot_dir']                                    # String		- Directory where the plots should be saved
		self.save_dir = arguments['save_dir']                                    # String		- Directory where the data should be saved
		self.load_file = arguments['load_file']                                  # List [1 or 2]	- List of strings with the names of the difference vector filename
		self.load_priors = arguments['load_priors']                              #			  or the CMB filename and LSS filename (contained in 'save_dir')
                                                                                         # String		- Name of the priors filename (contained in 'save_dir')
		self.smoothing = arguments['smoothing']                                  # Float		- Amount of smoothing the gaussian filter applies to the histograms
		self.interpolate_bins = arguments['interpolate_bins']                    # Int			- Number of bins to interpolate the CMB and LSS histograms on to
		self.integration_bounds = arguments['integration_bounds']                # Float		- The fraction of samples to integrate to form the integration bounds
                                                                                         #			  in the integral between intervals method
		self.plot_show = arguments['plot_show']                                  # Boolean		- Show plots when they are made (True) or just save them (False)
                self.sampling_method = arguments['sampling_method']                      # List	[sample params]	- Strings with 'uniform' or 'gaussian' for each sampling parameter
		self.sampling_parameter = arguments['sampling_parameter']                # List	[sample params]	- Strings containing each of the sampling parameters, the parameter 
                                                                                         #		          will be ignored if it's not in the MCMC chains
                self.sampling_constraints = arguments['sampling_constraints']            # List [2 for i in sample params]
                                                                                         #			- List containing list of lower and upper bounds for 'uniform' or
                                                                                         #			- mean and standard deviation for 'gaussian'
                self.sample_CMB = arguments['sample_CMB']            		         # Boolean		- Whether to importance sample the CMB chains
		self.sample_LSS = arguments['sample_LSS']				 # Boolean		- Whether to importance sample the LSS chains
		self.filename = arguments['filename']                                    # String		- Filename at the beginning of each saved files and plots
		if self.method != None:
			if not 'difference_vector' in self.method:
				if self.interpolate_bins != None:
					if self.interpolate_bins < 1:
						print 'Number of interpolate bins too low'
						sys.exit()
				if (('marshall' in self.method) or ('verde' in self.method)):
					if self.load_file != None:
						if self.load_priors == None:
							print 'Need to set the priors filename to load from'
							sys.exit()
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

	def plot(self, filename, plot_crosshairs = False):
		# Plots the triangle plots for each of the parameters
		# Used in differenc_vector() and __init__()
		self.rows = self.num_params-1	  				 # Int		       - Number of rows in the subplot
		self.columns = self.num_params-1  				 # Int		       - Number of columns in the subplot
		self.plot_setup()		  				 # List	               - List of the axes instances in the subplot
		self.get_labels()		  				 # List	               - List of the labels for each of the axes
		sci = self.sci_ranges()		  				 # List	               - List of the exponents for the ranges
		a = 0				  				 # Int		       - Axes tracker
		for row in xrange(self.rows):
			for column in xrange(row, self.columns):
				axes = tuple([k for k in xrange(self.num_params) if ((k != row) and (k != column+1))]) 
										  # Tuple	       - Tuple of the axes to sum over to give the 2D contour to plot in each subplot
				if row == column:
					if sci[a][0] != 0:
						self.ax[a].set_xlabel(self.labels[column+1] + ' [$\\times10^{'+ '{0:0d}'.format(int(sci[a][0])) + '}$]', labelpad=0)
										  #		       - Sets x label using parameter labels and the exponent
					else:
						self.ax[a].set_xlabel(self.labels[column+1], labelpad=0)
										  #		       - Sets x label using parameter labels without exponent if exponent is zero
						
					if sci[a][1] != 0:
						self.ax[a].set_ylabel(self.labels[row] + ' [$\\times10^{'+ '{0:0d}'.format(int(sci[a][1])) + '}$]', labelpad=0)
										  #		       - Sets y label using parameter labels and the exponent
					else:
						self.ax[a].set_ylabel(self.labels[row], labelpad=0)
										  #		       - Sets y label using parameter labels without exponent if exponent is zero
				domains = (self.ranges[column+1, 1] - self.ranges[column+1, 0])*(self.ranges[row, 1] - self.ranges[row, 0])
						  				  #		       - Gets the domains for the relevant parameter ranges
				colours = ['blue', 'red'] 			  # List [2]   - Colour of plots
				linestyles = ['-', '--'] 			  # List [2]   - Linestyles for the plots
				for k in xrange(len(self.hist)):
					P_plot = np.sum(self.hist[k]*self.domainsize/domains, axis = axes) 
						                                  # Array [bins, bins] - Sum each of the parameter dimensions to get two dimensional histogram
					x = self.ranges[column+1]/(10**sci[a][0]) # Array [bins]       - x axis values scaled by the exponent
					y = self.ranges[row]/(10**sci[a][1])      # Array [bins]       - y axis values scaled by the exponent
					P_levels = [self.sigma(P_plot, 0.997/domains), self.sigma(P_plot, 0.95/domains), self.sigma(P_plot, 0.68/domains)]
										  # List [3]	       - The  values of the histogram corresponding to the 1, 2 and 3 sigma contours
					if not all(q<r for q, r in zip(P_levels, P_levels[1:])):
						if len(self.hist) == 1: 
							self.ax[a].contourf(x, y, P_plot) #	       - Plots filled contour of the histogram when P_levels aren't in rising order
											  #	         (this can happen when there are not enough samples)
							if plot_crosshairs:
								self.crosshairs(a, x, y)  #	       - Plots lines passing through the origin (for difference vector)
						else:
							self.ax[a].contourf(x, y, P_plot) #	       - Plots filled contour of the histogram when P_levels aren't in rising order
											  #		 (this can happen whenthere are not enough samples)
					else:
						if len(self.hist) == 1: 
							self.ax[a].contour(x, y, P_plot, colors = 'purple', levels = P_levels)
										  #		       - Plots difference vector histograms with 1, 2 and 3 sigma contours
							if plot_crosshairs:
								self.crosshairs(a, x, y)  #	       - Plots lines passing through the origin (for difference vector)	
						else:
							self.ax[a].contour(x, y, P_plot, colors = colours[k], linestyles = linestyles[k], levels = P_levels)
										  #		       - Plots 1, 2 and 3 sigma contours for CMB and LSS histograms, blue solid lines
										  #			 for CMB and red dashed lines for LSS
					self.ax[a].locator_params(nbins=5)	  #		       - Sets the number of bins for each axis
				a += 1						  # Int		       - Adds one to the axis tracker
		self.plot_save(filename)					  #		       - Saves the plot

	def plot_save(self, filename):
		# Saves the plots
		# Used in plot() and plot_surprise()
		plt.savefig(self.plot_dir + '/' + filename + '.pdf', bbox_inches='tight')	#	- Saves the figure as a pdf and remove excess whitespace
		if self.plot_show:
			plt.show()								#	- Shows the plot
		plt.close()									#	- Closes the plot

	def plot_setup(self):
		# Sets up the axis subplots for plotting
		# Used in plot() and plot_surprise()
		from matplotlib.gridspec import GridSpec
		plot_params = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'legend.fontsize': 8, 'font.size': 8} 
											#	- Plot style parameters
		plt.rcParams.update(plot_params)					#	- Sets plot style parameters	
		fig = plt.figure(figsize=((3.5/3)*self.columns, 1.5*self.rows))		# 	- Makes figure big enough for each of the subplots
		gs = GridSpec(self.rows, self.columns)					#	- Sets the number of subplots
		gs.update(wspace = 0, hspace=0, bottom = 0.2)				#	- Removes space between plots
		if ((self.rows == 1) and (self.columns == 1)):			
			self.ax = [plt.subplot(gs[0, 0])]				# List	- Adds subplot to list
		else:
			self.ax = []							# List	- Empty list place holder
			for row in xrange(self.rows):
				for column in xrange(row, self.columns):
					self.ax.append(plt.subplot(gs[row, column]))	# List  - Adds each subplot to the list
					if row != column:
						self.ax[-1].set_xticks([])		#	- Removes x tick labels from internal subplots
						self.ax[-1].set_yticks([])		#	- Removes y tick labels from interpal subplots

	def plot_surprise(self, filename):
		# Plot the surprise on a single subplot
		# Used in __init__()
		self.columns = 1					 # Int       - Number of columns in subplot
		self.rows = 1						 # Int       - Number of columns in subplot
		self.plot_setup()					 # List      - List of the axes in the plot
		ax = self.ax[0]						 # Axis      - The subplot axis
		ax.plot([0, 0], [0, 3], color = 'black')		 #           - Plots the origin line
		ax.barh(2, self.D1_int, 0.8, color = 'blue')		 #           - Plots the relative entropy when P1 updates P2
		ax.barh(2.1, self.S1, 0.6,color = 'dodgerblue')		 #           - Plots the surprise when P1 updates P2
		ax.barh(1, self.D2_int, 0.8, color = 'firebrick')	 #           - Plots the relative entropy when P2 updates P1
		ax.barh(1.1, self.S2, 0.6, color = 'red')		 #           - Plots the surprise when P2 updates P1
		xmin = min(0, np.min(self.D1_int))			 # Float     - Finds the minimum value of the relative entropies and the surprises 
		xmin = min(xmin, np.min(self.D2_int))			 # 
		xmin = min(xmin, np.min(self.S1))			 #
		xmin = min(xmin, np.min(self.S2))			 #
		xmax = max(0, np.max(self.D1_int))			 # Float     - Finds the maximum value of the relative entropies and the surprises
		xmax = max(xmax, np.max(self.D2_int))			 #
		xmax = max(xmax, np.max(self.S1))			 #
		xmax = max(xmax, np.max(self.S2))			 #
		lim = []						 # List	     - Empty place holder for the x limits
		if xmin < 0:
			lim.append(xmin)				 # List [1]  - List with the minimum x limit
			xmin = xmin + 0.1 * xmin			 # Float     - Moves the minimum a little lower
		if xmax > 0:
			lim.append(xmax)				 # List [2]  - List with the minimum and maximum x limits
			xmax = xmax + 0.1 * xmax			 # Float     - Moves the maximum a little higher
		if (xmin is 0) and (xmax is 0):
			xmin = -1					 # Float     - Sets the minimum x limit to -1 if relative entropies and surprises are zero
			xmax = 1					 # Float     - Sets the maximum x limit to 1 if relative entropies and surprises are zero
			lim.append(0)					 # List [1]  - List with zero
		ax.set_xticks(lim)					 # 	     - Sets the x tick labels
		ax.set_xlim([xmin, xmax])				 #	     - Sets the x limits
		ax.set_ylim([0.8, 3])					 #	     - Sets the y limits
		ax.set_yticks([])					 #	     - Removes the y tick labels
		ax.set_xlabel('${\\rm Bits}\\times\log2$', labelpad = 0) #	     - Sets the x label
		self.plot_save(filename)				 #	     - Saves the plot

	def read_chains(self, dataset, dataset_num, columns):
		# Read in the samples from the chains
		# Used in get_samples() and importance_sampling()
		data = []										   # List		         - Empty place holder for the samples
		for num in xrange(int(dataset_num)):				
			filename = self.chain_dir + '/' + dataset + "_" + str(num+1) + ".txt"		   # String	                 - Reads in each of the MCMC chains
			num_lines = 0									   # Int		         - Line number tracker
			with open(filename) as f:
				for read_line in f:
					num_lines += 1							   # Int		         - Adds one to the number tracker
			skip_lines = int(num_lines/3.)							   # Int	 	         - Number of lines in the chain to skip
			num_lines = 0									   # Int			 - Resets the number tracker
			with open(filename) as f:
				for read_line in f:
					if num_lines > skip_lines:
						line = filter(None, read_line.replace('\n','').split(' ')) # String	                 - Line with newtab removed, split ' '
						repeat = int(float(line[0]))				   # Int			 - Gets the weight of each set of samples
						for i in xrange(repeat):
							data.append([float(line[j]) for j in columns])	   # List [num_samples, columns] - Adds the samples of each parameter 
					num_lines += 1							   # Int			 - Adds one to the number tracker
		return np.array(data)

	def save(self, filename, choice):
		# Saves the histograms, ranges, domainsize (and priors)
		# Used in combine() and difference_vector()
		if choice == 'difference_vector':
			np.save(self.save_dir + '/' + filename + '_difference_vector_hist.npy', self.hist)		# 	- Saves difference vector histogram
			np.save(self.save_dir + '/' + filename + '_difference_vector_ranges.npy', self.ranges)		# 	- Saves difference vector ranges
			np.save(self.save_dir + '/' + filename + '_difference_vector_domainsize.npy', self.domainsize)  # 	- Saves difference vector domainsize
		if choice == 'CMB':
			np.save(self.save_dir + '/' + filename + '_hist.npy', self.CMB_hist)				#	- Saves CMB histogram
		if choice == 'LSS':
			np.save(self.save_dir + '/' + filename + '_hist.npy', self.LSS_hist)				# 	- Saves LSS histogram
			np.save(self.save_dir + '/' + filename + '_ranges.npy', self.ranges)				#	- Saves LSS ranges
			np.save(self.save_dir + '/' + filename + '_domainsize.npy', self.domainsize)			#	- Saves LSS domainsize
		if choice == 'priors':
			np.save(self.save_dir + '/' + filename + '.npy', self.priors)					#	- Saves priors
		if choice == 'shifted':
			np.save(self.save_dir + '/' + filename + '.npy', self.shifted_hist)				#	- Saves shifted histogram


	def sci_lims(self, x, y):
		# Gets the power to which the axis labels are raised to and saves it so that the values can be made neater
		# Used in sci_ranges()
		fig = plt.figure()						#		- Makes new figure
		ax = fig.add_subplot(111)					#		- Adds a subplot
		ax.plot([x[0], x[-1]], [y[0], y[-1]])				#		- Plots a line from the minimum of each range to the maximum
	 	ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))	#		- Gets ticklabels in a scientific form
		plt.tight_layout()						#		- Sets tight_layout (seems to be a bug in matplotlib which means this is needed)
		xoffset = ax.get_xaxis().get_offset_text().get_text()		# String	- The x-axis exponent text
		xoffset = xoffset.replace('1e', '')				# String	- Removes the exponent from the beginning of the string
		xoffset = xoffset.replace(u'\u2212', '-')			# String	- Removes any unicode characters as they cannot be interpretted
		if xoffset == '':
			xoffset = '0'						# String	- Returns zero string if there is no exponent
		xoffset = float(xoffset)					# Float		- Turns the string to a float
		yoffset = ax.get_yaxis().get_offset_text().get_text()		# String	- The x-axis exponent text
		yoffset = yoffset.replace('1e', '')                  		# String	- Removes the exponent from the beginning of the string
		yoffset = yoffset.replace(u'\u2212', '-')            		# String	- Removes any unicode characters as they cannot be interpretted
		if yoffset == '':                                    		                                                                               
			yoffset = '0'                                		# String	- Returns zero string if there is no exponent
		yoffset = float(yoffset)                             		# Float		- Turns the string to a float
		plt.close()							#		- Closes the figure
		return [xoffset, yoffset]

	def sci_ranges(self):
		# Gets the power to which the axis labels are raised for each parameter in the plot
		# Used in plot()
		sci = []									# List		- Empty list for containing the axis powers
		for i in xrange(self.rows):
			for j in xrange(i, self.columns):
				sci.append(self.sci_lims(self.ranges[j+1], self.ranges[i]))	# List		- Adds the axis powers to the list for the ranges of each subplot
		return sci

	def sigma(self, P, l):
		# Calculates the value of the probability distribution above which the integral of the gives the integration bounds (l)
		# Used in get_ibi() and plot()
		P_s = np.sort(P.ravel())[::-1] # Array [bins**num_parms]	- Sorts the probability distribution from highest to lowest
		sigma = 0		       # Float				- Tracks the sum of the values from the probability distribution	
		for i in P_s:			
			sigma += i	       # Float				- For each value in the sorted probability distribution add it to the tracker
			if sigma >= l:
				return i       #				- When the value of the tracker reaches the integration bounds, return the value of the histogram
		return 1	

def difference_histogram(CMB, bins, num_bins, LSS):
	# Calculates the histogram of the difference vector
	# Used in difference_vector()
	diff = CMB-LSS								# Array [length of CMB, num_params] - Difference between all the CMB samples from a single LSS
	histogram, edges = np.histogramdd(diff, bins = num_bins, range = bins)  # Array [bins for i in num_params]  - Histogram of the difference
	return histogram

if __name__ == "__main__":
	# Used when using tension.py as a standalone script with parameter file in params/
	from datetime import datetime
	import os
	import psutil
	process = psutil.Process(os.getpid()) 				# Int			- Gets the process id
	startTime = datetime.now()					# Float 		- Saves the time when starting the script
	import argparse

	parser = argparse.ArgumentParser()				# argparse parser	- Creates a parser to read in the parameter file
	parser.add_argument('-file', type = str, default = 'params')	# 			- Gets file name to read in (default is params/params.py)
	sys.path.append('params/')					# 			- Adds params/ to the system path to load parameter file
	params = __import__(parser.parse_args().file)			# 			- Loads the parameter file as a module
	tension(params.parameters)					# 			- Initialises the class and run the tension script

	print "Time elapsed ", datetime.now() -startTime			#			- Prints the amount of time the script took
	total_mem = process.memory_info().rss >> 20			# Float			- Gets the amount of memory the script used 
	print "Memory used ", total_mem, 'MB'						#			- Prints the total memory in MB
	
