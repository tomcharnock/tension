from datetime import datetime
import os
import psutil
process = psutil.Process(os.getpid())
startTime = datetime.now()

import argparse
import numpy as np
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('-file', type = str, default = 'params')

args = initialise(parser.parse_args().file)

if args.plot_dir != None:
	from plot import *
	
if 'difference_vector' in args.method:
	hist, ranges, domainsize = difference_vector(args)
	C, tension = get_tension(hist, ranges, domainsize, args.num_params)
	print 'C = ',C
	print 'tension = ',tension,'sigma'
	if args.plot_dir != None:
		filename = get_filename(args.plot_dir, args.LSS, 'difference_vector', bins = args.bins, extra = str(tension))
		plot(args, ranges, domainsize, [hist], filename, plot_crosshairs = True)

else:
	CMB_hist, LSS_hist, shifted_hist, priors, ranges, domainsize = combine(args)
	output_string = ''
	if 'surprise' in args.method:
		D1_int, S1 = get_surprise(CMB_hist, LSS_hist, domainsize)
		D2_int, S2 = get_surprise(LSS_hist, CMB_hist, domainsize)
		print 'D1 = ', D1_int, ' and D2 = ', D2_int
		print 'S1 = ', S1, ' and S2 = ', S2	
		output_string += 'surprise_' + str(S1) + '_' + str(S2) + '_'
		if args.plot_dir != None:
			filename = get_filename(args.plot_dir, args.LSS, 'surprise', bins = args.interpolate_bins, extra = str(S1) + '_' + str(S2))
			plot_surprise(D1_int, S1, D2_int, S2, filename)
	if 'bhattacharyya' in args.method:
		P = np.sqrt(CMB_hist*LSS_hist)
        	B = np.sum(P*domainsize)
		print 'B = ',B
		output_string += 'bhattacharrya_' + str(B) + '_'
	if 'overlap_coefficient' in args.method:
		OVL = np.minimum(CMB_hist, LSS_hist)
		O = np.sum(OVL*domainsize)
		print 'O = ',O
		output_string += 'overlap_coefficient_' + str(O) + '_'
	if 'ibi' in args.method:
		I1 = get_ibi(CMB_hist, LSS_hist, domainsize, args.integration_bounds)
		I2 = get_ibi(LSS_hist, CMB_hist, domainsize, args.integration_bounds)
		print 'I1 = ', I1
		print 'I2 = ', I2
		output_string += 'ibi_' + str(I1) + '_' + str(I2) + '_'
	if ('marshall' in args.method):
		P = CMB_hist*LSS_hist
		logR = np.log(np.sum(P*priors*domainsize)/(np.sum(CMB_hist*priors*domainsize)*np.sum(LSS_hist*priors*domainsize)))
		print 'log R = ', logR
		output_string += 'logR_' + str(logR) + '_'
	if ('verde' in args.method):
		logT = np.log(np.sum(CMB_hist*shifted_hist*priors*domainsize)/np.sum(CMB_hist*LSS_hist*priors*domainsize))
		print 'log T = ', logT
		output_string += 'logT_' + str(logT) + '_'
	if args.plot_dir != None:
		output_string = output_string[:-1]
		filename = get_filename(args.plot_dir, args.LSS, output_string)
		plot(args, ranges, domainsize, [CMB_hist, LSS_hist], filename)

print datetime.now() - startTime, ' seconds.'
total_mem = process.memory_info().rss >> 20
print total_mem, 'MB'
