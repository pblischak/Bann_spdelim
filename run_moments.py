#!/usr/bin/env python2

"""
<< run_moments.py >>

Usage:

python2 run_moments.py [input-sfs-file]

There are a couple of oddities that had to be dealt with in hacky ways.

First, the data are read in a strange way, such that values that are supposed
to be sero are either given really high values (like, 10^60) or really small
values (~10^-60). Not sure why (data seems to be formatted correctly), but we 
fix it by replacing them with 0.0.

Second, the optimization can get stuck in some weird places in paramter space,
gicing estimates near the border of the upper and lower bounds. We've implemented
iterative runs through the models so that if one of the parameters is wonky the
model gets rerun.
"""

from __future__ import print_function, division
import moments as mts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
from sys import argv, exit, stdout
import demographic_models as dm

def AIC(llik, k):
	"""
	Calculate AIC.
	"""
	return (2.0 * float(k)) - (2.0 * llik)

if __name__ == '__main__':
	"""
	This is super repetitive. Why didn't I write a function? Who knows.
	"""

	data = mts.Spectrum.from_file(argv[1]) # Read in data
	# There is something weird going on with how 0.0's are read in
	# so we need to reset obscure values to 0.0.
	data[data  < 1e-10] = 0.0
	data[data  > 1e10] = 0.0
	print(np.max(data))
	print(np.min(data))
	ns    = data.sample_sizes # Get samples sizes
	mu    = 4.38e-8 # Set mutation rate
	gtime = 0.001  # Generation time in thousands of years
	pfx   = argv[1].split(".")[0]
	try:
		f_out = open(pfx+".results", 'w') # Open output file
	except:
		print("Error: Couldn't open results file.")
		exit(-1)
	print("Model", "iter", "llik", "AIC", "params", "theta", sep="\t", file=f_out)


	###########################
	# Analyzing Model 1       #
	# Parameters: nu1, nu2, T #
	###########################
	iters = 1
	while(iters <= 10):
		print("** Running Model 1 (Iteration {0}) **".format(iters))
		func = dm.model1
		m1_params0 = [1.0, 1.0, 1.0] # Initial guess at parameters
		m1_upper   = [100.0, 100.0, 100.0] # Set upper bound
		m1_lower   = [1.0e-3, 1.0e-3, 1.0e-3] # Set lower bound
		m1_param0  = mts.Misc.perturb_params(m1_params0, fold=1,
											upper_bound = m1_upper,
											lower_bound = m1_lower) # perturb initial guess

		m1_poptg = mts.Inference.optimize_log(m1_params0, data, func,
											lower_bound = m1_lower,
											upper_bound = m1_upper,
											verbose=len(m1_params0), maxiter=10)
		m1_model    = func(m1_poptg, ns)
		m1_ll_model = mts.Inference.ll_multinom(m1_model, data)
		m1_theta    = mts.Inference.optimal_sfs_scaling(m1_model, data)

		m1_plot_model = mts.ModelPlot.generate_model(func, m1_poptg, ns)
		mts.ModelPlot.plot_model(m1_plot_model, save_file=pfx+"-Model1.png", draw_scale=False, nref=m1_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("1", iters, m1_ll_model, AIC(m1_ll_model, len(m1_poptg)), m1_poptg, m1_theta, sep="\t", file=f_out)
		iters += 1
		
	#mts.Plotting.plot_2d_comp_multinom(m1_model, data, vmin=1, resid_range=3)
	#pylab.savefig(pfx+"-Residuals1.png")
	
	#########################################
	# Analyzing Model 2                     #
	# Parameters: nu1_0, nu1, nu2_0, nu2, T #
	#########################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 2 (Iteration {0}) **".format(iters))
		func = dm.model2
		m2_params0 = [1.0, 1.0, 1.0, 1.0, 1.0] # Initial guess at parameters
		m2_upper   = [100.0, 100.0, 100.0, 100.0, 100.0] # Set upper bound
		m2_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3] # Set lower bound
		m2_params0 = mts.Misc.perturb_params(m2_params0, fold=1,
											upper_bound = m2_upper,
											lower_bound = m2_lower) # perturb initial guess
		m2_poptg = mts.Inference.optimize_log(m2_params0, data, func,
											lower_bound = m2_lower,
											upper_bound = m2_upper,
											verbose=len(m2_params0), maxiter=10)
		m2_model    = func(m2_poptg, ns)
		m2_ll_model = mts.Inference.ll_multinom(m2_model, data)
		m2_theta    = mts.Inference.optimal_sfs_scaling(m2_model, data)

		m2_plot_model = mts.ModelPlot.generate_model(func, m2_poptg, ns)
		mts.ModelPlot.plot_model(m2_plot_model, save_file=pfx+"-Model2.png", draw_scale=False, nref=m2_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		iters += 1
		print("2", iters, m2_ll_model, AIC(m2_ll_model, len(m2_poptg)), m2_poptg, m2_theta, sep="\t", file=f_out)

	#########################################
	# Analyzing Model 3                     #
	# Parameters: nu1, nu2, T, m            #
	#########################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 3 (Iteration {0}) **".format(iters))
		func = dm.model3
		m3_params0 = [1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m3_upper   = [100.0, 100.0, 100.0, 10.0] # Set upper bound
		m3_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3] # Set lower bound
		m3_params0 = mts.Misc.perturb_params(m3_params0, fold=1,
											upper_bound = m3_upper,
											lower_bound = m3_lower) # perturb initial guess
		m3_poptg = mts.Inference.optimize_log(m3_params0, data, func,
											lower_bound = m3_lower,
											upper_bound = m3_upper,
											verbose=len(m3_params0), maxiter=10)
		m3_model    = func(m3_poptg, ns)
		m3_ll_model = mts.Inference.ll_multinom(m3_model, data)
		m3_theta    = mts.Inference.optimal_sfs_scaling(m3_model, data)

		m3_plot_model = mts.ModelPlot.generate_model(func, m3_poptg, ns)
		mts.ModelPlot.plot_model(m3_plot_model, save_file=pfx+"-Model3.png", draw_scale=False, nref=m3_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("3", iters, m3_ll_model, AIC(m3_ll_model, len(m3_poptg)), m3_poptg, m3_theta, sep="\t", file=f_out)
		iters += 1

	############################################
	# Analyzing Model 4                        #
	# Parameters: nu1_0, nu1, nu2_0, nu2, T, m #
	############################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 4 (Iteration {0}) **".format(iters))
		func = dm.model4
		m4_params0 = [1.0, 1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m4_upper   = [100.0, 100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m4_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1e-3] # Set lower bound
		m4_params0 = mts.Misc.perturb_params(m4_params0, fold=1,
											upper_bound = m4_upper,
											lower_bound = m4_lower) # perturb initial guess
		m4_poptg = mts.Inference.optimize_log(m4_params0, data, func,
											lower_bound = m4_lower,
											upper_bound = m4_upper,
											verbose=len(m4_params0), maxiter=10)
		m4_model    = func(m4_poptg, ns)
		m4_ll_model = mts.Inference.ll_multinom(m4_model, data)
		m4_theta    = mts.Inference.optimal_sfs_scaling(m4_model, data)

		m4_plot_model = mts.ModelPlot.generate_model(func, m4_poptg, ns)
		mts.ModelPlot.plot_model(m4_plot_model, save_file=pfx+"-Model4.png", draw_scale=False, nref=m4_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("4", iters, m4_ll_model, AIC(m4_ll_model, len(m4_poptg)), m4_poptg, m4_theta, sep="\t", file=f_out)
		iters += 1

	################################
	# Analyzing Model 5            #
	# Parameters: nu1, nu2, T, m12 #
	################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 5 (Iteration {0}) **".format(iters))
		func = dm.model5
		m5_params0 = [1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m5_upper   = [100.0, 100.0, 100.0, 10.0] # Set upper bound
		m5_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1e-3] # Set lower bound
		m5_params0 = mts.Misc.perturb_params(m5_params0, fold=1,
											upper_bound = m5_upper,
											lower_bound = m5_lower) # perturb initial guess
		m5_poptg = mts.Inference.optimize_log(m5_params0, data, func,
											lower_bound = m5_lower,
											upper_bound = m5_upper,
											verbose=len(m5_params0), maxiter=10)
		m5_model    = func(m5_poptg, ns)
		m5_ll_model = mts.Inference.ll_multinom(m5_model, data)
		m5_theta    = mts.Inference.optimal_sfs_scaling(m5_model, data)

		m5_plot_model = mts.ModelPlot.generate_model(func, m5_poptg, ns)
		mts.ModelPlot.plot_model(m5_plot_model, save_file=pfx+"-Model5.png", draw_scale=False, nref=m5_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("5", iters, m5_ll_model, AIC(m5_ll_model, len(m5_poptg)), m5_poptg, m5_theta, sep="\t", file=f_out)
		iters += 1

	################################
	# Analyzing Model 6            #
	# Parameters: nu1, nu2, T, m21 #
	################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 6 (Iteration {0}) **".format(iters))
		func = dm.model6
		m6_params0 = [1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m6_upper   = [100.0, 100.0, 100.0, 10.0] # Set upper bound
		m6_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1e-3] # Set lower bound
		m6_params0 = mts.Misc.perturb_params(m6_params0, fold=1,
											upper_bound = m6_upper,
											lower_bound = m6_lower) # perturb initial guess
		m6_poptg = mts.Inference.optimize_log(m6_params0, data, func,
											lower_bound = m6_lower,
											upper_bound = m6_upper,
											verbose=len(m6_params0), maxiter=10)
		m6_model    = func(m6_poptg, ns)
		m6_ll_model = mts.Inference.ll_multinom(m6_model, data)
		m6_theta    = mts.Inference.optimal_sfs_scaling(m6_model, data)

		m6_plot_model = mts.ModelPlot.generate_model(func, m6_poptg, ns)
		mts.ModelPlot.plot_model(m6_plot_model, save_file=pfx+"-Model6.png", draw_scale=False, nref=m6_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("6", iters, m6_ll_model, AIC(m6_ll_model, len(m6_poptg)), m6_poptg, m6_theta, sep="\t", file=f_out)
		iters += 1

	#######################################
	# Analyzing Model 7                   #
	# Parameters: nu1, nu2, Tdiv, Tiso, m #
	#######################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 7 (Iteration {0}) **".format(iters))
		func = dm.model7
		m7_params0 = [1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m7_upper   = [100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m7_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1e-3, 1.0e-3] # Set lower bound
		m7_params0 = mts.Misc.perturb_params(m7_params0, fold=1,
											upper_bound = m7_upper,
											lower_bound = m7_lower) # perturb initial guess
		m7_poptg = mts.Inference.optimize_log(m7_params0, data, func,
											lower_bound = m7_lower,
											upper_bound = m7_upper,
											verbose=len(m7_params0), maxiter=10)
		m7_model    = func(m7_poptg, ns)
		m7_ll_model = mts.Inference.ll_multinom(m7_model, data)
		m7_theta    = mts.Inference.optimal_sfs_scaling(m7_model, data)

		m7_plot_model = mts.ModelPlot.generate_model(func, m7_poptg, ns)
		mts.ModelPlot.plot_model(m7_plot_model, save_file=pfx+"-Model7.png", draw_scale=False, nref=m7_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("7", iters, m7_ll_model, AIC(m7_ll_model, len(m7_poptg)), m7_poptg, m7_theta, sep="\t", file=f_out)
		iters += 1

	#######################################
	# Analyzing Model 8                   #
	# Parameters: nu1, nu2, Tdiv, Tmig, m #
	#######################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 8 (Iteration {0}) **".format(iters))
		func = dm.model8
		m8_params0 = [1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m8_upper   = [100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m8_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1e-3, 1.0e-3] # Set lower bound
		m8_params0 = mts.Misc.perturb_params(m8_params0, fold=1,
											upper_bound = m8_upper,
											lower_bound = m8_lower) # perturb initial guess
		m8_poptg = mts.Inference.optimize_log(m8_params0, data, func,
											lower_bound = m8_lower,
											upper_bound = m8_upper,
											verbose=len(m8_params0), maxiter=10)
		m8_model    = func(m8_poptg, ns)
		m8_ll_model = mts.Inference.ll_multinom(m8_model, data)
		m8_theta    = mts.Inference.optimal_sfs_scaling(m8_model, data)

		m8_plot_model = mts.ModelPlot.generate_model(func, m8_poptg, ns)
		mts.ModelPlot.plot_model(m8_plot_model, save_file=pfx+"-Model8.png", draw_scale=False, nref=m8_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("8", iters, m8_ll_model, AIC(m8_ll_model, len(m8_poptg)), m8_poptg, m8_theta, sep="\t", file=f_out)
		iters += 1

	#########################################
	# Analyzing Model 9                     #
	# Parameters: nu1, nu2, Tdiv, Tiso, m12 #
	#########################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 9 (Iteration {0}) **".format(iters))
		func = dm.model9
		m9_params0 = [1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m9_upper   = [100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m9_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1e-3, 1.0e-3] # Set lower bound
		m9_params0 = mts.Misc.perturb_params(m9_params0, fold=1,
											upper_bound = m9_upper,
											lower_bound = m9_lower) # perturb initial guess
		m9_poptg = mts.Inference.optimize_log(m9_params0, data, func,
											lower_bound = m9_lower,
											upper_bound = m9_upper,
											verbose=len(m9_params0), maxiter=10)
		m9_model    = func(m9_poptg, ns)
		m9_ll_model = mts.Inference.ll_multinom(m9_model, data)
		m9_theta    = mts.Inference.optimal_sfs_scaling(m9_model, data)

		m9_plot_model = mts.ModelPlot.generate_model(func, m9_poptg, ns)
		mts.ModelPlot.plot_model(m9_plot_model, save_file=pfx+"-Model9.png", draw_scale=False, nref=m9_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("9", iters, m9_ll_model, AIC(m9_ll_model, len(m9_poptg)), m9_poptg, m9_theta, sep="\t", file=f_out)
		iters += 1

	#########################################
	# Analyzing Model 10                    #
	# Parameters: nu1, nu2, Tdiv, Tmig, m12 #
	#########################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 10 (Iteration {0}) **".format(iters))
		func = dm.model10
		m10_params0 = [1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m10_upper   = [100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m10_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3] # Set lower bound
		m10_params0 = mts.Misc.perturb_params(m10_params0, fold=1,
											upper_bound = m10_upper,
											lower_bound = m10_lower) # perturb initial guess
		m10_poptg = mts.Inference.optimize_log(m10_params0, data, func,
											lower_bound = m10_lower,
											upper_bound = m10_upper,
											verbose=len(m10_params0), maxiter=10)
		m10_model    = func(m10_poptg, ns)
		m10_ll_model = mts.Inference.ll_multinom(m10_model, data)
		m10_theta    = mts.Inference.optimal_sfs_scaling(m10_model, data)

		m10_plot_model = mts.ModelPlot.generate_model(func, m10_poptg, ns)
		mts.ModelPlot.plot_model(m10_plot_model, save_file=pfx+"-Model10.png", draw_scale=False, nref=m10_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("10", iters, m10_ll_model, AIC(m10_ll_model, len(m10_poptg)), m10_poptg, m10_theta, sep="\t", file=f_out)
		iters += 1

	#########################################
	# Analyzing Model 11                    #
	# Parameters: nu1, nu2, Tdiv, Tiso, m21 #
	#########################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 11 (Iteration {0}) **".format(iters))
		func = dm.model11
		m11_params0 = [1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m11_upper   = [100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m11_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1e-3, 1.0e-3] # Set lower bound
		m11_params0 = mts.Misc.perturb_params(m11_params0, fold=1,
											upper_bound = m11_upper,
											lower_bound = m11_lower) # perturb initial guess
		m11_poptg = mts.Inference.optimize_log(m11_params0, data, func,
											lower_bound = m11_lower,
											upper_bound = m11_upper,
											verbose=len(m11_params0), maxiter=10)
		m11_model    = func(m11_poptg, ns)
		m11_ll_model = mts.Inference.ll_multinom(m11_model, data)
		m11_theta    = mts.Inference.optimal_sfs_scaling(m11_model, data)

		m11_plot_model = mts.ModelPlot.generate_model(func, m11_poptg, ns)
		mts.ModelPlot.plot_model(m11_plot_model, save_file=pfx+"-Model11.png", draw_scale=False, nref=m11_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("11", iters, m11_ll_model, AIC(m11_ll_model, len(m11_poptg)), m11_poptg, m11_theta, sep="\t", file=f_out)
		iters += 1

	#########################################
	# Analyzing Model 12                    #
	# Parameters: nu1, nu2, Tdiv, Tmig, m21 #
	#########################################
	f_out.flush() # Make sure output from previous analysis is written to file.
	iters = 1
	while(iters <= 10):
		print("** Running Model 12 (Iteration {0}) **".format(iters))
		func = dm.model12
		m12_params0 = [1.0, 1.0, 1.0, 1.0, 0.1] # Initial guess at parameters
		m12_upper   = [100.0, 100.0, 100.0, 100.0, 10.0] # Set upper bound
		m12_lower   = [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3] # Set lower bound
		m12_params0 = mts.Misc.perturb_params(m12_params0, fold=1,
											upper_bound = m12_upper,
											lower_bound = m12_lower) # perturb initial guess
		m12_poptg = mts.Inference.optimize_log(m12_params0, data, func,
											lower_bound = m12_lower,
											upper_bound = m12_upper,
											verbose=len(m12_params0), maxiter=10)
		m12_model    = func(m12_poptg, ns)
		m12_ll_model = mts.Inference.ll_multinom(m12_model, data)
		m12_theta    = mts.Inference.optimal_sfs_scaling(m12_model, data)

		m12_plot_model = mts.ModelPlot.generate_model(func, m12_poptg, ns)
		mts.ModelPlot.plot_model(m12_plot_model, save_file=pfx+"-Model12.png", draw_scale=False, nref=m12_theta/(4*mu), gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

		print("12", m12_ll_model, AIC(m12_ll_model, len(m12_poptg)), m12_poptg, m12_theta, sep="\t", file=f_out)
		iters += 1
