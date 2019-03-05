# Library of demographic models to be imported and tested
# All models are numbered the same as in Titus et al. (2019).

import moments as mts
import numpy as np

def model1(params, ns):
	"""
	Isolation only.
	"""
	nu1, nu2, T = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], T)
	return fs

def model2(params, ns):
	"""
	Isolation with population expansion.
	"""
	nu1_0, nu1, nu2_0, nu2, T = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])

	nu1_func = lambda t: nu1_0 * (nu1 / nu1_0) ** (t / T)
	nu2_func = lambda t: nu2_0 * (nu2 / nu2_0) ** (t / T)
	nu_func  = lambda t: [nu1_func(t), nu2_func(t)]
	fs.integrate(nu_func, T)
	return fs


def model3(params, ns):
	"""
	Isolation with symmetric mygration.
	"""
	nu1, nu2, T, m = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], T, m = np.array([[0,m],[m,0]]))
	return fs

def model4(params, ns):
	"""
	Isolation with symmetric migration and growth.
	"""
	nu1_0, nu1, nu2_0, nu2, T, m = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])

	nu1_func = lambda t: nu1_0 * (nu1 / nu1_0) ** (t / T)
	nu2_func = lambda t: nu2_0 * (nu2 / nu2_0) ** (t / T)
	nu_func  = lambda t: [nu1_func(t), nu2_func(t)]
	fs.integrate(nu_func, T, m = np.array([[0,m],[m,0]]))
	return fs

def model5(params, ns):
	"""
	Isolation with one-way migration: pop1 into pop2.

	This was the best model with Fastsimcoal.
	"""
	nu1, nu2, T, m12 = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], T, m = np.array([[0,m12],[0,0]]))
	return fs

def model6(params, ns):
	"""
	Isolation with one-way migration: pop2 into pop1.
	"""
	nu1, nu2, T, m21 = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], T, m = np.array([[0,0],[m21,0]]))
	return fs

def model7(params, ns):
	"""
	Isolation with initial migration, followed by no gene flow.
	"""
	nu1, nu2, Tdiv, Tiso, m = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], Tdiv, m = np.array([[0,m],[m,0]]))
	fs.integrate([nu1,nu2], Tiso, m = np.array([[0,0],[0,0]]))
	return fs

def model8(params, ns):
	"""
	Isolation followed secondary contact with symmetric migration.
	"""
	nu1, nu2, Tdiv, Tmig, m = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], Tdiv, m = np.array([[0,0],[0,0]]))
	fs.integrate([nu1,nu2], Tmig, m = np.array([[0,m],[m,0]]))
	return fs

def model9(params, ns):
	"""
	Isolation w/ initial migration: pop1 into pop2.
	"""
	nu1, nu2, Tdiv, Tiso, m12 = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], Tdiv, m = np.array([[0,m12],[0,0]]))
	fs.integrate([nu1,nu2], Tiso, m = np.array([[0,0],[0,0]]))
	return fs

def model10(params, ns):
	"""
	Isolation w/ secondary contact: pop1 into pop2.
	"""
	nu1, nu2, Tdiv, Tmig, m12 = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], Tdiv, m = np.array([[0,0],[0,0]]))
	fs.integrate([nu1,nu2], Tmig, m = np.array([[0,m12],[0,0]]))
	return fs

def model11(params, ns):
	"""
	Isolation w/ initial migration: pop2 into pop1.
	"""
	nu1, nu2, Tdiv, Tiso, m21 = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], Tdiv, m = np.array([[0,0],[m21,0]]))
	fs.integrate([nu1,nu2], Tiso, m = np.array([[0,0],[0,0]]))
	return fs

def model12(params, ns):
	"""
	Isolation w/ secondary contact: pop2 into pop1.
	"""
	nu1, nu2, Tdiv, Tmig, m21 = params
	sts = mts.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
	fs  = mts.Spectrum(sts)
	fs  = mts.Manips.split_1D_to_2D(fs, ns[0], ns[1])
	fs.integrate([nu1,nu2], Tdiv, m = np.array([[0,0],[0,0]]))
	fs.integrate([nu1,nu2], Tmig, m = np.array([[0,0],[m21,0]]))
	return fs
