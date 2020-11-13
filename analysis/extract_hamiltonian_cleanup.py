# -*- coding: utf-8 -*-
"""See document on the wiki.

Created on Sun Jun 16 11:22:40 2019

@author: mccambria
"""


# %% Imports


import numpy
from numpy.linalg import eigvals
from numpy import pi
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import brute
from numpy import exp
import matplotlib.pyplot as plt


# %% Constants


d_gs = 2.87  # ground state zfs in GHz
gmuB = 2.8  # gyromagnetic ratio in MHz / G
gmuB_GHz = gmuB / 1000  # gyromagnetic ratio in GHz / G

# numbers
inv_sqrt_2 = 1/numpy.sqrt(2)
im = 0+1j


# %% Functions


def calc_single_hamiltonian(mag_B, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi):
    par_B = mag_B * numpy.cos(theta_B)
    perp_B = mag_B * numpy.sin(theta_B)
    hamiltonian = numpy.array([[d_gs + par_Pi + par_B,
                                inv_sqrt_2 * perp_B * exp(-1j * phi_B),
                                -perp_Pi * exp(1j * phi_Pi)],
                               [inv_sqrt_2 * perp_B * exp(1j * phi_B),
                                0,
                                inv_sqrt_2 * perp_B * exp(-1j * phi_B)],
                               [-perp_Pi * exp(-1j * phi_Pi),
                                inv_sqrt_2 * perp_B * exp(1j * phi_B),
                                d_gs + par_Pi - par_B]])
    return hamiltonian


def calc_hamiltonian(mag_B, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi):
    fit_vec = [theta_B, par_Pi, perp_Pi, phi_B, phi_Pi]
    if (type(mag_B) is list) or (type(mag_B) is numpy.ndarray):
        hamiltonian_list = [calc_single_hamiltonian(val, *fit_vec)
                            for val in mag_B]
        return hamiltonian_list
    else:
        return calc_single_hamiltonian(mag_B, *fit_vec)


def calc_res_pair(mag_B, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi):
    hamiltonian = calc_hamiltonian(mag_B, theta_B, par_Pi, perp_Pi,
                                   phi_B, phi_Pi)
    if (type(mag_B) is list) or (type(mag_B) is numpy.ndarray):
        vals = numpy.sort(eigvals(hamiltonian), axis=1)
        resonance_low = numpy.real(vals[:,1] - vals[:,0])
        resonance_high = numpy.real(vals[:,2] - vals[:,0])
    else:
        vals = numpy.sort(eigvals(hamiltonian))
        resonance_low = numpy.real(vals[1] - vals[0])
        resonance_high = numpy.real(vals[2] - vals[0])
    return resonance_low, resonance_high


def find_mag_B(res_desc, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi):
    # Just return the given mag_B if it's known
    if res_desc[0] is not None:
        return res_desc[0]
    # Otherwise we'll determine the most likely mag_B for this fit_vec by
    # finding the mag_B that minimizes the distance between the measured
    # resonances and the calculated resonances for a given fit_vec
    args = (res_desc, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi)
    result = minimize_scalar(find_mag_B_objective, bounds=(0, 1.0), args=args,
                             method='bounded')
    if result.success:
        mag_B = result.x
    else:
        # If we didn't find an optimal value, return something that will blow
        # up chisq and push us away from this fit_vec
        mag_B = 0.0
    return mag_B


def find_mag_B_objective(x, res_desc, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi):
    calculated_res_pair = calc_res_pair(x, theta_B, par_Pi, perp_Pi,
                                        phi_B, phi_Pi)
    diffs = numpy.array(calculated_res_pair) - numpy.array(res_desc[1:3])
    sum_squared_differences = numpy.sum(diffs**2)
    return sum_squared_differences


def plot_resonances(mag_B_range, theta_B, par_Pi, perp_Pi, phi_B, phi_Pi,
                    name='untitled'):

    smooth_mag_B = numpy.linspace(mag_B_range[0], mag_B_range[1], 1000)
    res_pairs = calc_res_pair(smooth_mag_B, theta_B, par_Pi, perp_Pi,
                              phi_B, phi_Pi)

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    fig.set_tight_layout(True)
    ax.set_title('Generating fit vector: {}'.format(name))
    ax.plot(smooth_mag_B, res_pairs[0])
    ax.plot(smooth_mag_B, res_pairs[1])
    ax.set_xlabel('B magnitude (GHz)')
    ax.set_ylabel('Resonance (GHz)')

    textstr = '\n'.join((
        r'$\theta_{B}=%.3f \ rad$' % (theta_B, ),
        r'$\Pi_{\parallel}=%.3f \ GHz$' % (par_Pi, ),
        r'$\Pi_{\perp}=%.3f \ GHz$' % (perp_Pi, )))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, textstr, fontsize=14, transform=ax.transAxes,
            verticalalignment='top', bbox=props)

    return fig, ax


def chisq_func_reduced(fit_vec, par_Pi, perp_Pi, phi_B, phi_Pi, res_descs):

    theta_B = fit_vec[0]
    fit_vec = [theta_B, par_Pi, perp_Pi]
    return chisq_func(fit_vec, phi_B, phi_Pi, res_descs)


def chisq_func(fit_vec, phi_B, phi_Pi, res_descs):

    num_resonance_descs = len(res_descs)
    mag_Bs = [find_mag_B(desc, *fit_vec, phi_B, phi_Pi) for desc in res_descs]

    # find_mag_B_objective returns the sum of squared residuals for a single
    # pair of resonances. We want to sum this over all pairs.
    squared_residuals = [find_mag_B_objective(mag_Bs[ind], res_descs[ind],
                         *fit_vec, phi_B, phi_Pi) for ind
                         in range(num_resonance_descs)]
    sum_squared_residuals = numpy.sum(squared_residuals)

    estimated_st_dev = 0.0001
    estimated_var = numpy.sqrt(estimated_st_dev)
    chisq = sum_squared_residuals / estimated_var

    return chisq


# %% Main


def main(name, res_descs):

    ############ Setup ############

    phi_B = 0.0
    phi_Pi = 0.0
    # phi_B = pi/3
    # phi_Pi = pi/3

    # fit_vec = [theta_B, par_Pi, perp_Pi,]
    param_bounds = ((0, pi/2), (-0.050, 0.050), (0, 0.050))

    res_descs = numpy.array(res_descs)
    for desc in res_descs:
        # Set degenerate resonances to the same value
        if desc[2] is None:
            desc[2] = desc[1]
        # Make sure resonances are sorted
        desc[1:3] = numpy.sort(desc[1:3])

    ############ Guess par_Pi and perp_Pi by zero field ############

    zero_field_res_desc = None
    par_Pi = None
    perp_Pi = None
    # See if we have zero-field resonances
    for desc in res_descs:
        if desc[0] == 0.0:
            zero_field_res_desc = desc
            break

    if zero_field_res_desc is not None:
        
        # Get the splitting and center_freq from the resonances
        zero_field_low = zero_field_res_desc[1]
        zero_field_high = zero_field_res_desc[2]

        zero_field_center = (zero_field_high + zero_field_low) / 2
        par_Pi = zero_field_center - d_gs

        # Similarly
        zero_field_splitting = zero_field_high - zero_field_low
        perp_Pi = zero_field_splitting / 2

    ############ Guess remaining parameters with brute force ############

    if (par_Pi is not None) and (perp_Pi is not None):
        # Just fit for theta_B
        param_ranges = param_bounds[0:1]
        args = (par_Pi, perp_Pi, phi_B, phi_Pi, res_descs)
        x0 = brute(chisq_func_reduced, param_ranges, args=args, Ns=20)
        guess_params = list(x0)
        guess_params.extend([par_Pi, perp_Pi])
    else:
        param_ranges = param_bounds
        args = (phi_B, phi_Pi, res_descs)
        x0 = brute(chisq_func, param_ranges, args=args, Ns=10)
        guess_params = list(x0)

    ############ Fine tuning with minimize ############

    args = (phi_B, phi_Pi, res_descs)
    res = minimize(chisq_func, guess_params, args=args,
                   bounds=param_bounds, method='SLSQP')
    if not res.success:
        print(res.message)
        return

    popt = res.x
    popt_full = numpy.append(popt, [phi_B, phi_Pi])
    
    chisq = res.fun
    print('Chi squared: {:.4g}'.format(chisq))
    degrees_of_freedom = len(res_descs) - len(x0)
    reduced_chisq = res.fun / degrees_of_freedom
    print('Reduced chi squared: {:.4g}'.format(reduced_chisq))

    ############ Plot the result ############

    # Get the mag_B for each pair of resonances with this fit_vec
    mag_Bs = [find_mag_B(desc, *popt_full) for desc in res_descs]

    # Plot the calculated resonances up to the max mag_B
    fig, ax = plot_resonances([0, max(mag_Bs)], *popt_full, name)

    # Plot the resonances
    ax.scatter(mag_Bs, res_descs[:,1])
    ax.scatter(mag_Bs, res_descs[:,2])

    # Return the full 5 parameters
    print('')
    print('theta_B, par_Pi, perp_Pi, phi_B, phi_Pi')
    popt_full_print = [round(el,3) for el in popt_full]
    print(popt_full_print)
    return popt_full


# %% Run the file


if __name__ == '__main__':
    
    # Plots with LaTeX
    plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{physics}',
        r'\usepackage{sfmath}',
        r'\usepackage{upgreek}',
        r'\usepackage{helvet}',
       ]  
    plt.rcParams.update({'font.size': 12})
    # plt.rcParams.update({'font.family': 'sans-serif'})
    # plt.rcParams.update({'font.sans-serif': ['Helvetica']})
    plt.rc('text', usetex=True)
    
    ############ Sample input ############
    
    # Name of the NV, sample, whatever
    name = 'test'
    
    # Enter the resonance descriptions as a list of lists. Each sublist should
    # have the form (all units GHz):
    # [magnetic field if known, lower resonance, higher resonance]
    res_descs = [[0.0, 2.8584, 2.8725],
                 [None, 2.8512, 2.8804],
                 [None, 2.8435, 2.8990],
                 [None, 2.8265, 2.9117],
                 [None, 2.7726, 3.0530],
                 [None, 2.7738, 3.4712],
                 ]
    
    ####################################

    # Run the script
    main(name, res_descs)