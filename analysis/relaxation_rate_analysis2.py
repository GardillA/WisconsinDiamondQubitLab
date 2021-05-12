# -*- coding: utf-8 -*-
"""
Created on Wed Sep 4 14:52:43 2019

This analysis script will plot and evaluate the omega and gamma rates for the
modified rate equations [(0,0) - (0,1) and (1,1) - (1,-1)] for the complete
data set. It calculates a standard error of each data point based on the
statistics over the number of runs. With the standard error on each point, the
subtracted data is then fit to a single exponential. From the (0,0) - (0,1)
exponential, we extact 3*Omega from the exponent, along with the standard
error on omega from the covariance of the fit.

From the (1,1) - (1,-1) exponential, we extract (2*gamma + Omega). Using the
Omega we just found, we calculate gamma and the associated standard error
from the covariance of the fit.

-User can specify if the offset should be a free parameter or if it should be
  set to 0. All our analysis of rates has been done without offset as a free
  param.

-If a value for omega and the omega uncertainty is passed, file will just
  evaluate gamma (t=with the omega provided).


@author: agardill
"""

# %% Imports

import numpy
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

import utils.tool_belt as tool_belt
from utils.tool_belt import States
from figures.relaxation_temp_dependence.revision1.orbach import omega_calc
from figures.relaxation_temp_dependence.revision1.orbach import gamma_calc

# %% Constants

manual_offset_gamma = 0.0

# %% Functions

# The exponential function without an offset
def exp_eq_omega(t, rate, amp):
    return  amp * exp(- rate * t)

def exp_eq_gamma(t, rate, amp):
    return  amp * exp(- rate * t) + manual_offset_gamma #+ 0.01 * exp(-3*0.040*t)

def biexp(t, omega, rate1, amp1, amp2):
    return  amp1 * exp(-rate1*t)

# The exponential function with an offset
def exp_eq_offset(t, rate, amp, offset):
    return  offset + amp * exp(- rate * t)

# A function to collect folders in mass analysis
def get_folder_list(keyword):
    path = 'E:/Shared drives/Kolkowitz Lab Group/nvdata/t1_double_quantum'

    folders = []

    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for folder in d:
                if keyword in folder:
                    folders.append(folder)

    return folders

# This function sorts the data from one folder of an experiment and passes it
# into main
def get_data_lists(folder_name):
    # Get the file list from this folder
    file_list = tool_belt.get_file_list(folder_name, '.txt')

    # Define booleans to be used later in putting data into arrays in the
    # correct order. This was mainly put in place for older data where we
    # took measurements in an inconsistent way (unlike we are now)
    zero_zero_bool = False
    zero_plus_bool = False
    plus_plus_bool = False
    plus_minus_bool = False

    # Initially create empty lists, so that if no data is recieved, a list is
    # still returned from this function
    zero_zero_counts = []
    zero_zero_ste = []
    zero_plus_counts = []
    zero_plus_ste = []
    zero_zero_time = []
    plus_plus_counts = []
    plus_plus_ste = []
    plus_minus_counts = []
    plus_minus_ste = []
    plus_plus_time = []

    # Unpack the data

    num_diffs = 0
    diffs = 0
    splitting_MHz = None

    # Unpack the data and sort into arrays. This allows multiple measurements of
    # the same type to be correctly sorted into one array
    for file in file_list:
        data = tool_belt.get_raw_data(folder_name, file[:-4])
        try:

            init_state_name = data['init_state']
            read_state_name = data['read_state']

            # older files still used 1,-1,0 convention. This will allow old
            # and new files to be evaluated
            if init_state_name == 1 or init_state_name == -1 or  \
                                    init_state_name == 0:
                high_state_name = 1
                low_state_name = -1
                zero_state_name = 0
            else:
                high_state_name = States.HIGH.name
                low_state_name = States.LOW.name
                zero_state_name = States.ZERO.name
            relaxation_time_range = numpy.array(data['relaxation_time_range'])
            num_steps = data['num_steps']

            num_runs = data['num_runs']
            sig_counts  = numpy.array(data['sig_counts'])
            ref_counts = numpy.array(data['ref_counts'])

            # Calculate time arrays in us
            min_relaxation_time, max_relaxation_time = \
                                        relaxation_time_range / 10**6
            time_array = numpy.linspace(min_relaxation_time,
                                        max_relaxation_time, num=num_steps)

            # if (init_state_name != States.ZERO.name) and (max_relaxation_time > 15):
            #     continue

            # Calculate the average signal counts over the runs, and st. error
#            print(sig_counts)
            avg_sig_counts = numpy.average(sig_counts[::], axis=0)
            ste_sig_counts = numpy.std(sig_counts[::], axis=0, ddof = 1) / numpy.sqrt(num_runs)
            avg_ref_counts = numpy.average(ref_counts[::], axis=0)
            ste_ref_counts = numpy.std(ref_counts[::], axis=0, ddof = 1) / numpy.sqrt(num_runs)
            # avg_sig_counts = numpy.average(sig_counts[::])
            # ste_sig_counts = numpy.std(sig_counts[::], ddof = 1) / numpy.sqrt(num_runs*num_steps)

            # Assume reference is constant and can be approximated to one value
            avg_ref = numpy.average(ref_counts[::])
            # avg_ref = 350
            avg_sig = numpy.average(sig_counts[::])
            num_diffs += 1
            diffs += avg_ref - (numpy.std(ref_counts[::], ddof = 1)**2)
            # print(avg_ref)
            # print(numpy.std(ref_counts[::], ddof = 1) / numpy.sqrt(num_runs*num_steps))
            # print(numpy.min(ref_counts))
            # print(avg_sig)
            # print(avg_sig_counts)
            # print()
            # print(avg_sig_counts)
            # if avg_sig > 280:
            #     print(avg_sig_counts)

            # Divide signal by reference to get normalized counts and st error
            # norm_avg_sig = avg_sig_counts / avg_ref_counts
            # norm_avg_sig_ste = ste_sig_counts / avg_ref_counts

            norm_avg_sig = avg_sig_counts / avg_ref
            norm_avg_sig_ste = ste_sig_counts / avg_ref
            # norm_avg_sig = avg_sig_counts
            # norm_avg_sig_ste = ste_sig_counts

            # if 1 not in time_array:
            #     print(init_state_name + ', ' + read_state_name)
            #     sig_33 = ref_counts[:,1]
            #     print(numpy.average(sig_33))
            #     print(numpy.std(sig_33, ddof = 1) / numpy.sqrt(num_runs))
            #     plt.hist(sig_33, 10, (270, 380), histtype='step')
            #     print()

            # time_array = numpy.array(range(0,num_runs*num_steps))
            # norm_avg_sig = sig_counts.flatten() #- avg_ref
            # norm_avg_sig_ste = time_array * 0

            # avg_ref_counts = numpy.average(ref_counts[::], axis=0)
            # ste_ref_counts = numpy.std(ref_counts[::], axis=0, ddof = 1) / numpy.sqrt(num_runs)
            # norm_avg_sig = avg_sig_counts / avg_ref_counts
            # norm_avg_sig_ste = norm_avg_sig * (numpy.sqrt(ste_sig_counts / avg_sig_counts) + numpy.sqrt(ste_ref_counts / avg_ref_counts))**2


            # Check to see which data set the file is for, and append the data
            # to the corresponding array
            if init_state_name == zero_state_name and \
                                read_state_name == zero_state_name:
                # Check to see if data has already been added to a list for
                #this experiment. If it hasn't, then create arrays of the data.
                if zero_zero_bool == False:
                    zero_zero_counts = norm_avg_sig
                    zero_zero_ste = norm_avg_sig_ste
                    zero_zero_time = time_array

                    zero_zero_ref_max_time = max_relaxation_time
                    zero_zero_bool = True
                # If data has already been sorted for this experiment, then check
                # to see if this current data is the shorter or longer measurement,
                # and either append before or after the prexisting data
                else:

                    if max_relaxation_time > zero_zero_ref_max_time:
                        zero_zero_counts = numpy.concatenate((zero_zero_counts,
                                                        norm_avg_sig))
                        zero_zero_ste = numpy.concatenate((zero_zero_ste,
                                                        norm_avg_sig_ste))
                        zero_zero_time = numpy.concatenate((zero_zero_time, time_array))

                    elif max_relaxation_time < zero_zero_ref_max_time:
                        zero_zero_counts = numpy.concatenate((norm_avg_sig,
                                              zero_zero_counts))
                        zero_zero_ste = numpy.concatenate((norm_avg_sig_ste,
                                              zero_zero_ste))
                        zero_zero_time = numpy.concatenate((time_array, zero_zero_time))


            # if init_state_name == zero_state_name and \
            #                     read_state_name == high_state_name:
            # if init_state_name == zero_state_name and \
            #                     read_state_name == low_state_name:
            if (init_state_name == zero_state_name and read_state_name == high_state_name) or \
                (init_state_name == zero_state_name and read_state_name == low_state_name):
                # norm_avg_sig *= 1.01
                if zero_plus_bool == False:
                    zero_plus_counts = norm_avg_sig
                    zero_plus_ste = norm_avg_sig_ste
                    zero_plus_time = time_array

                    zero_plus_ref_max_time = max_relaxation_time
                    zero_plus_bool = True
                else:

                    if max_relaxation_time > zero_plus_ref_max_time:
                        zero_plus_counts = numpy.concatenate((zero_plus_counts,
                                                        norm_avg_sig))
                        zero_plus_ste = numpy.concatenate((zero_plus_ste,
                                                        norm_avg_sig_ste))

                        zero_plus_time = numpy.concatenate((zero_plus_time, time_array))

                    elif max_relaxation_time < zero_plus_ref_max_time:
                        zero_plus_counts = numpy.concatenate((norm_avg_sig,
                                              zero_plus_counts))
                        zero_plus_ste = numpy.concatenate((norm_avg_sig_ste,
                                              zero_plus_ste))

                        zero_plus_time = numpy.concatenate(time_array, zero_plus_time)


            # if (init_state_name == high_state_name) and \
            #     (read_state_name == high_state_name):
            # if (init_state_name == low_state_name) and \
            #     (read_state_name == low_state_name):
            if (init_state_name == high_state_name and read_state_name == high_state_name) or \
                (init_state_name == low_state_name and read_state_name == low_state_name):
                # norm_avg_sig *= 1.0
                if plus_plus_bool == False:
                    plus_plus_counts = norm_avg_sig
                    plus_plus_ste = norm_avg_sig_ste
                    plus_plus_time = time_array

                    plus_plus_ref_max_time = max_relaxation_time
                    plus_plus_bool = True
                else:

                    if max_relaxation_time > plus_plus_ref_max_time:
                        plus_plus_counts = numpy.concatenate((plus_plus_counts,
                                                        norm_avg_sig))
                        plus_plus_ste = numpy.concatenate((plus_plus_ste,
                                                        norm_avg_sig_ste))
                        plus_plus_time = numpy.concatenate((plus_plus_time, time_array))

                    elif max_relaxation_time < plus_plus_ref_max_time:
                        plus_plus_counts = numpy.concatenate((norm_avg_sig,
                                                          plus_plus_counts))
                        plus_plus_ste = numpy.concatenate((norm_avg_sig_ste,
                                                          plus_plus_ste))
                        plus_plus_time = numpy.concatenate((time_array, plus_plus_time))

            # if init_state_name == high_state_name and \
            #                     read_state_name == low_state_name:
            # if init_state_name == low_state_name and \
            #                     read_state_name == high_state_name:
            if (init_state_name == high_state_name and read_state_name == low_state_name) or \
                (init_state_name == low_state_name and read_state_name == high_state_name):
                # We will want to put the MHz splitting in the file metadata
                uwave_freq_init = data['uwave_freq_init']
                uwave_freq_read = data['uwave_freq_read']
                # norm_avg_sig -= 0.0075
                if plus_minus_bool == False:
                    plus_minus_counts = norm_avg_sig
                    plus_minus_ste = norm_avg_sig_ste
                    plus_minus_time = time_array

                    plus_minus_ref_max_time = max_relaxation_time
                    plus_minus_bool = True
                else:

                    if max_relaxation_time > plus_minus_ref_max_time:
                        plus_minus_counts = numpy.concatenate((plus_minus_counts,
                                                        norm_avg_sig))
                        plus_minus_ste = numpy.concatenate((plus_minus_ste,
                                                        norm_avg_sig_ste))
                        plus_minus_time = numpy.concatenate((plus_minus_time, time_array))


                    elif max_relaxation_time < plus_minus_ref_max_time:
                        plus_minus_counts = numpy.concatenate((norm_avg_sig,
                                              plus_minus_counts))
                        plus_minus_ste = numpy.concatenate((norm_avg_sig_ste,
                                              plus_minus_ste))
                        plus_minus_time = numpy.concatenate((time_array, plus_minus_time))


                splitting_MHz = abs(uwave_freq_init - uwave_freq_read) * 10**3

        except Exception:
            print('Skipping {}'.format(str(file)))
            continue

    omega_exp_list = [zero_zero_counts, zero_zero_ste, \
                      zero_plus_counts, zero_plus_ste, \
                      zero_zero_time]
    gamma_exp_list = [plus_plus_counts, plus_plus_ste,  \
                      plus_minus_counts, plus_minus_ste, \
                      plus_plus_time]
    # print(diffs/num_diffs)
    return omega_exp_list, gamma_exp_list, num_runs, splitting_MHz
# %% Main

def main(path, folder, omega = None, omega_ste = None, doPlot = False, offset = True):

    rates_to_zero = False

    path_folder = path + folder
    # Get the file list from the folder
    omega_exp_list, gamma_exp_list, \
                num_runs, splitting_MHz  = get_data_lists(path_folder)

    # %% Fit the data

    if doPlot:
        fig, axes_pack = plt.subplots(1, 2, figsize=(17, 8))
        fig.set_tight_layout(True)

    omega_fit_failed = False
    gamma_fit_failed = False

    ax = None

    # If omega value is passed into the function, skip the omega fitting.
    if omega is not None and omega_ste is not None:
        omega_opti_params = numpy.array([None])
        zero_relaxation_counts = numpy.array([None])
        zero_relaxation_ste = numpy.array([None])
        zero_zero_time = numpy.array([None])
    else:
        #Fit to the (0,0) - (0,1) data to find Omega

        zero_zero_counts = omega_exp_list[0]
        zero_zero_ste = omega_exp_list[1]
        zero_plus_counts = omega_exp_list[2]
        zero_plus_ste = omega_exp_list[3]
        zero_zero_time = omega_exp_list[4]

        zero_relaxation_counts =  zero_zero_counts - zero_plus_counts
        # zero_relaxation_counts = zero_plus_counts
        zero_relaxation_ste = numpy.sqrt(zero_zero_ste**2 + zero_plus_ste**2)

        init_params_list = [0.1, 0.3]

        try:
            if offset:
                init_params_list.append(0)
                init_params = tuple(init_params_list)
                omega_opti_params, cov_arr = curve_fit(exp_eq_offset, zero_zero_time,
                                             zero_relaxation_counts, p0 = init_params,
                                             sigma = zero_relaxation_ste,
                                             absolute_sigma=True)

            else:
                if rates_to_zero:
                    # omega_opti_params = numpy.array([0.045*3,0.175])
                    omega_opti_params = numpy.array([0.0,0.0])
                    cov_arr = numpy.array([[0,0],[0,0]])
                else:
                    init_params = tuple(init_params_list)
                    omega_opti_params, cov_arr = curve_fit(exp_eq_omega, zero_zero_time,
                                                  zero_relaxation_counts, p0 = init_params,
                                                  sigma = zero_relaxation_ste,
                                                  absolute_sigma=True)#,
                                                 # bounds=((0,0),(1,0.5)),
                                                 # loss='soft_l1')

        except Exception:

            omega_fit_failed = True

            if doPlot:
                ax = axes_pack[0]
                ax.errorbar(zero_zero_time, zero_relaxation_counts,
                            yerr = zero_relaxation_ste,
                            label = 'data', fmt = 'o', color = 'blue')
                ax.set_xlabel('Relaxation time (ms)')
                ax.set_ylabel('Normalized signal Counts')
                ax.legend()

        if not omega_fit_failed:
            # Calculate omega nad its ste
            omega = omega_opti_params[0] / 3.0
            omega_ste = numpy.sqrt(cov_arr[0,0]) / 3.0

            print('Omega: {} +/- {} kHz'.format('%.3f'%omega,
                      '%.3f'%omega_ste))
            # Plotting the data
            if doPlot:
                zero_time_linspace = numpy.linspace(0, zero_zero_time[-1], num=1000)
                ax = axes_pack[0]
                ax.errorbar(zero_zero_time, zero_relaxation_counts,
                            yerr = zero_relaxation_ste,
                            label = 'data', fmt = 'o', color = 'blue')
                if offset:
                    ax.plot(zero_time_linspace,
                        exp_eq_offset(zero_time_linspace, *omega_opti_params),
                        'r', label = 'fit')
                else:
                    ax.plot(zero_time_linspace,
                        exp_eq_omega(zero_time_linspace, *omega_opti_params),
                        'r', label = 'fit')
                ax.set_xlabel('Relaxation time (ms)')
                ax.set_ylabel('Normalized signal Counts')
                ax.legend()
                text = r'$\Omega = $ {} $\pm$ {} kHz'.format('%.3f'%omega,
                      '%.3f'%omega_ste)

                props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
                ax.text(0.55, 0.9, text, transform=ax.transAxes, fontsize=12,
                        verticalalignment='top', bbox=props)

    if ax is not None:
        ax.set_title('Omega')
        # ax.set_title('(0,0) - (0,-1)')
        # ax.set_title('(0,0) - (0,+1)')
        # ax.set_yscale('log')

    # %% Fit to the (1,1) - (1,-1) data to find Gamma, only if Omega waas able
    # to fit

    plus_plus_counts = gamma_exp_list[0]
    plus_plus_ste = gamma_exp_list[1]
    plus_minus_counts = gamma_exp_list[2]
    plus_minus_ste = gamma_exp_list[3]
    plus_plus_time = gamma_exp_list[4]



    if len(plus_plus_counts) == 0:
        gamma = None
        gamma_ste = None
        plus_relaxation_counts = numpy.array([])
        plus_relaxation_ste = numpy.array([])
        plus_plus_time = numpy.array([])
        gamma_opti_params = numpy.array([])

    else:

        # Define the counts for the plus relaxation equation
        plus_relaxation_counts =  plus_plus_counts - plus_minus_counts
        # plus_relaxation_counts = plus_minus_counts
        plus_relaxation_ste = numpy.sqrt(plus_plus_ste**2 + plus_minus_ste**2)

        # Skip values at t=0 to get rid of pi pulse decoherence systematic
        # See wiki March 31st, 2021
        inds_to_remove = []
        for ind in range(len(plus_plus_time)):
            t = plus_plus_time[ind]
            if t == 0:
                inds_to_remove.append(ind)
        plus_plus_time = numpy.delete(plus_plus_time, inds_to_remove)
        plus_relaxation_counts = numpy.delete(plus_relaxation_counts, inds_to_remove)
        plus_relaxation_ste = numpy.delete(plus_relaxation_ste, inds_to_remove)
        ax = None
        init_params_list = [2*omega, 0.40]
        try:
            if offset:

                init_params_list.append(0)
                init_params = tuple(init_params_list)
                gamma_opti_params, cov_arr = curve_fit(exp_eq_offset,
                                 plus_plus_time, plus_relaxation_counts,
                                 p0 = init_params, sigma = plus_relaxation_ste,
                                 absolute_sigma=True)


            else:
                if rates_to_zero:
                    gamma_fit_func = lambda t, rate1, amp1, amp2: biexp(t, omega, rate1, amp1, amp2)
                    gamma_opti_params = numpy.array([0.0,0.0,0])
                    cov_arr = numpy.array([[0,0,0],[0,0,0],[0,0,0]])
                else:
                    # MCC
                    init_params = tuple(init_params_list)
                    gamma_fit_func = exp_eq_gamma
                    gamma_opti_params, cov_arr = curve_fit(exp_eq_gamma,
                                      plus_plus_time, plus_relaxation_counts,
                                      p0 = init_params, sigma = plus_relaxation_ste,
                                      absolute_sigma=True)#,
                                      # bounds=((0,0),(1,0.5)),
                                      # loss='soft_l1')
                    # init_params = (0.04, 0.22, 0.17, 0.0)
                    # gamma_fit_func = biexp
                    # init_params = (0.22, 0.17)
                    # gamma_fit_func = lambda t, rate1, amp1: biexp(t, omega, rate1, amp1, -0.005)
                    # init_params = (0.22, 0.17, 0.0)
                    # gamma_fit_func = lambda t, rate1, amp1, amp2: biexp(t, omega, rate1, amp1, amp2)
                    # gamma_opti_params, cov_arr = curve_fit(gamma_fit_func,
                    #                   plus_plus_time, plus_relaxation_counts,
                    #                   p0 = init_params, sigma = plus_relaxation_ste,
                    #                   absolute_sigma=True)
                    # print(gamma_opti_params)

        except Exception as e:
            gamma_fit_failed = True
            print(e)

            if doPlot:
                ax = axes_pack[1]
                ax.errorbar(plus_plus_time, plus_relaxation_counts,
                        yerr = plus_relaxation_ste,
                        label = 'data', fmt = 'o', color = 'blue')
                ax.set_xlabel('Relaxation time (ms)')
                ax.set_ylabel('Normalized signal Counts')

        if not gamma_fit_failed:

            # Calculate gamma and its ste
            gamma = (gamma_opti_params[0] - omega)/ 2.0
            gamma_ste = 0.5 * numpy.sqrt(cov_arr[0,0]+omega_ste**2)

            # Test MCC
            # gamma = 0.070
            # gamma_opti_params[0] = (2 * gamma) + omega
            # gamma_opti_params[1] = 0.20

            print('Gamma: {} +/- {} kHz'.format('%.3f'%gamma,
                      '%.3f'%gamma_ste))

            # Plotting
            if doPlot:
                plus_time_linspace = numpy.linspace(0, plus_plus_time[-1], num=1000)
                ax = axes_pack[1]
                ax.errorbar(plus_plus_time, plus_relaxation_counts,
                        yerr = plus_relaxation_ste,
                        label = 'data', fmt = 'o', color = 'blue')
                if offset:
                    ax.plot(plus_time_linspace,
                        exp_eq_offset(plus_time_linspace, *gamma_opti_params),
                        'r', label = 'fit')
                else:
                    ax.plot(plus_time_linspace,
                        # exp_eq_gamma(plus_time_linspace, *gamma_opti_params),  # MCC
                        gamma_fit_func(plus_time_linspace, *gamma_opti_params),
                        'r', label = 'fit')
                ax.set_xlabel('Relaxation time (ms)')
                ax.set_ylabel('Normalized signal Counts')
                ax.legend()
                text = r'$\gamma = $ {} $\pm$ {} kHz'.format('%.3f'%gamma,
                      '%.3f'%gamma_ste)
    #            ax.set_xlim([-0.001, 0.05])

                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
                ax.text(0.55, 0.90, text, transform=ax.transAxes, fontsize=12,
                        verticalalignment='top', bbox=props)

        if ax is not None:
            ax.set_title('gamma')
            # ax.set_title('(+1,+1) - (+1,-1)')
            # ax.set_title('(-1,-1) - (-1,+1)')
            # ax.set_yscale('log')

    if doPlot:
        fig.canvas.draw()
        fig.canvas.flush_events()

        # Saving the data

        data_dir='E:/Shared drives/Kolkowitz Lab Group/nvdata'

        time_stamp = tool_belt.get_time_stamp()
        raw_data = {'time_stamp': time_stamp,
                    'splitting_MHz': splitting_MHz,
                    'splitting_MHz-units': 'MHz',
#                    'offset_free_param?': offset,
                    'manual_offset_gamma': manual_offset_gamma,
                    'omega': omega,
                    'omega-units': 'kHz',
                    'omega_ste': omega_ste,
                    'omega_ste-units': 'khz',
                    'gamma': gamma,
                    'gamma-units': 'kHz',
                    'gamma_ste': gamma_ste,
                    'gamma_ste-units': 'khz',
                    'zero_relaxation_counts': zero_relaxation_counts.tolist(),
                    'zero_relaxation_counts-units': 'counts',
                    'zero_relaxation_ste': zero_relaxation_ste.tolist(),
                    'zero_relaxation_ste-units': 'counts',
                    'zero_zero_time': zero_zero_time.tolist(),
                    'zero_zero_time-units': 'ms',
                    'plus_relaxation_counts': plus_relaxation_counts.tolist(),
                    'plus_relaxation_counts-units': 'counts',
                    'plus_relaxation_ste': plus_relaxation_ste.tolist(),
                    'plus_relaxation_ste-units': 'counts',
                    'plus_plus_time': plus_plus_time.tolist(),
                    'plus_plus_time-units': 'ms',
                    'omega_opti_params': omega_opti_params.tolist(),
                    'gamma_opti_params': gamma_opti_params.tolist(),
                    }

        file_name = '{}-analysis'.format(folder)
        file_path = '{}/{}/{}'.format(data_dir, path_folder, file_name)
        tool_belt.save_raw_data(raw_data, file_path)
        tool_belt.save_figure(fig, file_path)

        return gamma, gamma_ste
# %% Run the file

if __name__ == '__main__':

    temp = 275

    # path = 'pc_hahn\\branch_cryo-setup\\t1_double_quantum\\data_collections\\'
    # path = 'pc_hahn\\branch_cryo-setup\\t1_dq_knill\\data_collections\\'
    # folder = 'hopper-nv1_2021_03_16-275K-5-gamma_plus_1-long3'.format(temp)

    # est_omega = omega_calc(temp)
    # est_gamma = gamma_calc(temp)
    # print('good times in ms')
    # print('Omega: {}'.format(4000/(3*est_omega)))
    # print('gamma: {}'.format(4000/(2*est_gamma + est_omega)))

    # gamma, ste = main(path, folder, omega=0.0, omega_ste=0.0,
    #                     doPlot=True, offset=False)
    # gamma, ste = main(path, folder, omega=None, omega_ste=None,
    #                   doPlot=True, offset=False)

    # %%

    path = 'pc_hahn\\branch_cryo-setup\\t1_interleave_knill\\data_collections\\trial_data\\'
    folders = [
                # 'hopper-nv1_2021_03_16-275K-50-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-50-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-51-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-51-gamma_plus_1'.format(temp),
                'hopper-nv1_2021_03_16-275K-54-gamma_minus_1'.format(temp),
                'hopper-nv1_2021_03_16-275K-54-gamma_plus_1'.format(temp),
                ]

    for folder in folders:
        # gamma, ste = main(path, folder, omega=0.042, omega_ste=0.0,
        gamma, ste = main(path, folder, omega=None, omega_ste=None,
                          doPlot=True, offset=False)

    path = 'pc_hahn\\branch_cryo-setup\\t1_dq_knill\\data_collections\\trial_data\\'
    folders = [
                # 'hopper-nv1_2021_03_16-275K-3-omega_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-3-omega_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-{}K-4'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-5-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-5-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-7-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-7-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-8-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-8-gamma_plus_1'.format(temp),
                ]

    for folder in folders:
        gamma, ste = main(path, folder, omega=None, omega_ste=None,
                          doPlot=True, offset=False)

    folders = [
                # 'hopper-nv1_2021_03_16-275K-9-gamma_minus_1'.format(temp),  # Just gamma
                # 'hopper-nv1_2021_03_16-275K-9-gamma_plus_1'.format(temp),
                ]

    for folder in folders:
        gamma, ste = main(path, folder, omega=0.040, omega_ste=0.0,
                          doPlot=True, offset=False)

    folders = [
                # 'hopper-nv1_2021_03_16-275K-10-gamma_minus_1'.format(temp),  # No rf, long, zero pulses
                # 'hopper-nv1_2021_03_16-275K-10-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-11-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-11-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-12-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-12-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-13-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-13-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-14-gamma_minus_1'.format(temp),  # Nicest
                # 'hopper-nv1_2021_03_16-275K-14-gamma_plus_1'.format(temp),
                'hopper-nv1_2021_03_16-275K-15-gamma_minus_1'.format(temp),  # Full
                'hopper-nv1_2021_03_16-275K-15-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-16-gamma_minus_1'.format(temp),  # No rf, short, finite pulses
                # 'hopper-nv1_2021_03_16-275K-16-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-17-gamma_minus_1'.format(temp),  # Zero pulses
                # 'hopper-nv1_2021_03_16-275K-17-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-18-omega'.format(temp),  # Zero pulses, just Omega
                # 'hopper-nv1_2021_03_16-275K-19-omega'.format(temp),  # A: Finite pulses, just Omega
                # 'hopper-nv1_2021_03_16-275K-20-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-21-omega'.format(temp),  # B: Finite pulses, just Omega, IQ LOW
                # 'hopper-nv1_2021_03_16-275K-22-omega'.format(temp),  # B
                # 'hopper-nv1_2021_03_16-275K-23-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-24-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-25-omega'.format(temp),  # C: Finite pulses, just Omega, IQ and RF commented out
                # 'hopper-nv1_2021_03_16-275K-26-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-27-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-28-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-29-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-30-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-31-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-32-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-33-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-34-omega'.format(temp),  # C, reverse order
                # 'hopper-nv1_2021_03_16-275K-35-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-36-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-37-omega'.format(temp),  # C
                # 'hopper-nv1_2021_03_16-275K-38-omega'.format(temp),  # A, reverse order
                # 'hopper-nv1_2021_03_16-275K-39-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-40-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-41-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-42-omega'.format(temp),  # A, original order
                # 'hopper-nv1_2021_03_16-275K-43-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-44-omega'.format(temp),  # A
                # 'hopper-nv1_2021_03_16-275K-45-omega'.format(temp),  # A
                ]

    for folder in folders:
        gamma, ste = main(path, folder, omega=None, omega_ste=None,
                          doPlot=True, offset=False)

    # path = 'pc_hahn\\branch_cryo-setup\\t1_double_quantum\\data_collections\\'
    # folders = ['hopper-nv1_2021_03_16-275K-6-gamma_minus_1'.format(temp),
    #             'hopper-nv1_2021_03_16-275K-6-gamma_plus_1'.format(temp),]

    path = 'pc_hahn\\branch_cryo-setup\\t1_double_quantum\\data_collections\\trial_data\\'
    folders = [
                # 'hopper-nv1_2021_03_16-275K-52-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-52-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-53-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-53-gamma_plus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-55-gamma_minus_1'.format(temp),
                # 'hopper-nv1_2021_03_16-275K-55-gamma_plus_1'.format(temp),
                ]

    for folder in folders:
        gamma, ste = main(path, folder, omega=None, omega_ste=None,
                          doPlot=True, offset=False)