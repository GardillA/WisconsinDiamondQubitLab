# -*- coding: utf-8 -*-
"""
Spin echo.

Polarize the nv to 0, then applies a pi/2 pulse to send the state to the
equator. Allow to precess for some time, then apply a pi pulse and allow to
precess for the same amount of time, cancelling the previous precession and
resulting in an echo. Finally readout after a second pi/s pulse.

Created on Wed Apr 24 15:01:04 2019

@author: mccambria
"""

# %% Imports


import utils.tool_belt as tool_belt
from scipy.optimize import minimize_scalar
from numpy import pi
import utils.kplotlib as kpl
from utils.kplotlib import KplColors
import numpy
import time
import matplotlib.pyplot as plt
from random import shuffle
import labrad
from utils.tool_belt import States
from scipy.optimize import curve_fit
from scipy.special import j0
from numpy.linalg import eigvals
import majorroutines.optimize as optimize
from utils.tool_belt import NormStyle


# %% Constants


im = 0 + 1j
inv_sqrt_2 = 1 / numpy.sqrt(2)
gmuB = 2.8e-3  # gyromagnetic ratio in GHz / G


def create_raw_data_figure(
    precession_dur_range,
    num_steps,
    avg_sig_counts=None,
    avg_ref_counts=None,
    norm_avg_sig=None,
):
    
    min_precession_dur = precession_dur_range[0]
    max_precession_dur = precession_dur_range[1]
    taus = numpy.linspace(
        min_precession_dur,
        max_precession_dur,
        num=num_steps,
        dtype=numpy.int32,
    )

    T = 2*taus
    
    num_steps = len(taus)
    # Plot setup
    fig, axes_pack = plt.subplots(1, 2, figsize=kpl.figsize_extralarge)
    ax_sig_ref, ax_norm = axes_pack
    ax_sig_ref.set_xlabel(r"$Precession time, T = 2\tau$ ($\mathrm{\mu s}$)")
    ax_norm.set_xlabel(r"Precession time, $T = 2\tau$ ($\mathrm{\mu s}$)")
    ax_sig_ref.set_ylabel(r"Fluorescence rate (counts / s $\times 10^3$)")
    ax_norm.set_ylabel("Normalized fluorescence")
    fig.suptitle("Spin Echo")

    # Plotting
    if avg_sig_counts is None:
        avg_sig_counts = numpy.empty(num_steps)
        avg_sig_counts[:] = numpy.nan
    kpl.plot_line(
        ax_sig_ref, T/1000, avg_sig_counts, label="Signal", color=KplColors.GREEN
    )
    if avg_ref_counts is None:
        avg_ref_counts = numpy.empty(num_steps)
        avg_ref_counts[:] = numpy.nan
    kpl.plot_line(
        ax_sig_ref, T/1000, avg_ref_counts, label="Reference", color=KplColors.RED
    )
    ax_sig_ref.legend(loc=kpl.Loc.LOWER_RIGHT)
    if norm_avg_sig is None:
        norm_avg_sig = numpy.empty(num_steps)
        norm_avg_sig[:] = numpy.nan
    kpl.plot_line(ax_norm, T/1000, norm_avg_sig, color=KplColors.BLUE)
    fig.suptitle('Spin Echo experiment')
    
    return fig, ax_sig_ref, ax_norm

# %% Functionality to calculate the expected magnetic field from the 13C revivals, 
#    and compare to what we expect based on the splitting

def calc_single_hamiltonian(theta_B, center_freq, mag_B):
    # Get parallel and perpendicular components of B field in
    # units of frequency
    par_B = gmuB * mag_B * numpy.cos(theta_B)
    perp_B = gmuB * mag_B * numpy.sin(theta_B)
    hamiltonian = numpy.array(
        [
            [center_freq + par_B, inv_sqrt_2 * perp_B, 0],
            [inv_sqrt_2 * perp_B, 0, inv_sqrt_2 * perp_B],
            [0, inv_sqrt_2 * perp_B, center_freq - par_B],
        ]
    )
    return hamiltonian


def calc_hamiltonian(theta_B, center_freq, mag_B):
    fit_vec = [center_freq, mag_B]
    if (type(theta_B) is list) or (type(theta_B) is numpy.ndarray):
        hamiltonian_list = [
            calc_single_hamiltonian(val, *fit_vec) for val in theta_B
        ]
        return hamiltonian_list
    else:
        return calc_single_hamiltonian(theta_B, *fit_vec)


def calc_res_pair(theta_B, center_freq, mag_B):
    hamiltonian = calc_hamiltonian(theta_B, center_freq, mag_B)
    if (type(theta_B) is list) or (type(theta_B) is numpy.ndarray):
        vals = numpy.sort(eigvals(hamiltonian), axis=1)
        resonance_low = numpy.real(vals[:, 1] - vals[:, 0])
        resonance_high = numpy.real(vals[:, 2] - vals[:, 0])
    else:
        vals = numpy.sort(eigvals(hamiltonian))
        resonance_low = numpy.real(vals[1] - vals[0])
        resonance_high = numpy.real(vals[2] - vals[0])
    return resonance_low, resonance_high


def zfs_cost_func(center_freq, mag_B, theta_B, meas_res_low, meas_res_high):
    calc_res_low, calc_res_high = calc_res_pair(theta_B, center_freq, mag_B)
    diff_low = calc_res_low - meas_res_low
    diff_high = calc_res_high - meas_res_high
    return numpy.sqrt(diff_low ** 2 + diff_high ** 2)


def theta_B_cost_func(
    theta_B, center_freq, mag_B, meas_res_low, meas_res_high
):
    calc_res_low, calc_res_high = calc_res_pair(theta_B, center_freq, mag_B)
    diff_low = calc_res_low - meas_res_low
    diff_high = calc_res_high - meas_res_high
    return numpy.sqrt(diff_low ** 2 + diff_high ** 2)


def plot_resonances_vs_theta_B(data, center_freq=None,revival_time_guess=None,num_revivals_guess=None):

    # %% Setup

    fit_func, popt, stes, fit_fig = fit_data(data,revival_time_guess,num_revivals_guess)
    # print(popt)
    if (fit_func is None) or (popt is None):
        print("Fit failed!")
        return

    nv_sig = data["nv_sig"]
    resonance_LOW = nv_sig["resonance_LOW"]
    resonance_HIGH = nv_sig["resonance_HIGH"]
    # resonance_LOW = 2.7979
    # resonance_HIGH = 2.9456
    # print('test',popt)

    revival_time = popt[1]
    revival_time_ste = stes[1]
    mag_B, mag_B_ste = mag_B_from_revival_time(revival_time, revival_time_ste)

    # %% Angle matching

    # Find the angle that minimizes the distances of the predicted resonances
    # from the measured resonances
    theta_B = None
    if center_freq is None:
        center_freq = (resonance_LOW + resonance_HIGH) / 2
        # print(center_freq)
        # center_freq = 2.8718356422016003
        # print(center_freq)
    args = (center_freq, mag_B, resonance_LOW, resonance_HIGH)
    result = minimize_scalar(
        theta_B_cost_func, bounds=(0, pi / 2), args=args, method="bounded"
    )
    if result.success:
        theta_B = result.x
        theta_B_deg = theta_B * 180 / pi
        print(
            "theta_B = {:.4f} radians, {:.3f} degrees".format(
                theta_B, theta_B_deg
            )
        )
        print("cost = {:.3e}".format(result.fun))
    else:
        print("minimize_scalar failed to find theta_B")

    # %% Plotting

    num_steps = 1000
    linspace_theta_B = numpy.linspace(0, pi / 2, num_steps)

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    fig.set_tight_layout(True)
    res_pairs = calc_res_pair(linspace_theta_B, center_freq, mag_B)
    # res_pairs_high = calc_res_pair(linspace_theta_B, center_freq, mag_B+mag_B_ste)
    # res_pairs_low = calc_res_pair(linspace_theta_B, center_freq, mag_B-mag_B_ste)
    linspace_theta_B_deg = linspace_theta_B * (180 / pi)
    ax.plot(linspace_theta_B_deg, res_pairs[0], label="Calculated low")
    # ax.fill_between(linspace_theta_B_deg, res_pairs_high[0], res_pairs_low[0],
    #                 alpha=0.5)
    ax.plot(linspace_theta_B_deg, res_pairs[1], label="Calculated high")
    # ax.fill_between(linspace_theta_B_deg, res_pairs_high[1], res_pairs_low[1],
    #                 alpha=0.5)

    const = [resonance_LOW for el in range(0, num_steps)]
    ax.plot(linspace_theta_B_deg, const, label="Measured low")
    const = [resonance_HIGH for el in range(0, num_steps)]
    ax.plot(linspace_theta_B_deg, const, label="Measured high")
    

    if theta_B is not None:
        text = r"$\theta_{B} = $%.3f" % (theta_B_deg)
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        ax.text(
            0.05,
            0.65,
            text,
            fontsize=14,
            transform=ax.transAxes,
            verticalalignment="top",
            bbox=props,
        )

    ax.set_xlabel(r"$\theta_{B}$ (deg)")
    ax.set_ylabel("Resonances (GHz)")
    ax.legend()

    fig.canvas.draw()
    fig.set_tight_layout(True)
    fig.canvas.flush_events()

    return fit_func, popt, stes, fit_fig, theta_B_deg, fig


# %% Functions


def mag_B_from_revival_time(revival_time, revival_time_ste=None):
    # 1071 Hz/G is the C13 Larmor precession frequency
    mag_B = ((revival_time / 10 ** 9) * 1071) ** -1
    if revival_time_ste is not None:
        mag_B_ste = mag_B * (revival_time_ste / revival_time)
        return mag_B, mag_B_ste
    else:
        return mag_B


def spin_echo_fit(tau, offset, revival_time, decay_time, T2, n_revivals):
    '''
    Fit from supplement of B. J. Shields et al. "Efficient Readout of a Single Spin State in Diamond 
    via Spin-to-Charge Conversion" PRL (2015)
    '''
    summed_term = 0
    amplitude = 1 - offset
    # print(n_revivals)
    for ind in range(0, int(n_revivals)):
        term = numpy.exp(-(((tau - ind * revival_time) / decay_time) ** 2))
        summed_term += term
    return offset + amplitude * numpy.exp(-(tau/T2)**3) * summed_term


def fit_data(data,revival_time_guess=None,num_revivals_guess=None):

    precession_dur_range = data["precession_time_range"]
    num_steps = data["num_steps"]
    nv_sig = data["nv_sig"]

    # %% Set up

    min_precession_dur = precession_dur_range[0]
    max_precession_dur = precession_dur_range[1]
    taus, tau_step = numpy.linspace(
        min_precession_dur,
        max_precession_dur,
        num=num_steps,
        dtype=numpy.int32,
        retstep=True,
    )

    T = 2*taus
    T_step = 2*tau_step

    fit_func = spin_echo_fit

    # Normalization and uncertainty
    
    try:
        norm_avg_sig = data['norm_avg_sig']
        norm_avg_sig_ste = data['norm_avg_sig_ste']
    except Exception:
        sig_counts = data["sig_counts"]
        ref_counts = data["ref_counts"]
        num_reps = data['num_reps']
        spin_readout_dur = nv_sig['spin_readout_dur']
        norm_style = nv_sig['norm_style']
        
        ret_vals = tool_belt.process_counts(sig_counts, ref_counts, num_reps, spin_readout_dur, norm_style)
        (
            sig_counts_avg_kcps,
            ref_counts_avg_kcps,
            norm_avg_sig,
            norm_avg_sig_ste,
        ) = ret_vals


    # Estimated fit parameters

    # Assume that the bulk of points are the floor and that revivals take
    # us back to 1.0
    offset = numpy.average(norm_avg_sig)
    decay_time = 1e3
    T2 = 100e3

    # To estimate the revival frequency let's find the highest peak in the FFT
    transform = numpy.fft.rfft(norm_avg_sig)
    freqs = numpy.fft.rfftfreq(num_steps, d=T_step)
    transform_mag = numpy.absolute(transform)
    
    fig, ax = plt.subplots()
    ax.plot(freqs[1:]*1e6, transform_mag[1:])
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('FFT Magnitude')

    # For a nice spin echo there may be two dominant frequencies of similar
    # magnitudes. We'll find the right one by brute force
    sorted_inds = numpy.argsort(transform_mag[2:])
    dominant_freqs = [freqs[sorted_inds[-1] + 1], freqs[sorted_inds[-2] + 1]]
    
    
    # Fit

    # The fit doesn't like dealing with vary large numbers. We'll convert to
    # us here and then convert back to ns after the fit for consistency.

    T_us = T / 1000
    decay_time_us = decay_time / 1000
    T2_us = T2 / 1000
    max_precession_dur_us = 2*max_precession_dur / 1000

    # Get the init params we want to test and try them out. Compare them with
    # a scaled chi squared: the sum of squared residuals times the number
    # of degrees of freedom to account for overfitting.
    init_params_tests = []
    min_bounds_tests = []
    max_bounds_tests = []
    best_scaled_chi_sq = None
    best_num_revivals = None
    best_popt = None
    
    for freq in dominant_freqs:

        if revival_time_guess == None:
            revival_time = 1 / freq
        else:
            revival_time = revival_time_guess
            

        if num_revivals_guess == None:
            num_revivals = round(2*max_precession_dur / revival_time) + 1
        else:
            num_revivals = num_revivals_guess
        
        revival_time_us = revival_time / 1000
        init_params = [
            offset,
            revival_time_us,
            decay_time_us,
            T2_us, 
        ]
        min_bounds = (0.5, 0.0, 0.0, 0.0)
        max_bounds = (
            1.0,
            max_precession_dur_us,
            max_precession_dur_us,
            max_precession_dur_us,
        )
        min_bounds_tests.append(min_bounds)
        max_bounds_tests.append(max_bounds)
        init_params_tests.append(init_params)
        
        try:
            fit_func_n = lambda tau, offset, revival_time, decay_time, T2: fit_func(tau, offset, revival_time, 
                                    decay_time, T2, num_revivals)
            # try:
            popt, pcov = curve_fit(
                fit_func_n,
                T_us,
                norm_avg_sig,
                sigma=norm_avg_sig_ste,
                absolute_sigma=True,
                p0=init_params,
                # bounds=(min_bounds, max_bounds),
            )
    
            fit_func_lambda = lambda tau: fit_func_n(tau, *popt)
            residuals = fit_func_lambda(T_us) - norm_avg_sig
            chi_sq = numpy.sum((numpy.array(residuals) ** 2) / (numpy.array(norm_avg_sig_ste) ** 2))
            scaled_chi_sq = chi_sq * len(popt)
            if best_scaled_chi_sq is None or (
                scaled_chi_sq < best_scaled_chi_sq
            ):
                best_scaled_chi_sq = scaled_chi_sq
                best_num_revivals = num_revivals
                best_popt = popt
    
        except Exception as e:
            print(e)
            
    popt = best_popt
    popt = numpy.append(popt, best_num_revivals)
    # print(numpy.append(popt, best_num_revivals))
    popt[1] *= 1000
    popt[2] *= 1000
    popt[3] *= 1000
    print(popt)
    revival_time = popt[1]
    stes = numpy.sqrt(numpy.diag(pcov))
    if (fit_func is not None) and (popt is not None):
        fit_fig = create_fit_figure(
            precession_dur_range,
            num_steps,
            norm_avg_sig,
            norm_avg_sig_ste,
            fit_func,
            popt,
        )
        
    return fit_func, popt, stes, fit_fig


def create_fit_figure(
    precession_dur_range,
    num_steps,
    norm_avg_sig,
    norm_avg_sig_ste,
    fit_func,
    popt,
):

    min_precession_dur = precession_dur_range[0]
    max_precession_dur = precession_dur_range[1]
    taus = numpy.linspace(
        min_precession_dur,
        max_precession_dur,
        num=num_steps,
        dtype=numpy.int32,
    )

    # Account for the pi/2 pulse on each side of a tau
    T = 2*taus

    linspace_T = numpy.linspace(
        min(T), max(T), num=1000)

    fit_fig, ax = plt.subplots(1,1,figsize=kpl.figsize_large)

    kpl.plot_points(ax, T / 1000, norm_avg_sig, color=KplColors.BLUE, label="data")

    kpl.plot_line(ax,
        linspace_T / 1000,
        fit_func(linspace_T, *popt),
        color=KplColors.RED,
        label="fit",
    )
    ax.set_xlabel(r"Precession time, $T = 2\tau$ ($\mathrm{\mu s}$)")
    ax.set_ylabel("Normalized fluorescence")
    fit_fig.suptitle('Spin Echo experiment')
    # ax.legend()

    revival_time = popt[1]
    decay_time = popt[2]
    T_2 = popt[3]
    
    # uni_nu = "\u03BD"
    eq_text = r"$(1 - A) + A e^{-(T / T_2)^3} \sum_0^n e^{-(( T - n T_r) / T_d)^2}$"
    # size = kpl.Size.SMALL
    # if decay > 2*max_uwave_time:
    #     base_text = "A = {:.3f} \n1/{} = {:.1f} ns \nd >> {:.0f} ns"
    #     text = base_text.format(Amp,uni_nu, 1/popt[1], max_uwave_time)
    # else:
    #     base_text = "A = {:.3f} \n1/{} = {:.1f} ns \nd = {:.1f} us"
    #     text = base_text.format(Amp,uni_nu, 1/popt[1], decay/1e3)
    
    
    text_popt = "\n".join(
        (
            r"$A=$%.2f" % (popt[0]),
            r"$T_{d}=$%.3f $\mathrm{\mu s}$" % (decay_time / 1000),
            r"$T_{r}=$%.3f $\mathrm{\mu s}$" % (revival_time / 1000),
            r"$T_{2}=$%.3f $\mathrm{\mu s}$" % (T_2 / 1000),
        ))

    # text_B = "\n".join(
    #     (
    #         "Estimated D.C magnetic field",
    #         r"from $T_r$ is $B=$%.3f G" % (mag_B_from_revival_time(revival_time)),
    #     ))
    
    kpl.anchored_text(ax, eq_text, kpl.Loc.UPPER_RIGHT, size=kpl.Size.SMALL)
    kpl.anchored_text(ax, text_popt, kpl.Loc.LOWER_LEFT, size=kpl.Size.SMALL)
    # kpl.anchored_text(ax, text_B, kpl.Loc.LOWER_RIGHT, size=kpl.Size.SMALL)

    return fit_fig


# %% Main


def main(
    nv_sig,
    precession_dur_range,
    num_steps,
    num_reps,
    num_runs,
    state=States.LOW,
    close_plot=False,
    calc_theta_B = False,
    widqol = False
):

    tool_belt.check_exp_lock()
    tool_belt.set_exp_lock()
    
    with labrad.connect() as cxn:
        angle = main_with_cxn(
            cxn,
            nv_sig,
            precession_dur_range,
            num_steps,
            num_reps,
            num_runs,
            state,
            close_plot,
            calc_theta_B,
            widqol
        )
        return angle


def main_with_cxn(
    cxn,
    nv_sig,
    precession_time_range,
    num_steps,
    num_reps,
    num_runs,
    state=States.LOW,
    close_plot=False,
    calc_theta_B = False,
    widqol = False
):
    
    tool_belt.reset_cfm(cxn)
    kpl.init_kplotlib()
    
    counter_server = tool_belt.get_server_counter(cxn)
    pulsegen_server = tool_belt.get_server_pulse_gen(cxn)
    

    # %% Sequence setup

    laser_key = "spin_laser"
    laser_name = nv_sig[laser_key]
    tool_belt.set_filter(cxn, nv_sig, laser_key)
    laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)
    polarization_time = nv_sig["spin_pol_dur"]
    spin_readout_dur = nv_sig["spin_readout_dur"]
    norm_style = nv_sig['norm_style']

    rabi_period = nv_sig["rabi_{}".format(state.name)]
    uwave_freq = nv_sig["resonance_{}".format(state.name)]
    uwave_power = nv_sig["uwave_power_{}".format(state.name)]

    # Get pulse frequencies
    uwave_pi_pulse = tool_belt.get_pi_pulse_dur(rabi_period)
    uwave_pi_on_2_pulse = tool_belt.get_pi_on_2_pulse_dur(rabi_period)

    seq_file_name = "spin_echo.py"

    # %% Create the array of relaxation times

    # Array of times to sweep through
    # Must be ints
    min_precession_time = int(precession_time_range[0])
    max_precession_time = int(precession_time_range[1])

    taus = numpy.linspace(
        min_precession_time,
        max_precession_time,
        num=num_steps,
        dtype=numpy.int32,
    )
    # print(taus)
    # Account for the pi/2 pulse on each side of a tau
    # plot_taus = (taus + uwave_pi_pulse) / 1000
    # plot_taus = (2*taus) / 1000

    # %% Fix the length of the sequence to account for odd amount of elements

    # Our sequence pairs the longest time with the shortest time, and steps
    # toward the middle. This means we only step through half of the length
    # of the time array.

    # That is a problem if the number of elements is odd. To fix this, we add
    # one to the length of the array. When this number is halfed and turned
    # into an integer, it will step through the middle element.

    if len(taus) % 2 == 0:
        half_length_taus = int(len(taus) / 2)
    elif len(taus) % 2 == 1:
        half_length_taus = int((len(taus) + 1) / 2)

    # Then we must use this half length to calculate the list of integers to be
    # shuffled for each run

    tau_ind_list = list(range(0, half_length_taus))

    # %% Create data structure to save the counts

    # We create an array of NaNs that we'll fill
    # incrementally for the signal and reference counts.
    # NaNs are ignored by matplotlib, which is why they're useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.

    sig_counts = numpy.zeros([num_runs, num_steps])
    sig_counts[:] = numpy.nan
    ref_counts = numpy.copy(sig_counts)

    # %% Make some lists and variables to save at the end

    opti_coords_list = []
    tau_index_master_list = [[] for i in range(num_runs)]

    # %% Analyze the sequence
    
    num_reps = int(num_reps)
    seq_args = [
        min_precession_time,
        polarization_time,
        spin_readout_dur,
        uwave_pi_pulse,
        uwave_pi_on_2_pulse,
        max_precession_time,
        state.value,
        laser_name,
        laser_power,
    ]
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    # print(seq_args)
    ret_vals = pulsegen_server.stream_load(seq_file_name, seq_args_string)
    seq_time = ret_vals[0]
    
    
    # create figure
    raw_fig, ax_sig_ref, ax_norm = create_raw_data_figure(precession_time_range, num_steps)
    # raw_fig, axes_pack = plt.subplots(1, 2, figsize=kpl.figsize_extralarge)
    run_indicator_text = "Run #{}/{}"
    text = run_indicator_text.format(0, num_runs)
    run_indicator_obj = kpl.anchored_text(ax_norm, text, loc=kpl.Loc.UPPER_RIGHT)
    
    print('')
    print(tool_belt.get_expected_run_time_string(cxn,'spin_echo',seq_time,num_steps/2,num_reps,num_runs))
    print('')
    # return
    
    # %% Get the starting time of the function, to be used to calculate run time

    startFunctionTime = time.time()
    start_timestamp = tool_belt.get_time_stamp()

    # %% Collect the data

    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()

    for run_ind in range(num_runs):

        print(" \nRun index: {}".format(run_ind))

        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break

        # Optimize
        opti_coords = optimize.main_with_cxn(cxn, nv_sig)
        opti_coords_list.append(opti_coords)

        # Set up the microwaves
        sig_gen_cxn = tool_belt.get_server_sig_gen(cxn, state)
        sig_gen_cxn.set_freq(uwave_freq)
        sig_gen_cxn.set_amp(uwave_power)
        sig_gen_cxn.uwave_on()
            
        # Set up the laser
        tool_belt.set_filter(cxn, nv_sig, laser_key)
        laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)

        # Load the APD
        counter_server.start_tag_stream()

        # Shuffle the list of tau indices so that it steps thru them randomly
        shuffle(tau_ind_list)

        for tau_ind in tau_ind_list:

            # 'Flip a coin' to determine which tau (long/shrt) is used first
            rand_boolean = numpy.random.randint(0, high=2)

            if rand_boolean == 1:
                tau_ind_first = tau_ind
                tau_ind_second = -tau_ind - 1
            elif rand_boolean == 0:
                tau_ind_first = -tau_ind - 1
                tau_ind_second = tau_ind

            # add the tau indexxes used to a list to save at the end
            tau_index_master_list[run_ind].append(tau_ind_first)
            tau_index_master_list[run_ind].append(tau_ind_second)

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break

            # print(" \nFirst relaxation time: {}".format(taus[tau_ind_first]))
            # print("Second relaxation time: {}".format(taus[tau_ind_second]))

            seq_args = [
                taus[tau_ind_first],
                polarization_time,
                spin_readout_dur,
                uwave_pi_pulse,
                uwave_pi_on_2_pulse,
                taus[tau_ind_second],
                state.value,
                laser_name,
                laser_power,
            ]
            seq_args_string = tool_belt.encode_seq_args(seq_args)
            # Clear the tagger buffer of any excess counts
            counter_server.clear_buffer()
            
            if 'daq' in counter_server.name:
                counter_server.load_stream_reader(0, seq_time,  int(4*num_reps))
                n_apd_samples = int(4*num_reps)
                
            pulsegen_server.stream_immediate(
                seq_file_name, num_reps, seq_args_string
            )
            
            # new_counts = counter_server.read_counter_modulo_gates(4, 1)
            # sample_counts = new_counts[0]
            
            # sig_counts[run_ind, tau_ind_first] = sample_counts[0]
            # ref_counts[run_ind, tau_ind_first] = sample_counts[1]
            # sig_counts[run_ind, tau_ind_second] = sample_counts[2]
            # ref_counts[run_ind, tau_ind_second] = sample_counts[3]
            
            new_counts = counter_server.read_counter_separate_gates(n_apd_samples)
            sample_counts = new_counts[0]

            count = sum(sample_counts[0::4])
            sig_counts[run_ind, tau_ind_first] = count
            # print("First signal = " + str(count))

            count = sum(sample_counts[1::4])
            ref_counts[run_ind, tau_ind_first] = count
            # print("First Reference = " + str(count))

            count = sum(sample_counts[2::4])
            sig_counts[run_ind, tau_ind_second] = count
            # print("Second Signal = " + str(count))

            count = sum(sample_counts[3::4])
            ref_counts[run_ind, tau_ind_second] = count
            # print("Second Reference = " + str(count))

        counter_server.stop_tag_stream()

        # %% incremental plotting
        
        # Update the run indicator
        text = run_indicator_text.format(run_ind + 1, num_runs)
        run_indicator_obj.txt.set_text(text)
        
        # Average the counts over the iterations
        inc_sig_counts = sig_counts[: run_ind + 1]
        inc_ref_counts = ref_counts[: run_ind + 1]
        ret_vals = tool_belt.process_counts(
            inc_sig_counts, inc_ref_counts, num_reps, spin_readout_dur, norm_style
        )
        (
            sig_counts_avg_kcps,
            ref_counts_avg_kcps,
            norm_avg_sig,
            norm_avg_sig_ste,
        ) = ret_vals
        
        kpl.plot_line_update(ax_sig_ref, line_ind=0, y=sig_counts_avg_kcps)
        kpl.plot_line_update(ax_sig_ref, line_ind=1, y=ref_counts_avg_kcps)
        kpl.plot_line_update(ax_norm, y=norm_avg_sig)
        
        
        # %% Save the data we have incrementally for long T1s

        raw_data = {
            "start_timestamp": start_timestamp,
            "nv_sig": nv_sig,
            "nv_sig-units": tool_belt.get_nv_sig_units(cxn),
            "spin_readout_dur": spin_readout_dur,
            "spin_readout_dur-units": "ns",
            "uwave_freq": uwave_freq,
            "uwave_freq-units": "GHz",
            "uwave_power": uwave_power,
            "uwave_power-units": "dBm",
            "rabi_period": rabi_period,
            "rabi_period-units": "ns",
            "uwave_pi_pulse": uwave_pi_pulse,
            "uwave_pi_pulse-units": "ns",
            "uwave_pi_on_2_pulse": uwave_pi_on_2_pulse,
            "uwave_pi_on_2_pulse-units": "ns",
            "precession_time_range": precession_time_range,
            "precession_time_range-units": "ns",
            "state": state.name,
            "num_steps": num_steps,
            "num_reps": num_reps,
            "run_ind": run_ind,
            "taus": taus.tolist(),
            "tau_index_master_list": tau_index_master_list,
            "opti_coords_list": opti_coords_list,
            "opti_coords_list-units": "V",
            "sig_counts": sig_counts.astype(int).tolist(),
            "sig_counts-units": "counts",
            "ref_counts": ref_counts.astype(int).tolist(),
            "ref_counts-units": "counts",
        }

        # This will continuously be the same file path so we will overwrite
        # the existing file with the latest version
        file_path = tool_belt.get_file_path(
            __file__, start_timestamp, nv_sig["name"], "incremental"
        )
        tool_belt.save_raw_data(raw_data, file_path)
        tool_belt.save_figure(raw_fig, file_path)


    # %% Plot the data

    ret_vals = tool_belt.process_counts(sig_counts, ref_counts, num_reps, spin_readout_dur, norm_style)
    (
        sig_counts_avg_kcps,
        ref_counts_avg_kcps,
        norm_avg_sig,
        norm_avg_sig_ste,
    ) = ret_vals
    
    kpl.plot_line_update(ax_sig_ref, line_ind=0, y=sig_counts_avg_kcps)
    kpl.plot_line_update(ax_sig_ref, line_ind=1, y=ref_counts_avg_kcps)
    kpl.plot_line_update(ax_norm, y=norm_avg_sig)
    run_indicator_obj.remove()
    

    # %% Save the data

    endFunctionTime = time.time()

    timeElapsed = endFunctionTime - startFunctionTime

    timestamp = tool_belt.get_time_stamp()

    raw_data = {
        "timestamp": timestamp,
        "timeElapsed": timeElapsed,
        "nv_sig": nv_sig,
        "uwave_freq": uwave_freq,
        "uwave_freq-units": "GHz",
        "uwave_power": uwave_power,
        "uwave_power-units": "dBm",
        "rabi_period": rabi_period,
        "rabi_period-units": "ns",
        "uwave_pi_pulse": uwave_pi_pulse,
        "uwave_pi_pulse-units": "ns",
        "uwave_pi_on_2_pulse": uwave_pi_on_2_pulse,
        "uwave_pi_on_2_pulse-units": "ns",
        "precession_time_range": precession_time_range,
        "precession_time_range-units": "ns",
        "taus": taus.tolist(),
        "state": state.name,
        "num_steps": num_steps,
        "num_reps": num_reps,
        "num_runs": num_runs,
        "tau_index_master_list": tau_index_master_list,
        "opti_coords_list": opti_coords_list,
        "sig_counts": sig_counts.astype(int).tolist(),
        "sig_counts-units": "counts",
        "ref_counts": ref_counts.astype(int).tolist(),
        "ref_counts-units": "counts",
        "norm_avg_sig": norm_avg_sig.astype(float).tolist(),
        "norm_avg_sig-units": "arb",
        "norm_avg_sig_ste": norm_avg_sig_ste.tolist(),
    }

    nv_name = nv_sig["name"]
    file_path = tool_belt.get_file_path(__file__, timestamp, nv_name)
    tool_belt.save_figure(raw_fig, file_path)
    if not widqol:
        tool_belt.save_raw_data(raw_data, file_path)
    tool_belt.save_data_csv(file_path, taus*2/1000, norm_avg_sig, 'Free precession times, T (us)', 'Normalized fluorescence' )

    # %% Fit and save figs
    try:
        ret_vals = fit_data(raw_data, revival_time_guess=60 * 1000)
        fit_func, popt, stes, fit_fig = ret_vals
        theta_B_deg = None
        
        if calc_theta_B:
            ret_vals = plot_resonances_vs_theta_B(raw_data)
            fit_func, popt, stes, fit_fig, theta_B_deg, angle_fig = ret_vals

        file_path_fit = tool_belt.get_file_path(__file__, timestamp, nv_name + "-fit")
        tool_belt.save_figure(fit_fig, file_path_fit)
        if calc_theta_B:
            file_path_angle = tool_belt.get_file_path(__file__, timestamp, nv_name + "-angle")
            tool_belt.save_figure(angle_fig, file_path_angle)
    except Exception:
        print("Fit Failed")
        theta_B_deg = None
        
    if close_plot:
        plt.close()
        
    tool_belt.reset_cfm(cxn)
    tool_belt.set_exp_unlock()
    
    return theta_B_deg


# %% Run the file


if __name__ == "__main__":

    file_name = "2023_02_28-15_38_15-E6-nv1"
    file_name = '2023_01_31-21_37_45-E6-nv1'
    
    data = tool_belt.get_raw_data(file_name)
    nv_name = data['nv_sig']["name"]
    timestamp = data['timestamp']

    revival_time_guess = 60 * 1000
    
    ret_vals = fit_data(data, revival_time_guess=revival_time_guess)
    
    
        
    


