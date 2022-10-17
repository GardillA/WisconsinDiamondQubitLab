# -*- coding: utf-8 -*-
"""
Electron spin resonance routine. Scans the microwave frequency, taking counts
at each point.

Created on Thu Apr 11 15:39:23 2019

@author: mccambria
"""

# %% Imports


import utils.tool_belt as tool_belt
import majorroutines.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import labrad
from utils.tool_belt import States
from random import shuffle
import sys
import utils.kplotlib as kpl


# %% Figure functions


def create_fit_figure(
    freq_range,
    freq_center,
    num_steps,
    norm_avg_sig,
    fit_func,
    popt,
    norm_avg_sig_ste=None,
):

    freqs = calculate_freqs(freq_range, freq_center, num_steps)
    smooth_freqs = calculate_freqs(freq_range, freq_center, 1000)

    fig, ax = plt.subplots()
    kpl.plot_line(ax, freqs, norm_avg_sig)
    kpl.plot_line(ax, smooth_freqs, fit_func(smooth_freqs, *popt))
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Normalized fluorescence")
    # ax.legend(loc="lower right")

    text = "\n".join(
        (
            "A = {:.3f}",
            r"$\sigma$ = {:.4f} GHz",
            "f = {:.4f} GHz",
        )
    )
    if fit_func == single_gaussian_dip:
        low_text = text.format(*popt[0:3])
        high_text = None
    elif fit_func == double_gaussian_dip:
        low_text = text.format(*popt[0:3])
        high_text = text.format(*popt[3:6])

    kpl.text(ax, 0.05, 0.05, low_text)
    if high_text is not None:
        kpl.text(ax, 0.74, 0.05, high_text)

    kpl.tight_layout(fig)
    fig.canvas.draw()
    fig.canvas.flush_events()

    return fig


# %%


def return_res_with_error(data):
    """
    Returns the frequency/error of the deepest resonance in a spectrum.
    Intended for extracting the frequency/error of a single resonance.
    """

    freq_center = data["freq_center"]
    freq_range = data["freq_range"]
    num_steps = data["num_steps"]
    ref_counts = data["ref_counts"]
    sig_counts = data["sig_counts"]
    num_runs = data["num_runs"]
    ret_vals = process_counts(ref_counts, sig_counts, num_runs)
    (
        avg_ref_counts,
        avg_sig_counts,
        norm_avg_sig,
        ste_ref_counts,
        ste_sig_counts,
        norm_avg_sig_ste,
    ) = ret_vals

    fit_func, popt, pcov = fit_resonance(
        freq_range, freq_center, num_steps, norm_avg_sig, norm_avg_sig_ste
    )
    if len(popt) == 6:
        print("Double resonance detected!")
        low_res_depth = popt[0]
        high_res_depth = popt[3]
        if low_res_depth > high_res_depth:
            res_ind = 2
        else:
            res_ind = 5
    else:
        res_ind = 2
    res = popt[res_ind]
    res_err = np.sqrt(pcov[res_ind, res_ind])
    return res, res_err


def calculate_freqs(freq_range, freq_center, num_steps):
    half_freq_range = freq_range / 2
    freq_low = freq_center - half_freq_range
    freq_high = freq_center + half_freq_range
    return np.linspace(freq_low, freq_high, num_steps)


def gaussian(freq, constrast, sigma, center):
    return constrast * np.exp(-((freq - center) ** 2) / (2 * (sigma ** 2)))


def double_gaussian_dip(
    freq,
    low_constrast,
    low_sigma,
    low_center,
    high_constrast,
    high_sigma,
    high_center,
):
    low_gauss = gaussian(freq, low_constrast, low_sigma, low_center)
    high_gauss = gaussian(freq, high_constrast, high_sigma, high_center)
    return 1.0 - low_gauss - high_gauss


def single_gaussian_dip(freq, constrast, sigma, center):
    return 1.0 - gaussian(freq, constrast, sigma, center)


# def get_guess_params(freqs, norm_avg_sig, ref_counts):
def get_guess_params(
    freq_range, freq_center, num_steps, norm_avg_sig, ref_counts=None
):

    # %% Guess the locations of the minimums

    freqs = calculate_freqs(freq_range, freq_center, num_steps)

    contrast = 0.15  # Arb
    sigma = 0.003  # GHz
    #    sigma = 0.010  # MHz
    fwhm = 2.355 * sigma

    # Convert to index space
    fwhm_ind = fwhm * (num_steps / freq_range)
    if fwhm_ind < 1:
        fwhm_ind = 1

    # Bit of processing
    inverted_norm_avg_sig = 1 - norm_avg_sig
    if ref_counts is not None:
        # ref_counts contains a list of lists. Each list is a single run.
        # Each point is a single freq in that run. We want to know the
        # noise we should expect to see on a point averaged over the runs.
        ref_ste = np.std(ref_counts) / np.sqrt(len(ref_counts))
        rel_ref_ste = ref_ste / np.average(ref_counts)
        height = 5 * rel_ref_ste
        # print(height)
    else:
        height = 0.06

    # Peaks must be separated from each other by the estimated fwhm (rayleigh
    # criteria), have a contrast of at least the noise or 5% (whichever is
    # greater), and have a width of at least two points
    peak_inds, details = find_peaks(
        inverted_norm_avg_sig, distance=fwhm_ind, height=height, width=2
    )
    peak_inds = peak_inds.tolist()
    peak_heights = details["peak_heights"].tolist()

    low_freq_guess = None
    high_freq_guess = None
    if len(peak_inds) > 1:
        # Find the location of the highest peak
        max_peak_height = max(peak_heights)
        max_peak_peak_inds = peak_heights.index(max_peak_height)
        max_peak_freqs = peak_inds[max_peak_peak_inds]

        # Remove what we just found so we can find the second highest peak
        peak_inds.pop(max_peak_peak_inds)
        peak_heights.pop(max_peak_peak_inds)

        # Find the location of the next highest peak
        next_max_peak_height = max(peak_heights)
        next_max_peak_peak_inds = peak_heights.index(
            next_max_peak_height
        )  # Index in peak_inds
        next_max_peak_freqs = peak_inds[
            next_max_peak_peak_inds
        ]  # Index in freqs

        # List of higest peak then next highest peak
        peaks = [max_peak_freqs, next_max_peak_freqs]

        # Only keep the smaller peak if it's at least 50% the height of the larger peak
        if next_max_peak_height > 0.5 * max_peak_height:
            # Sort by frequency
            peaks.sort()
            low_freq_guess = freqs[peaks[0]]
            high_freq_guess = freqs[peaks[1]]
        else:
            low_freq_guess = freqs[peaks[0]]
            high_freq_guess = None

    elif len(peak_inds) == 1:
        low_freq_guess = freqs[peak_inds[0]]
        high_freq_guess = None
    else:
        print("Could not locate peaks, using center frequency")
        low_freq_guess = freq_center
        high_freq_guess = None

    low_freq_guess = 2.861
    high_freq_guess = 2.882

    if low_freq_guess is None:
        return None, None

    # %% Fit!

    if high_freq_guess is None:
        fit_func = single_gaussian_dip
        guess_params = [contrast, sigma, low_freq_guess]
    else:
        fit_func = double_gaussian_dip
        guess_params = [
            contrast,
            sigma,
            low_freq_guess,
            contrast,
            sigma,
            high_freq_guess,
        ]

    return fit_func, guess_params


def fit_resonance(
    freq_range,
    freq_center,
    num_steps,
    norm_avg_sig,
    norm_avg_sig_ste=None,
    ref_counts=None,
):

    freqs = calculate_freqs(freq_range, freq_center, num_steps)

    fit_func, guess_params = get_guess_params(
        freq_range, freq_center, num_steps, norm_avg_sig, ref_counts
    )

    # fit_func = single_gaussian_dip
    # guess_params = [0.2, 0.004, freq_center]

    # try:
    if norm_avg_sig_ste is not None:
        popt, pcov = curve_fit(
            fit_func,
            freqs,
            norm_avg_sig,
            p0=guess_params,
            sigma=norm_avg_sig_ste,
            absolute_sigma=True,
        )
        # popt = guess_params
        # if len(popt) == 6:
        #     zfs = (popt[2] + popt[5]) / 2
        #     low_res_err = np.sqrt(pcov[2,2])
        #     hig_res_err = np.sqrt(pcov[5,5])
        #     zfs_err = np.sqrt(low_res_err**2 + hig_res_err**2) / 2
        # else:
        #     zfs = popt[2]
        #     zfs_err = np.sqrt(pcov[2,2])

        # print(zfs)
        # print(zfs_err)
        # temp_from_resonances.main(zfs, zfs_err)

    else:
        if fit_func == single_gaussian_dip:
            fit_bounds = (0, np.infty)
        else:
            fit_bounds = (0, np.infty)
        print(guess_params)
        popt, pcov = curve_fit(
            fit_func, freqs, norm_avg_sig, p0=guess_params, bounds=fit_bounds
        )
    # except Exception as e:
    #     print(e)
    #     popt = guess_params
    #     pcov = None

    return fit_func, popt, pcov


def simulate(res_freq, freq_range, contrast, rabi_period, uwave_pulse_dur):

    rabi_freq = rabi_period ** -1

    smooth_freqs = calculate_freqs(freq_range, res_freq, 1000)

    omega = np.sqrt((smooth_freqs - res_freq) ** 2 + rabi_freq ** 2)
    amp = (rabi_freq / omega) ** 2
    angle = (
        omega * 2 * np.pi * uwave_pulse_dur / 2
    )  # we use frequencies, so we have to convert by 2 pi here
    prob = amp * (np.sin(angle)) ** 2

    rel_counts = 1.0 + (contrast * prob)

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    ax.plot(smooth_freqs, rel_counts)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Contrast (arb. units)")

    return smooth_freqs, rel_counts


def process_counts(ref_counts, sig_counts, num_runs):
    """
    Extract the normalized average signal at each data point.
    Assume the reference measurements are independent of the signal
    measurements (as they should be) so we can use a single averaged
    value for the reference. Since we sometimes don't do many runs
    (<10), we often will have an insufficient sample size to run
    stats on for norm_avg_sig calculation. We assume Poisson
    statistics instead.
    """

    # Find the averages across runs
    sig_counts_avg = np.average(sig_counts, axis=0)
    single_ref_avg = np.average(ref_counts)
    ref_counts_avg = np.average(ref_counts, axis=0)

    sig_counts_ste = np.sqrt(sig_counts_avg) / np.sqrt(num_runs)
    single_ref_ste = np.sqrt(single_ref_avg) / np.sqrt(num_runs)
    ref_counts_ste = np.sqrt(ref_counts_avg) / np.sqrt(num_runs)

    # New style, single reference
    if True:
        norm_avg_sig = sig_counts_avg / single_ref_avg
        norm_avg_sig_ste = norm_avg_sig * np.sqrt(
            (sig_counts_ste / sig_counts_avg) ** 2
            + (single_ref_ste / single_ref_avg) ** 2
        )

    # Old style, point-by-point reference
    else:
        norm_avg_sig = sig_counts_avg / ref_counts_avg
        norm_avg_sig_ste = norm_avg_sig * np.sqrt(
            (sig_counts_ste / sig_counts_avg) ** 2
            + (ref_counts_ste / ref_counts_avg) ** 2
        )

    return (
        ref_counts_avg,
        sig_counts_avg,
        norm_avg_sig,
        ref_counts_ste,
        sig_counts_ste,
        norm_avg_sig_ste,
    )


# %% User functions


def state(
    nv_sig,
    apd_indices,
    state,
    freq_range,
    num_steps,
    num_reps,
    num_runs,
    composite=False,
    opti_nv_sig=None,
    ret_file_name=False,
):

    freq_center = nv_sig["resonance_{}".format(state.name)]
    uwave_power = nv_sig["uwave_power_{}".format(state.name)]
    uwave_pulse_dur = nv_sig["rabi_{}".format(state.name)] // 2

    resonance_list = main(
        nv_sig,
        apd_indices,
        freq_center,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
        uwave_power,
        uwave_pulse_dur,
        state,
        composite,
        opti_nv_sig,
        ret_file_name,
    )

    return resonance_list
    # return resonance_list, nv_sig


# %% Main


def main(
    nv_sig,
    apd_indices,
    freq_center,
    freq_range,
    num_steps,
    num_reps,
    num_runs,
    uwave_power,
    uwave_pulse_dur,
    state=States.HIGH,
    composite=False,
    opti_nv_sig=None,
    ret_file_name=False,
):

    with labrad.connect() as cxn:
        resonance_list = main_with_cxn(
            cxn,
            nv_sig,
            apd_indices,
            freq_center,
            freq_range,
            num_steps,
            num_reps,
            num_runs,
            uwave_power,
            uwave_pulse_dur,
            state,
            composite,
            opti_nv_sig,
            ret_file_name,
        )
    return resonance_list


def main_with_cxn(
    cxn,
    nv_sig,
    apd_indices,
    freq_center,
    freq_range,
    num_steps,
    num_reps,
    num_runs,
    uwave_power,
    uwave_pulse_dur,
    state=States.HIGH,
    composite=False,
    opti_nv_sig=None,
    ret_file_name=False,
):

    # %% Initial calculations and setup
    
    kpl.init_kplotlib(font_size="small")

    tool_belt.reset_cfm(cxn)

    # Calculate the frequencies we need to set
    half_freq_range = freq_range / 2
    freq_low = freq_center - half_freq_range
    freq_high = freq_center + half_freq_range
    freqs = np.linspace(freq_low, freq_high, num_steps)

    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.
    ref_counts = np.empty([num_runs, num_steps])
    ref_counts[:] = np.nan
    sig_counts = np.copy(ref_counts)

    laser_key = "spin_laser"
    laser_name = nv_sig[laser_key]
    laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)

    polarization_time = nv_sig["spin_pol_dur"]
    readout = nv_sig["spin_readout_dur"]
    readout_sec = readout / (10 ** 9)
    if composite:
        uwave_pi_pulse = round(nv_sig["rabi_{}".format(state.name)] / 2)
        uwave_pi_on_2_pulse = round(nv_sig["rabi_{}".format(state.name)] / 4)
        seq_args = [
            polarization_time,
            readout,
            uwave_pi_pulse,
            uwave_pi_on_2_pulse,
            1,
            1,
            apd_indices[0],
            state.value,
            laser_name,
            laser_power,
        ]
        seq_args = [int(el) for el in seq_args]
    else:
        seq_args = [
            uwave_pulse_dur,
            polarization_time,
            readout,
            uwave_pulse_dur,
            apd_indices[0],
            state.value,
            laser_name,
            laser_power,
        ]
    # print(seq_args)
    # return
    seq_args_string = tool_belt.encode_seq_args(seq_args)

    opti_coords_list = []

    # create figure
    fig, axes_pack = plt.subplots(1, 2, figsize=(17, 8.5))
    ax = axes_pack[0]
    ax.plot([], [])
    ax.set_title("Non-normalized Count Rate Versus Frequency")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Count rate (kcps)")

    ax = axes_pack[1]
    ax.plot([], [])
    ax.set_title("Normalized Count Rate vs Frequency")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Contrast (arb. units)")

    # %% Get the starting time of the function

    start_timestamp = tool_belt.get_time_stamp()

    # %% Collect the data

    # Create a list of indices to step through the freqs. This will be shuffled
    freq_index_master_list = [[] for i in range(num_runs)]
    freq_ind_list = list(range(0, num_steps))

    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()

    for run_ind in range(num_runs):
        print("Run index: {}".format(run_ind))

        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break

        # Optimize and save the coords we found
        if opti_nv_sig:
            opti_coords = optimize.main_with_cxn(cxn, opti_nv_sig, apd_indices)
            drift = tool_belt.get_drift()
            adj_coords = nv_sig["coords"] + np.array(drift)
            tool_belt.set_xyz(cxn, adj_coords)
        else:
            opti_coords = optimize.main_with_cxn(cxn, nv_sig, apd_indices)
        opti_coords_list.append(opti_coords)

        # Set up the microwaves and laser. Then load the pulse streamer
        # (must happen after optimize and iq_switch since run their
        # own sequences)
        sig_gen_cxn = tool_belt.get_signal_generator_cxn(cxn, state)
        sig_gen_cxn.set_amp(uwave_power)
        if composite:
            sig_gen_cxn.load_iq()
            cxn.arbitrary_waveform_generator.load_knill()
        tool_belt.set_filter(cxn, nv_sig, laser_key)
        laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)
        if composite:
            ret_vals = cxn.pulse_streamer.stream_load(
                "discrete_rabi2.py", seq_args_string
            )
        else:
            ret_vals = cxn.pulse_streamer.stream_load(
                "rabi.py", seq_args_string
            )

        # Start the tagger stream
        cxn.apd_tagger.start_tag_stream(apd_indices)

        # Take a sample and step through the shuffled frequencies
        shuffle(freq_ind_list)
        for freq_ind in freq_ind_list:

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break

            freq_index_master_list[run_ind].append(freq_ind)
            sig_gen_cxn.set_freq(freqs[freq_ind])
            sig_gen_cxn.uwave_on()

            # It takes 400 us from receipt of the command to
            # switch frequencies so allow 1 ms total
            #            time.sleep(0.001)
            # Clear the tagger buffer of any excess counts
            cxn.apd_tagger.clear_buffer()
            # Start the timing stream
            cxn.pulse_streamer.stream_start(int(num_reps))

            # Get the counts
            new_counts = cxn.apd_tagger.read_counter_separate_gates(1)

            sample_counts = new_counts[0]

            # signal counts are even - get every second element starting from 0
            sig_gate_counts = sample_counts[0::2]
            sig_counts[run_ind, freq_ind] = sum(sig_gate_counts)
            # print(sum(sig_gate_counts))

            # ref counts are odd - sample_counts every second element starting from 1
            ref_gate_counts = sample_counts[1::2]
            ref_counts[run_ind, freq_ind] = sum(ref_gate_counts)
            # print(sum(ref_gate_counts))

        cxn.apd_tagger.stop_tag_stream()

        # %% incremental plotting

        # Average the counts over the iterations
        single_ref_avg = np.average(ref_counts[: (run_ind + 1)])
        ref_counts_avg = np.average(ref_counts[: (run_ind + 1)], axis=0)
        sig_counts_avg = np.average(sig_counts[: (run_ind + 1)], axis=0)
        norm_avg_sig = sig_counts_avg / single_ref_avg

        kcps_uwave_off_avg = (ref_counts_avg / (num_reps * 1000)) / readout_sec
        kcpsc_uwave_on_avg = (sig_counts_avg / (num_reps * 1000)) / readout_sec

        ax = axes_pack[0]
        ax.cla()
        ax.plot(freqs, kcps_uwave_off_avg, "r-", label="Reference")
        ax.plot(freqs, kcpsc_uwave_on_avg, "g-", label="Signal")

        ax.legend()

        ax = axes_pack[1]
        ax.cla()
        ax.plot(freqs, norm_avg_sig, "b-")

        text_popt = f"Run # {run_ind + 1}/{num_runs}"
        kpl.text(ax, 0.8, 0.9, text_popt)

        kpl.tight_layout(fig)
        fig.canvas.draw()
        fig.canvas.flush_events()

        # %% Save the data we have incrementally for long measurements

        rawData = {
            "start_timestamp": start_timestamp,
            "nv_sig": nv_sig,
            "nv_sig-units": tool_belt.get_nv_sig_units(),
            "freq_center": freq_center,
            "freq_center-units": "GHz",
            "freq_range": freq_range,
            "freq_range-units": "GHz",
            "uwave_pulse_dur": uwave_pulse_dur,
            "uwave_pulse_dur-units": "ns",
            "state": state.name,
            "num_steps": num_steps,
            "run_ind": run_ind,
            "uwave_power": uwave_power,
            "uwave_power-units": "dBm",
            "readout": readout,
            "readout-units": "ns",
            "freq_index_master_list": freq_index_master_list,
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
        tool_belt.save_raw_data(rawData, file_path)
        tool_belt.save_figure(fig, file_path)

    # %% Process and plot the data

    ret_vals = process_counts(ref_counts, sig_counts, num_runs)
    (
        avg_ref_counts,
        avg_sig_counts,
        norm_avg_sig,
        ste_ref_counts,
        ste_sig_counts,
        norm_avg_sig_ste,
    ) = ret_vals

    # Convert to kilocounts per second
    kcps_uwave_off_avg = (avg_ref_counts / (num_reps * 1000)) / readout_sec
    kcpsc_uwave_on_avg = (avg_sig_counts / (num_reps * 1000)) / readout_sec

    # The first plot will display both the uwave_off and uwave_off counts
    ax = axes_pack[0]
    ax.cla()
    ax.plot(freqs, kcps_uwave_off_avg, "r-", label="Reference")
    ax.plot(freqs, kcpsc_uwave_on_avg, "g-", label="Signal")
    ax.set_title("Non-normalized Count Rate Versus Frequency")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Count rate (kcps)")
    ax.legend()
    # The second plot will show their subtracted values
    ax = axes_pack[1]
    ax.cla()
    ax.plot(freqs, norm_avg_sig, "b-")
    ax.set_title("Normalized Count Rate vs Frequency")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Contrast (arb. units)")

    fig.canvas.draw()
    fig.tight_layout()
    fig.canvas.flush_events()

    # %% Clean up and save the data

    tool_belt.reset_cfm(cxn)

    timestamp = tool_belt.get_time_stamp()

    rawData = {
        "timestamp": timestamp,
        "nv_sig": nv_sig,
        "nv_sig-units": tool_belt.get_nv_sig_units(),
        "opti_coords_list": opti_coords_list,
        "opti_coords_list-units": "V",
        "freq_center": freq_center,
        "freq_center-units": "GHz",
        "freq_range": freq_range,
        "freq_range-units": "GHz",
        "uwave_pulse_dur": uwave_pulse_dur,
        "uwave_pulse_dur-units": "ns",
        "state": state.name,
        "num_steps": num_steps,
        "num_reps": num_reps,
        "num_runs": num_runs,
        "uwave_power": uwave_power,
        "uwave_power-units": "dBm",
        "readout": readout,
        "readout-units": "ns",
        "freq_index_master_list": freq_index_master_list,
        "opti_coords_list": opti_coords_list,
        "opti_coords_list-units": "V",
        "sig_counts": sig_counts.astype(int).tolist(),
        "sig_counts-units": "counts",
        "ref_counts": ref_counts.astype(int).tolist(),
        "ref_counts-units": "counts",
        "norm_avg_sig": norm_avg_sig.astype(float).tolist(),
        "norm_avg_sig-units": "arb",
        "norm_avg_sig_ste": norm_avg_sig_ste.astype(float).tolist(),
        "norm_avg_sig_ste-units": "arb",
    }

    nv_name = nv_sig["name"]
    filePath = tool_belt.get_file_path(__file__, timestamp, nv_name)
    # print(filePath)
    if ret_file_name:
        raw_file_name = filePath.stem()
    tool_belt.save_figure(fig, filePath)
    tool_belt.save_raw_data(rawData, filePath)

    # %% Fit the data

    fit_func, popt, pcov = fit_resonance(
        freq_range, freq_center, num_steps, norm_avg_sig, norm_avg_sig_ste
    )
    if (fit_func is not None) and (popt is not None):
        fit_fig = create_fit_figure(
            freq_range, freq_center, num_steps, norm_avg_sig, fit_func, popt
        )
    else:
        fit_fig = None

    filePath = tool_belt.get_file_path(__file__, timestamp, nv_name + "-fit")
    if fit_fig is not None:
        tool_belt.save_figure(fit_fig, filePath)

    # %% Return

    if fit_func == single_gaussian_dip:
        print("Single resonance at {:.4f} GHz".format(popt[2]))
        print("\n")
        ret_vals = [popt[2], None]
    elif fit_func == double_gaussian_dip:
        print(
            "Resonances at {:.4f} GHz and {:.4f} GHz".format(popt[2], popt[5])
        )
        print("Splitting of {:d} MHz".format(int((popt[5] - popt[2]) * 1000)))
        print("\n")
        ret_vals = [popt[2], popt[5]]
    else:
        print("No resonances found")
        print("\n")
        ret_vals = [None, None]

    if ret_file_name:
        ret_vals.append(raw_file_name)
    # print(ret_vals)
    return ret_vals


# %% Run the file


if __name__ == "__main__":

    # print(__file__)
    # sys.exit()

    # file = "2022_10_13-16_05_00-hopper-nv4_2022_10_13"

    # data = tool_belt.get_raw_data(file)

    # freq_center = data["freq_center"]
    # freq_range = data["freq_range"]
    # num_steps = data["num_steps"]
    # num_runs = data["num_runs"]
    # norm_avg_sig = np.array(data["norm_avg_sig"])
    # norm_avg_sig_ste = np.array(data["norm_avg_sig_ste"])

    # freqs = calculate_freqs(freq_range, freq_center, num_steps)

    # # ax.plot(freqs, norm_avg_sig, label=label_list[f])
    # # ax.set_xlabel("Frequency (GHz)")
    # # ax.set_ylabel("Contrast (arb. units)")
    # # ax.legend(loc="lower right")

    # fit_func, popt, pcov = fit_resonance(
    #     freq_range, freq_center, num_steps, norm_avg_sig, norm_avg_sig_ste
    # )
    # create_fit_figure(
    #     freq_range,
    #     freq_center,
    #     num_steps,
    #     norm_avg_sig,
    #     fit_func,
    #     popt,
    # )
    # print(popt)
    # print(pcov)

    kpl.init_kplotlib(font_size="small")
    # matplotlib.rcParams["axes.linewidth"] = 1.0

    file = "2022_10_13-18_23_57-hopper-nv4_2022_10_13"
    data = tool_belt.get_raw_data(file)
    freq_center = data["freq_center"]
    freq_range = data["freq_range"]
    num_steps = data["num_steps"]
    ref_counts = data["ref_counts"]
    sig_counts = data["sig_counts"]
    num_runs = data["num_runs"]
    ret_vals = process_counts(ref_counts, sig_counts, num_runs)
    (
        avg_ref_counts,
        avg_sig_counts,
        norm_avg_sig,
        ste_ref_counts,
        ste_sig_counts,
        norm_avg_sig_ste,
    ) = ret_vals

    fit_func, popt, pcov = fit_resonance(
        freq_range,
        freq_center,
        num_steps,
        norm_avg_sig,
        norm_avg_sig_ste,
        ref_counts,
    )

    # popt[2] -= np.sqrt(pcov[2, 2])

    create_fit_figure(
        freq_range,
        freq_center,
        num_steps,
        norm_avg_sig,
        fit_func,
        popt,
        norm_avg_sig_ste=norm_avg_sig_ste,
    )

    # plt.show(block=True)

    # res_freq, freq_range, contrast, rabi_period, uwave_pulse_dur
    # simulate(2.8351, 0.035, 0.02, 170, 170/2)
