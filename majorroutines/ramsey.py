# -*- coding: utf-8 -*-
"""
Ramsey measruement.

This routine polarizes the nv state into 0, then applies a pi/2 pulse to
put the state into a superposition between the 0 and + or - 1 state. The state
then evolves for a time, tau, of free precesion, and then a second pi/2 pulse
is applied. The amount of population in 0 is read out by collecting the
fluorescence during a readout.

It then takes a fast fourier transform of the time data to attempt to extract
the frequencies in the ramsey experiment. If the funtion can't determine the
peaks in the fft, then a detuning is used.

Lastly, this file curve_fits the data to a triple sum of cosines using the
found frequencies.

Created on Wed Apr 24 15:01:04 2019

@author: agardill
"""

# %% Imports


import utils.tool_belt as tool_belt
import utils.positioning as positioning
import utils.kplotlib as kpl
from utils.kplotlib import KplColors

from scipy.signal import find_peaks
# from numpy import pi
import numpy
import time
import matplotlib.pyplot as plt
from random import shuffle
import labrad
from utils.tool_belt import States
from scipy.optimize import curve_fit
import majorroutines.optimize as optimize

# %% fit

def create_theory_figure(taus_us,offset,decay,amp_1,freq_1, amp_2, freq_2, amp_3, freq_3):
    
    func = tool_belt.cosine_sum
    
    theory_vals = func(taus_us,offset,decay,amp_1,freq_1, amp_2, freq_2, amp_3, freq_3)
    
    fig_fit, ax = plt.subplots(1, 1, figsize=kpl.figsize_large)
    kpl.plot_line(ax,taus_us, theory_vals,color=KplColors.BLUE,label='theory')
    ax.set_xlabel(r'Free precession time ($\mu$s)')
    ax.set_ylabel('Contrast (arb. units)')
    ax.legend()
    text1 = "\n".join((#r'$C + e^{-t/d} [a_1 \mathrm{cos}(2 \pi \nu_1 t) + a_2 \mathrm{cos}(2 \pi \nu_2 t) + a_3 \mathrm{cos}(2 \pi \nu_3 t)]$',
                       r'$C + e^{-t/d\sum_{i=1}^3} a_i\mathrm{cos}(2 \pi \nu_i t)$',
                       r'$C = $' + '%.2f'%(offset),
                        r'$d = $' + '%.2f'%(abs(decay)) + ' us',
                        r'$\nu_1 = $' + '%.2f'%(freq_1) + ' MHz',
                        r'$\nu_2 = $' + '%.2f'%(freq_2) + ' MHz',
                        r'$\nu_3 = $' + '%.2f'%(freq_3) + ' MHz'
                        ))

    size = kpl.Size.SMALL

    kpl.anchored_text(ax, text1, kpl.Loc.LOWER_RIGHT, size=size)


def create_raw_data_figure(
    taus,
    avg_sig_counts=None,
    avg_ref_counts=None,
    norm_avg_sig=None,
):
    num_steps = len(taus)
    # Plot setup
    fig, axes_pack = plt.subplots(1, 2, figsize=kpl.double_figsize)
    ax_sig_ref, ax_norm = axes_pack
    ax_sig_ref.set_xlabel(r'Microwave duration, $\tau$ (ns)')
    ax_sig_ref.set_xlabel(r"Free precession time, $\tau$ ($\mathrm{\mu s}$)")
    ax_sig_ref.set_ylabel(r"Fluorescence rate (counts / s $\times 10^3$)")
    ax_norm.set_xlabel(r"Free precession time, $\tau$ ($\mathrm{\mu s}$)")
    ax_norm.set_ylabel("Normalized fluorescence")
    fig.suptitle('Ramsey experiment')

    # Plotting
    if avg_sig_counts is None:
        avg_sig_counts = numpy.empty(num_steps)
        avg_sig_counts[:] = numpy.nan
    kpl.plot_line(
        ax_sig_ref, taus, avg_sig_counts, label="Signal", color=KplColors.GREEN
    )
    if avg_ref_counts is None:
        avg_ref_counts = numpy.empty(num_steps)
        avg_ref_counts[:] = numpy.nan
    kpl.plot_line(
        ax_sig_ref, taus, avg_ref_counts, label="Reference", color=KplColors.RED
    )
    ax_sig_ref.legend(loc=kpl.Loc.LOWER_RIGHT)
    if norm_avg_sig is None:
        norm_avg_sig = numpy.empty(num_steps)
        norm_avg_sig[:] = numpy.nan
    kpl.plot_line(ax_norm, taus, norm_avg_sig, color=KplColors.BLUE)

    return fig, ax_sig_ref, ax_norm



def extract_oscillations(norm_avg_sig, precession_time_range, num_steps, detuning):
    # Create an empty array for the frequency arrays
    FreqParams = numpy.empty([3])

    # Perform the fft
    time_step = (precession_time_range[1]/1e3 - precession_time_range[0]/1e3) \
                                                    / (num_steps - 1)

    transform = numpy.fft.rfft(norm_avg_sig)
#    window = max_precession_time - min_precession_time
    freqs = numpy.fft.rfftfreq(num_steps, d=time_step)
    transform_mag = numpy.absolute(transform)

    # Plot the fft
    fig_fft, ax= plt.subplots(1, 1, figsize=kpl.figsize_large)
    ax.plot(freqs[1:], transform_mag[1:])  # excludes frequency 0 (DC component)
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('FFT magnitude')
    ax.set_title('Ramsey FFT')
    fig_fft.canvas.draw()
    fig_fft.canvas.flush_events()


    # Guess the peaks in the fft. There are parameters that can be used to make
    # this more efficient
    freq_guesses_ind = find_peaks(transform_mag[1:]
                                  , prominence = 0.5
#                                  , height = 0.8
#                                  , distance = 2.2 / freq_step
                                  )

#    print(freq_guesses_ind[0])

    # Check to see if there are three peaks. If not, try the detuning passed in
    if len(freq_guesses_ind[0]) != 3:
        print('Number of frequencies found: {}'.format(len(freq_guesses_ind[0])))
#        detuning = 3 # MHz

        FreqParams[0] = detuning - 2.2
        FreqParams[1] = detuning
        FreqParams[2] = detuning + 2.2
    else:
        FreqParams[0] = freqs[freq_guesses_ind[0][0]]
        FreqParams[1] = freqs[freq_guesses_ind[0][1]]
        FreqParams[2] = freqs[freq_guesses_ind[0][2]]

    return fig_fft, FreqParams

def fit_ramsey(norm_avg_sig,taus,  precession_time_range, FreqParams):

    taus_us = numpy.array(taus)/1e3
    # Guess the other params for fitting
    amp_1 = -0.04
    amp_2 = amp_1
    amp_3 = amp_1
    decay = 1.85
    offset = .92

    guess_params = (offset, decay, amp_1, FreqParams[0],
                        amp_2, FreqParams[1],
                        amp_3, FreqParams[2])
    
    fit_func = tool_belt.cosine_sum
    init_params = guess_params

    try:
        popt,pcov = curve_fit(fit_func, taus_us, norm_avg_sig,
                      p0=init_params,
                        bounds=([0, 0, -1, -15, -1, -15, -1, -15, ]
                                , [1, numpy.infty, 1, 15, 1, 15, 1, 15, ]) )
    except Exception:
        print('Something went wrong!')
        popt = guess_params

    taus_us_linspace = numpy.linspace(precession_time_range[0]/1e3, precession_time_range[1]/1e3, num=1000)

    fig_fit, ax = plt.subplots(1, 1, figsize=kpl.figsize)
    kpl.plot_line(ax,taus_us, norm_avg_sig,color=KplColors.BLUE,label='data' )
    kpl.plot_line(ax,taus_us_linspace, fit_func(taus_us_linspace,*popt),color=KplColors.RED,label='fit')
    ax.set_xlabel(r'Free precession time, $\tau$ ($\mu$s)')
    ax.set_ylabel("Normalized fluorescence")
    fig_fit.suptitle('Ramsey experiment')
    # ax.legend()
    freqs_sorted = numpy.sort([popt[3],popt[5],popt[7]])
    
    eq_text = r'$C + e^{-\tau/d\sum_{i=1}^3} a_i\mathrm{cos}(2 \pi \nu_i \tau)$'
    

    guess_params = (offset, decay, amp_1, FreqParams[0],
                        amp_2, FreqParams[1],
                        amp_3, FreqParams[2])
    
    text1 = "\n".join((
                       r'$C = $' + '%.2f'%(popt[0]),
                        r'$a_1 = $' + '%.2f'%(popt[2]),
                         r'$a_2 = $' + '%.2f'%(popt[4]),
                          r'$a_3 = $' + '%.2f'%(popt[6]),
                        ))
    
    text2 = "\n".join((
                        r'$d = $' + '%.2f'%(abs(popt[1])) + ' $\mu$s',
                        r'$\nu_1 = $' + '%.2f'%(freqs_sorted[0]) + ' MHz',
                        r'$\nu_2 = $' + '%.2f'%(freqs_sorted[1]) + ' MHz',
                        r'$\nu_3 = $' + '%.2f'%(freqs_sorted[2]) + ' MHz'
                        ))
    size = kpl.Size.SMALL

    kpl.anchored_text(ax, text1, kpl.Loc.LOWER_LEFT, size=size)
    kpl.anchored_text(ax, text2, kpl.Loc.LOWER_RIGHT, size=size)
    kpl.anchored_text(ax, eq_text, kpl.Loc.UPPER_RIGHT, size=size)



#  Plot the data itself and the fitted curve

    fig_fit.canvas.draw()
#    fig.set_tight_layout(True)
    fig_fit.canvas.flush_events()

    return fig_fit

# %% Main


def main(
    nv_sig,
    detuning,
    precession_dur_range,
    num_steps,
    num_reps,
    num_runs,
    state=States.LOW,
    opti_nv_sig = None,
    one_precession_time = False,
    close_plot=False,
    do_fft = True,
    widqol = False
):
    tool_belt.check_exp_lock()
    tool_belt.set_exp_lock()

    with labrad.connect() as cxn:
        angle = main_with_cxn(
            cxn,
            nv_sig,
            detuning,
            precession_dur_range,
            num_steps,
            num_reps,
            num_runs,
            state,
            opti_nv_sig,
            one_precession_time,
            close_plot,
            do_fft,
            widqol
        )
        return angle


def main_with_cxn(
    cxn,
    nv_sig,
    detuning,
    precession_time_range,
    num_steps,
    num_reps,
    num_runs,
    state=States.LOW,
    opti_nv_sig = None,
    one_precession_time = False,
    close_plot=False,
    do_fft = True,
    widqol =False
):

    tool_belt.reset_cfm(cxn)
    
    counter_server = tool_belt.get_server_counter(cxn)
    pulsegen_server = tool_belt.get_server_pulse_gen(cxn)
    # arbwavegen_server = tool_belt.get_server_arb_wave_gen(cxn)


    kpl.init_kplotlib()

    # %% Sequence setup

    laser_key = "spin_laser"
    laser_name = nv_sig[laser_key]
    tool_belt.set_filter(cxn, nv_sig, laser_key)
    laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)
    polarization_time = nv_sig["spin_pol_dur"]
    spin_readout_dur = nv_sig["spin_readout_dur"]
    norm_style = nv_sig["norm_style"]

    rabi_period = nv_sig["rabi_{}".format(state.name)]
    uwave_freq = nv_sig["resonance_{}".format(state.name)]
    uwave_power = nv_sig["uwave_power_{}".format(state.name)]
    # Detune the pi/2 pulse frequency
    uwave_freq_detuned = uwave_freq + detuning / 10**3

    # Get pulse frequencies
    uwave_pi_pulse = 0
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
    # plot_taus = (taus + uwave_pi_pulse) / 1000

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

    # set up to drive transition through zero
    seq_args = [
        min_precession_time/2,
        polarization_time,
        spin_readout_dur,
        uwave_pi_pulse,
        uwave_pi_on_2_pulse,
        max_precession_time/2,
        state.value,
        laser_name,
        laser_power,
    ]
    # print(seq_args)
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    ret_vals = pulsegen_server.stream_load(seq_file_name, seq_args_string)
    seq_time = ret_vals[0]

    # Create raw data figure for incremental plotting
    raw_fig, ax_sig_ref, ax_norm = create_raw_data_figure(
        taus/1000
    )
    # Set up a run indicator for incremental plotting
    run_indicator_text = "Run #{}/{}"
    text = run_indicator_text.format(0, num_runs)
    run_indicator_obj = kpl.anchored_text(ax_norm, text, loc=kpl.Loc.UPPER_RIGHT)
    
    print('')
    print(tool_belt.get_expected_run_time_string(cxn,'ramsey',seq_time,num_steps/2,num_reps,num_runs))
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

        # Optimize and save the coords we found
        if opti_nv_sig:
            opti_coords = optimize.main_with_cxn(cxn, opti_nv_sig)
            drift = positioning.get_drift(cxn)
            adj_coords = nv_sig['coords'] + numpy.array(drift)
            positioning.set_xyz(cxn, adj_coords)
        else:
            opti_coords = optimize.main_with_cxn(cxn, nv_sig)
        opti_coords_list.append(opti_coords)

        # Set up the microwaves
        sig_gen_cxn = tool_belt.get_server_sig_gen(cxn, state)
        sig_gen_cxn.set_freq(uwave_freq_detuned)
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

            if one_precession_time:
                tau_ind_first = 0
                tau_ind_second = 0

            # add the tau indexxes used to a list to save at the end
            tau_index_master_list[run_ind].append(tau_ind_first)
            tau_index_master_list[run_ind].append(tau_ind_second)

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break

            # print(" \nFirst relaxation time: {}".format(taus[tau_ind_first]))
            # print("Second relaxation time: {}".format(taus[tau_ind_second]))

            seq_args = [
                taus[tau_ind_first]/2,
                polarization_time,
                spin_readout_dur,
                uwave_pi_pulse,
                uwave_pi_on_2_pulse,
                taus[tau_ind_second]/2,
                state.value,
                laser_name,
                laser_power,
            ]
            seq_args_string = tool_belt.encode_seq_args(seq_args)
            # Clear the counter/tagger buffer of any excess counts
            counter_server.clear_buffer()
            
            if 'daq' in counter_server.name: #####I PUT THIS BEFORE THE NEXT PARAGRAPH. IT USED TO BE AFTER
                counter_server.load_stream_reader(0, seq_time,  int(4*num_reps))
                n_apd_samples = int(4*num_reps)
            # print(seq_args)
            pulsegen_server.stream_immediate(
                seq_file_name, num_reps, seq_args_string
            )

            # Each sample is of the form [*(<sig_shrt>, <ref_shrt>, <sig_long>, <ref_long>)]
            # So we can sum on the values for similar index modulus 4 to
            # parse the returned list into what we want.
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

        ### Incremental plotting

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
            'detuning': detuning,
            'detuning-units': 'MHz',
            "spin_readout_dur": spin_readout_dur,
            "spin_readout_dur-units": "ns",
            "uwave_freq": uwave_freq_detuned,
            "uwave_freq-units": "GHz",
            "uwave_power": uwave_power,
            "uwave_power-units": "dBm",
            "rabi_period": rabi_period,
            "rabi_period-units": "ns",
            "uwave_pi_on_2_pulse": uwave_pi_on_2_pulse,
            "uwave_pi_on_2_pulse-units": "ns",
            "precession_time_range": precession_time_range,
            "precession_time_range-units": "ns",
            "state": state.name,
            "num_steps": num_steps,
            "num_reps": num_reps,
            "run_ind": run_ind,
            "tau_index_master_list": tau_index_master_list,
            "opti_coords_list": opti_coords_list,
            "opti_coords_list-units": "V",
            "taus": taus.tolist(),
            "taus-units": 'ns',
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

    # %% Hardware clean up

    ### Process and plot the data

    ret_vals = tool_belt.process_counts(sig_counts, ref_counts, num_reps, spin_readout_dur, norm_style)
    (
        sig_counts_avg_kcps,
        ref_counts_avg_kcps,
        norm_avg_sig,
        norm_avg_sig_ste,
    ) = ret_vals

    # Raw data
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
        'detuning': detuning,
        'detuning-units': 'MHz',
        "uwave_freq": uwave_freq_detuned,
        "uwave_freq-units": "GHz",
        "uwave_power": uwave_power,
        "uwave_power-units": "dBm",
        "rabi_period": rabi_period,
        "rabi_period-units": "ns",
        "uwave_pi_on_2_pulse": uwave_pi_on_2_pulse,
        "uwave_pi_on_2_pulse-units": "ns",
        "precession_time_range": precession_time_range,
        "precession_time_range-units": "ns",
        "taus": taus.tolist(),
        "taus-units": 'ns',
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
        "norm_avg_sig_ste": norm_avg_sig_ste.astype(float).tolist(),
        "norm_avg_sig_ste-units": "arb",
    }

    file_path = tool_belt.get_file_path(__file__, timestamp, nv_sig["name"])
    tool_belt.save_figure(raw_fig, file_path)
    if not widqol:
        tool_belt.save_raw_data(raw_data, file_path)
    tool_belt.save_data_csv(file_path, taus/1000, norm_avg_sig, 'Free precession times (us)', 'Normalized fluorescence' )

    # %% Fit and save figs

    # Fourier transform
    fig_fft, FreqParams = extract_oscillations(norm_avg_sig,
                               precession_time_range, num_steps, detuning)
    
    # Save the fft figure
    if do_fft:
        file_path_fft = tool_belt.get_file_path(__file__, timestamp, nv_sig["name"] + '_fft')
        tool_belt.save_figure(fig_fft, file_path_fft)

    # Fit actual data
    fig_fit = fit_ramsey(norm_avg_sig,taus,  precession_time_range, FreqParams)

    # Save the file in the same file directory
    file_path_fit = tool_belt.get_file_path(__file__, timestamp, nv_sig["name"] + '_fit')
    tool_belt.save_figure(fig_fit, file_path_fit)
    
    if close_plot:
        plt.close()

    tool_belt.reset_cfm(cxn)
    tool_belt.set_exp_unlock()
    
    return


# %% Run the file


if __name__ == "__main__":
    
    kpl.init_kplotlib()
    
    file = '2023_01_19-11_14_34-E6test-nv1'

    data = tool_belt.get_raw_data(file)
    detuning= data['detuning']
    nv_sig = data['nv_sig']
    sig_counts = data['sig_counts']
    ref_counts = data['ref_counts']
    norm_avg_sig= numpy.average(sig_counts,axis=0)/numpy.average(ref_counts)
    # norm_avg_sig = data['norm_avg_sig']
    precession_time_range = data['precession_time_range']
    num_steps = data['num_steps']

    taus = data['taus']
    taus = numpy.array(taus)

    analysis= True
    analytics = False
    
    if analysis:
        
        freq_offset = 0.  # our measurement of the resonance frequency could be slightly off because the linewidth is several MHz
        
        freq_1 = detuning + freq_offset
        freq_2 = detuning - 2.2 + freq_offset
        freq_3 = detuning + 2.2 + freq_offset
        
        FreqParams = [freq_1,freq_2,freq_3]
        fit_ramsey(norm_avg_sig,taus,  precession_time_range, FreqParams)

    if analytics:

        func = tool_belt.cosine_sum#(t, offset, decay, amp_1, freq_1, amp_2, freq_2, amp_3, freq_3)
        
        taus_us = taus/1000
        offset=.91
        decay = 7
        amp_1 = -.04
        amp_2 = amp_1
        amp_3 = amp_1
        detuning = 4
        freq_offset = -1
        freq_1 = detuning-2.2 + freq_offset
        freq_2 = detuning + freq_offset
        freq_3 = detuning+2.2 + freq_offset
        
        create_theory_figure(taus_us, offset, decay, amp_1, freq_1, amp_2, freq_2, amp_3, freq_3)
