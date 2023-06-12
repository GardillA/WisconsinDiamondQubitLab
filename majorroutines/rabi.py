# -*- coding: utf-8 -*-
"""
Rabi flopping routine. Sweeps the pulse duration of a fixed uwave frequency.

Created on Tue Apr 23 11:49:23 2019

@author: mccambria
"""


# %% Imports


import utils.tool_belt as tool_belt
import utils.kplotlib as kpl
from utils.kplotlib import KplColors
import utils.positioning as positioning
import numpy
import os
import time
import matplotlib.pyplot as plt
from random import shuffle
from scipy.optimize import curve_fit
import labrad
import majorroutines.optimize as optimize
from utils.tool_belt import NormStyle



# %% Functions

def fit_data(uwave_time_range, num_steps, fit_func, norm_avg_sig, norm_avg_sig_ste = None):

    #  Set up

    min_uwave_time = uwave_time_range[0]
    max_uwave_time = uwave_time_range[1]
    taus, tau_step = numpy.linspace(min_uwave_time, max_uwave_time,
                            num=num_steps, dtype=numpy.int32, retstep=True)

    # fit_func = tool_belt.cosexp_1_at_0

    #  Estimated fit parameters

    offset = numpy.average(norm_avg_sig)
    decay = 1000

    # To estimate the frequency let's find the highest peak in the FFT
    transform = numpy.fft.rfft(norm_avg_sig)
    freqs = numpy.fft.rfftfreq(num_steps, d=tau_step)
    transform_mag = numpy.absolute(transform)
    # [1:] excludes frequency 0 (DC component)
    max_ind = numpy.argmax(transform_mag[1:])
    frequency = freqs[max_ind + 1]

    #  Fit

    init_params = [offset, frequency, decay]

    try:
        popt, pcov = curve_fit(fit_func, taus, norm_avg_sig,
                            p0=init_params,
                            sigma=norm_avg_sig_ste,
                            absolute_sigma=True)
    except Exception as e:
        print(e)
        popt = None

    return fit_func, popt, pcov

def create_fit_figure(uwave_time_range, num_steps, uwave_freq, norm_avg_sig,
                      norm_avg_sig_ste, fit_func=None, popt=None):

    min_uwave_time = uwave_time_range[0]
    max_uwave_time = uwave_time_range[1]
    taus, tau_step = numpy.linspace(min_uwave_time, max_uwave_time,
                            num=num_steps, dtype=numpy.int32, retstep=True)
    smooth_taus = numpy.linspace(min_uwave_time, max_uwave_time, num=1000)
    
    # Fitting
    if (fit_func is None) or (popt is None):
        fit_func, popt, pcov = fit_data(
            uwave_time_range,
            num_steps,
            fit_func,
            norm_avg_sig,
            norm_avg_sig_ste,
        )
    # Plot setup
    fig, ax = plt.subplots(1,1,figsize=kpl.figsize)
    ax.set_xlabel(r'Microwave duration, $\tau$  (ns)')
    ax.set_ylabel("Normalized fluorescence")
    fig.suptitle('Rabi Measurement')
    
    # Plotting
    if norm_avg_sig_ste is not None:
        kpl.plot_points(ax, taus, norm_avg_sig, yerr=norm_avg_sig_ste)
    else:
        kpl.plot_line(ax, taus, norm_avg_sig)
    kpl.plot_line(
        ax,
        smooth_taus,
        fit_func(smooth_taus, *popt),
        color=KplColors.RED,
    )
    Amp = 1- popt[0]
    decay = abs(popt[2])
    # two_pi = 2 * np.pi
    # amp = 1 - offset
    # return offset + (np.exp(-t / abs(decay)) * abs(amp) * np.cos((two_pi * freq * t)))

    uni_nu = "\u03BD"
    eq_text = r"$(1 - A) + A \cos ( 2 \pi \nu \tau ) e^{-\tau / d}$"
    freq_text = "Set Frequency {} GHz".format(uwave_freq)
    size = kpl.Size.SMALL
    if decay > 2*max_uwave_time:
        base_text = "A = {:.3f} \n1/{} = {:.1f} ns \nd >> {:.0f} ns"
        text = base_text.format(Amp,uni_nu, 1/popt[1], max_uwave_time)
    else:
        base_text = "A = {:.3f} \n1/{} = {:.1f} ns \nd = {:.1f} us"
        text = base_text.format(Amp,uni_nu, 1/popt[1], decay/1e3)
    kpl.anchored_text(ax, text, kpl.Loc.LOWER_LEFT, size=size)
    kpl.anchored_text(ax, freq_text, kpl.Loc.LOWER_RIGHT, size=size)
    kpl.anchored_text(ax, eq_text, kpl.Loc.UPPER_RIGHT, size=size)
    
    return fig, ax, fit_func, popt, pcov

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
    ax_sig_ref.set_ylabel(r"Fluorescence rate (counts / s $\times 10^3$)")
    ax_norm.set_xlabel(r'Microwave duration, $\tau$  (ns)')
    ax_norm.set_ylabel("Normalized fluorescence")
    fig.suptitle('Rabi Measurement')

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


def simulate(uwave_time_range, freq, resonant_freq, contrast,
             measured_rabi_period=None, resonant_rabi_period=None):

    if measured_rabi_period is None:
        resonant_rabi_freq = resonant_rabi_period**-1
        res_dev = freq - resonant_freq
        measured_rabi_freq = numpy.sqrt(res_dev**2 + resonant_rabi_freq**2)
        measured_rabi_period = measured_rabi_freq**-1
        print('measured_rabi_period: {} ns'.format(measured_rabi_period))
    elif resonant_rabi_period is None:
        measured_rabi_freq = measured_rabi_period**-1
        res_dev = freq-resonant_freq
        resonant_rabi_freq = numpy.sqrt(measured_rabi_freq**2 - res_dev**2)
        resonant_rabi_period = resonant_rabi_freq**-1
        print('resonant_rabi_period: {} ns'.format(resonant_rabi_period))
    else:
        raise RuntimeError('Pass either a measured_rabi_period or a ' \
                           'resonant_rabi_period, not both/neither.')

    min_uwave_time = uwave_time_range[0]
    max_uwave_time = uwave_time_range[1]
    smooth_taus = numpy.linspace(min_uwave_time, max_uwave_time,
                          num=1000, dtype=numpy.int32)
    amp = (resonant_rabi_freq / measured_rabi_freq)**2
    angle = measured_rabi_freq * 2 * numpy.pi * smooth_taus / 2
    prob = amp * (numpy.sin(angle))**2

    rel_counts = 1.0 - (contrast * prob)

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    ax.plot(smooth_taus, rel_counts)
    ax.set_xlabel('Tau (ns)')
    ax.set_ylabel('Contrast (arb. units)')

# def simulate_split(uwave_time_range, freq,
#                    res_freq_low, res_freq_high, contrast, rabi_period):

#     rabi_freq = rabi_period**-1

#     min_uwave_time = uwave_time_range[0]
#     max_uwave_time = uwave_time_range[1]
#     smooth_taus = numpy.linspace(min_uwave_time, max_uwave_time,
#                           num=1000, dtype=numpy.int32)

#     omega = numpy.sqrt((freq-res_freq)**2 + rabi_freq**2)
#     amp = (rabi_freq / omega)**2
#     angle = omega * 2 * numpy.pi * smooth_taus / 2
#     prob = amp * (numpy.sin(angle))**2

#     rel_counts = 1.0 - (contrast * prob)

#     fig, ax = plt.subplots(figsize=(8.5, 8.5))
#     ax.plot(smooth_taus, rel_counts)
#     ax.set_xlabel('Tau (ns)')
#     ax.set_ylabel('Contrast (arb. units)')


# %% Main


def main(nv_sig, uwave_time_range, state,
         num_steps, num_reps, num_runs,
         opti_nv_sig = None,
         return_popt=False,close_plot=False):

    with labrad.connect() as cxn:
        rabi_per, sig_counts, ref_counts, popt = main_with_cxn(cxn, nv_sig,
                                         uwave_time_range, state,
                                         num_steps, num_reps, num_runs,
                                         opti_nv_sig,close_plot)

        if return_popt:
            return rabi_per, popt
        if not return_popt:
            return rabi_per


def main_with_cxn(cxn, nv_sig,  uwave_time_range, state,
                  num_steps, num_reps, num_runs,
                  opti_nv_sig = None, close_plot=False):

    counter_server = tool_belt.get_server_counter(cxn)
    pulsegen_server = tool_belt.get_server_pulse_gen(cxn)
    arbwavegen_server = tool_belt.get_server_arb_wave_gen(cxn)

    tool_belt.reset_cfm(cxn)
    kpl.init_kplotlib()

    # %% Get the starting time of the function, to be used to calculate run time

    startFunctionTime = time.time()
    start_timestamp = tool_belt.get_time_stamp()

    # %% Initial calculations and setup
    print(nv_sig['resonance_{}'.format(state.name)])
    uwave_freq = nv_sig['resonance_{}'.format(state.name)]
    uwave_power = nv_sig['uwave_power_{}'.format(state.name)]

    laser_key = 'spin_laser'
    laser_name = nv_sig[laser_key]
    tool_belt.set_filter(cxn, nv_sig, laser_key)
    laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)

    norm_style = nv_sig["norm_style"]
    polarization_time = nv_sig['spin_pol_dur']
    spin_readout_dur = nv_sig['spin_readout_dur']
    readout_sec = spin_readout_dur / (10**9)

    # Array of times to sweep through
    # Must be ints since the pulse streamer only works with int64s
    min_uwave_time = uwave_time_range[0]
    max_uwave_time = uwave_time_range[1]
    taus = numpy.linspace(min_uwave_time, max_uwave_time,
                          num=num_steps, dtype=numpy.int32)

    # check if running external iq_mod with SRS
    iq_key = False
    if 'uwave_iq_{}'.format(state.name) in nv_sig:
        iq_key = nv_sig['uwave_iq_{}'.format(state.name)]

    # Analyze the sequence
    num_reps = int(num_reps)
    # file_name = os.path.basename(__file__)
    seq_args = [taus[0], polarization_time,
                spin_readout_dur, max_uwave_time,
                state.value, laser_name, laser_power]
#    for arg in seq_args:
#        print(type(arg))
    # print(seq_args)
    # return
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    file_name = 'rabi.py'
    ret_vals = pulsegen_server.stream_load(file_name, seq_args_string)
    period = ret_vals[0]
    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.
    sig_counts = numpy.empty([num_runs, num_steps], dtype=numpy.float32)
    sig_counts[:] = numpy.nan
    ref_counts = numpy.copy(sig_counts)
    # norm_avg_sig = numpy.empty([num_runs, num_steps])

    # %% Make some lists and variables to save at the end

    opti_coords_list = []
    tau_index_master_list = [[] for i in range(num_runs)]

    # Create a list of indices to step through the taus. This will be shuffled
    tau_ind_list = list(range(0, num_steps))

    # Create raw data figure for incremental plotting
    raw_fig, ax_sig_ref, ax_norm = create_raw_data_figure(
        taus
    )
    # Set up a run indicator for incremental plotting
    run_indicator_text = "Run #{}/{}"
    text = run_indicator_text.format(0, num_runs)
    run_indicator_obj = kpl.anchored_text(ax_norm, text, loc=kpl.Loc.UPPER_RIGHT)
    
    print('')
    print(tool_belt.get_expected_run_time_string(cxn,'rabi',
                                         period,num_steps,num_reps,num_runs))
    print('')
    # %% Collect the data

    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()
    for run_ind in range(num_runs):

        print('Run index: {}'. format(run_ind))

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

        tool_belt.set_filter(cxn, nv_sig, "spin_laser")
        laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)

        # Apply the microwaves
        sig_gen_cxn = tool_belt.get_server_sig_gen(cxn, state)
        sig_gen_cxn.set_freq(uwave_freq)
        sig_gen_cxn.set_amp(uwave_power)
        if iq_key:
            sig_gen_cxn.load_iq()
            # arbwavegen_server.load_arb_phases([0])
        sig_gen_cxn.uwave_on()

        
        counter_server.start_tag_stream()

        # Shuffle the list of indices to use for stepping through the taus
        shuffle(tau_ind_list)

        # start_time = time.time()
        for tau_ind in tau_ind_list:
            if tool_belt.safe_stop():
                break
            
            if 'daq' in counter_server.name:
                counter_server.load_stream_reader(0, period,  int(2*num_reps))
                n_apd_samples = int(2*num_reps)
#            print(taus[tau_ind])
            # add the tau indexxes used to a list to save at the end
            tau_index_master_list[run_ind].append(tau_ind)
            # Stream the sequence
            seq_args = [taus[tau_ind], polarization_time,
                        spin_readout_dur, max_uwave_time,
                        state.value, laser_name, laser_power]
            seq_args_string = tool_belt.encode_seq_args(seq_args)
            # print(seq_args)
            # Clear the tagger buffer of any excess counts
            counter_server.clear_buffer()

            # start_time = time.time()
            counter_server.clear_buffer()
            pulsegen_server.stream_immediate(file_name, num_reps,
                                             seq_args_string)
            # Get the counts
            new_counts = counter_server.read_counter_separate_gates(n_apd_samples)
            sample_counts = new_counts[0]

            # signal counts are even - get every second element starting from 0
            sig_gate_counts = sample_counts[0::2]
            sig_counts[run_ind, tau_ind] = sum(sig_gate_counts)
#            print('Sig counts: {}'.format(sum(sig_gate_counts)))

            # ref counts are odd - sample_counts every second element starting from 1
            ref_gate_counts = sample_counts[1::2]
            ref_counts[run_ind, tau_ind] = sum(ref_gate_counts)
#            print('Ref counts: {}'.format(sum(ref_gate_counts)))

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


        # %% Save the data we have incrementally for long measurements

        raw_data = {'start_timestamp': start_timestamp,
                    'nv_sig': nv_sig,
                    # 'nv_sig-units': tool_belt.get_nv_sig_units(),
                    'uwave_freq': uwave_freq,
                    'uwave_freq-units': 'GHz',
                    'uwave_power': uwave_power,
                    'uwave_power-units': 'dBm',
                    'uwave_time_range': uwave_time_range,
                    'uwave_time_range-units': 'ns',
                    'state': state.name,
                    'num_steps': num_steps,
                    'num_reps': num_reps,
                    'num_runs': num_runs,
                    'tau_index_master_list':tau_index_master_list,
                    'opti_coords_list': opti_coords_list,
                    'opti_coords_list-units': 'V',
                    'sig_counts': sig_counts.astype(int).tolist(),
                    'sig_counts-units': 'counts',
                    'ref_counts': ref_counts.astype(int).tolist(),
                    'ref_counts-units': 'counts'}

        # This will continuously be the same file path so we will overwrite
        # the existing file with the latest version
        file_path = tool_belt.get_file_path(__file__, start_timestamp,
                                            nv_sig['name'], 'incremental')
        tool_belt.save_raw_data(raw_data, file_path)
        tool_belt.save_figure(raw_fig, file_path)

   
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


    #  Plot the data itself and the fitted curve
    fit_success = True
    try:
        fit_func = tool_belt.cosexp_1_at_0
        fit_fig, ax, fit_func, popt, pcov = create_fit_figure(
            uwave_time_range, num_steps, uwave_freq, norm_avg_sig, norm_avg_sig_ste,
            fit_func 
        )
        rabi_period = 1/popt[1]
        v_unc = numpy.sqrt(pcov[1][1])
        # print(v_unc)
        rabi_period_unc = rabi_period**2 * v_unc
        print('Rabi period measured: {} +/- {} ns\n'.format('%.2f'%rabi_period, '%.2f'%rabi_period_unc))
    except Exception:
        fit_success = False
        print('Unsuccessful fit to data')
        

    # %% Clean up and save the data

    tool_belt.reset_cfm(cxn)

    endFunctionTime = time.time()

    timeElapsed = endFunctionTime - startFunctionTime

    timestamp = tool_belt.get_time_stamp()

    raw_data = {'timestamp': timestamp,
                'timeElapsed': timeElapsed,
                'nv_sig': nv_sig,
                'uwave_freq': uwave_freq,
                'uwave_freq-units': 'GHz',
                'uwave_power': uwave_power,
                'uwave_power-units': 'dBm',
                'uwave_time_range': uwave_time_range,
                'uwave_time_range-units': 'ns',
                'state': state.name,
                'num_steps': num_steps,
                'num_reps': num_reps,
                'num_runs': num_runs,
                'tau_index_master_list':tau_index_master_list,
                'opti_coords_list': opti_coords_list,
                'sig_counts': sig_counts.astype(int).tolist(),
                'sig_counts-units': 'counts',
                'ref_counts': ref_counts.astype(int).tolist(),
                'ref_counts-units': 'counts',
                'norm_avg_sig': norm_avg_sig.astype(float).tolist(),
                'norm_avg_sig-units': 'arb',
                "norm_avg_sig_ste": norm_avg_sig_ste.astype(float).tolist(),
                "norm_avg_sig_ste-units": "arb",}

    nv_name = nv_sig["name"]
    file_path = tool_belt.get_file_path(__file__, timestamp, nv_name)
    tool_belt.save_figure(raw_fig, file_path)
    if fit_success:
        file_path_fit = tool_belt.get_file_path(__file__, timestamp, nv_name + "-fit")
        tool_belt.save_figure(fit_fig, file_path_fit)
    tool_belt.save_raw_data(raw_data, file_path)
    
    tool_belt.save_data_csv(file_path, taus, norm_avg_sig, 'Microwave duration (ns)', 'Normalized fluorescence' )
    
    if close_plot:
        plt.close()

    if fit_success:
        return rabi_period, sig_counts, ref_counts, popt
    else:
        return 0, sig_counts, ref_counts, []



def replot(file):
    kpl.init_kplotlib()

    data = tool_belt.get_raw_data(file)
    
    uwave_freq = data['uwave_freq']
    uwave_time_range = data['uwave_time_range']
    num_steps = data['num_steps']
    min_uwave_time = uwave_time_range[0]
    max_uwave_time = uwave_time_range[1]
    taus = numpy.linspace(min_uwave_time, max_uwave_time, num=num_steps, dtype=numpy.int32)

    sig_counts = data['sig_counts']
    ref_counts = data['ref_counts']
    num_runs = data['num_runs']
    num_reps = data['num_reps']
    nv_sig = data['nv_sig']
    spin_readout_dur = nv_sig['spin_readout_dur']
    norm_style = NormStyle.SINGLE_VALUED
    
    ret_vals = tool_belt.process_counts(sig_counts, ref_counts, num_reps, spin_readout_dur, norm_style)
    (sig_counts_avg_kcps, ref_counts_avg_kcps,
        norm_avg_sig, norm_avg_sig_ste,    ) = ret_vals
    
    #  Plot the data itself and the fitted curve
    fit_func = tool_belt.cosexp_1_at_0
    fit_fig, ax, fit_func, popt, pcov = create_fit_figure( uwave_time_range, num_steps, uwave_freq, norm_avg_sig, norm_avg_sig_ste, fit_func  )
    rabi_period = 1/popt[1]
    v_unc = numpy.sqrt(pcov[1][1])

    rabi_period_unc = rabi_period**2 * v_unc
    print('Rabi period measured: {} +/- {} ns\n'.format('%.2f'%rabi_period, '%.2f'%rabi_period_unc))
    
# %% Run the file


if __name__ == '__main__':
    

    file = '2023_06_01-12_07_43-E6-nv1'
    
    replot(file)
