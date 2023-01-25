# -*- coding: utf-8 -*-
"""
Electron spin resonance routine. Scans the microwave frequency, taking counts
at each point.

Created on Thu Apr 11 15:39:23 2019

@author: mccambria
"""


# %% Imports


import utils.positioning as positioning
import utils.tool_belt as tool_belt
import numpy as np
import matplotlib.pyplot as plt
import labrad
from utils.tool_belt import States, NormStyle
from majorroutines import pulsed_resonance 
from random import shuffle
import majorroutines.optimize as optimize
import utils.kplotlib as kpl
from utils.kplotlib import KplColors
import time

def process_counts(ref_counts, sig_counts, norm_style=NormStyle.SINGLE_VALUED):
    """Extract the normalized average signal at each data point.
    Since we sometimes don't do many runs (<10), we often will have an
    insufficient sample size to run stats on for norm_avg_sig calculation.
    We assume Poisson statistics instead.
    """

    ref_counts = np.array(ref_counts)
    sig_counts = np.array(sig_counts)

    num_runs, num_points = ref_counts.shape

    # Find the averages across runs
    sig_counts_avg = np.average(sig_counts, axis=0)
    single_ref_avg = np.average(ref_counts)
    ref_counts_avg = np.average(ref_counts, axis=0)

    sig_counts_ste = np.sqrt(sig_counts_avg) / np.sqrt(num_runs)
    single_ref_ste = np.sqrt(single_ref_avg) / np.sqrt(num_runs * num_points)
    ref_counts_ste = np.sqrt(ref_counts_avg) / np.sqrt(num_runs)

    if norm_style == NormStyle.SINGLE_VALUED:
        norm_avg_sig = sig_counts_avg / single_ref_avg
        norm_avg_sig_ste = norm_avg_sig * np.sqrt(
            (sig_counts_ste / sig_counts_avg) ** 2
            + (single_ref_ste / single_ref_avg) ** 2
        )
    elif norm_style == NormStyle.POINT_TO_POINT:
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

# %% Main


def main(nv_sig, freq_center, freq_range,
         num_steps, num_runs, uwave_power, state=States.LOW, opti_nv_sig = None, close_plot=False):

    with labrad.connect() as cxn:
        main_with_cxn(cxn, nv_sig,  freq_center, freq_range,
                      num_steps, num_runs, uwave_power, state, opti_nv_sig, close_plot)

def main_with_cxn(cxn, nv_sig,  freq_center, freq_range,
                  num_steps, num_runs, uwave_power, state=States.LOW, opti_nv_sig = None, close_plot=False):

    kpl.init_kplotlib()
    
    # %% Initial calculations and setup

    tool_belt.reset_cfm(cxn)
    
    counter_server = tool_belt.get_server_counter(cxn)
    pulsegen_server = tool_belt.get_server_pulse_gen(cxn)
    
    # Set up the laser
    laser_key = 'spin_laser'
    laser_name = nv_sig[laser_key]
    laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)

    # Since this is CW we need the imaging readout rather than the spin 
    # readout typically used for state detection
    readout = nv_sig['imaging_readout_dur']  
    readout_sec = readout / (10**9)
    norm_style = nv_sig["norm_style"]
    
    file_name = 'resonance.py'
    seq_args = [readout, state.value, laser_name, laser_power, ]
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    # print(seq_args)
    # return

    # Calculate the frequencies we need to set
    half_freq_range = freq_range / 2
    freq_low = freq_center - half_freq_range
    freq_high = freq_center + half_freq_range
    freqs = np.linspace(freq_low, freq_high, num_steps)
    freq_ind_list = list(range(num_steps))
    freq_ind_master_list = []

    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    # counts = np.empty(num_steps)
    # counts[:] = np.nan

    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.
    ref_counts = np.empty([num_runs, num_steps])
    ref_counts[:] = np.nan
    sig_counts = np.copy(ref_counts)

    opti_coords_list = []
    
    ret_vals = pulsegen_server.stream_load(file_name, seq_args_string)
    period = ret_vals[0]
    
    print('')
    print(tool_belt.get_expected_run_time_string(cxn,'resonance',period,num_steps,1,num_runs))
    print('')
    
    # Create raw data figure for incremental plotting
    raw_fig, ax_sig_ref, ax_norm = pulsed_resonance.create_raw_data_figure(
        freq_center, freq_range, num_steps
    )
    # Set up a run indicator for incremental plotting
    run_indicator_text = "Run #{}/{}"
    text = run_indicator_text.format(0, num_runs)
    run_indicator_obj = kpl.anchored_text(ax_norm, text, loc=kpl.Loc.UPPER_RIGHT)

    # %% Get the starting time of the function

    start_timestamp = tool_belt.get_time_stamp()

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
            adj_coords = nv_sig['coords'] + np.array(drift)
            positioning.set_xyz(cxn, adj_coords)
        else:
            opti_coords = optimize.main_with_cxn(cxn, nv_sig)
        opti_coords_list.append(opti_coords)
        
        # Laser setup
        tool_belt.set_filter(cxn, nv_sig, laser_key)
        laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)
        # Start the laser now to get rid of transient effects
        # tool_belt.laser_on(cxn, laser_name, laser_power)
    
        sig_gen_cxn = tool_belt.get_server_sig_gen(cxn, state)
        sig_gen_cxn.set_amp(uwave_power)
        sig_gen_cxn.uwave_on()

        # Load the APD task with two samples for each frequency step
        ret_vals = pulsegen_server.stream_load(file_name, seq_args_string)
        period = ret_vals[0]
        
        if 'daq' in counter_server.name:
            counter_server.load_stream_reader(0, period,  2*num_steps)#put the total number of samples you expect for this run
        else:
            counter_server.start_tag_stream()
        
        # Shuffle the list of frequency indices so that we step through
        # them randomly
        shuffle(freq_ind_list)
        freq_ind_master_list.append(freq_ind_list)

        # Take a sample and increment the frequency
        for step_ind in range(num_steps):

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break

            freq_ind = freq_ind_list[step_ind]
            sig_gen_cxn.set_freq(freqs[freq_ind])

            # Start the timing stream
            counter_server.clear_buffer()
            pulsegen_server.stream_start() 

            new_counts = counter_server.read_counter_separate_gates(2) #originally 1
            sample_counts = new_counts[0]
            ref_gate_counts = sample_counts[0::2]
            ref_counts[run_ind, freq_ind]  = sum(ref_gate_counts)

            sig_gate_counts = sample_counts[1::2]
            sig_counts[run_ind, freq_ind] = sum(sig_gate_counts)

        counter_server.stop_tag_stream()
        
        ### Incremental plotting

        # Update the run indicator
        text = run_indicator_text.format(run_ind + 1, num_runs)
        run_indicator_obj.txt.set_text(text)

        # Average the counts over the iterations
        inc_sig_counts = sig_counts[: run_ind + 1]
        inc_ref_counts = ref_counts[: run_ind + 1]
        ret_vals = tool_belt.process_counts(
            inc_sig_counts, inc_ref_counts, 1, readout, norm_style
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

        rawData = {'start_timestamp': start_timestamp,
                   'nv_sig': nv_sig,
                   # 'nv_sig-units': tool_belt.get_nv_sig_units(),
                   'opti_coords_list': opti_coords_list,
                   'opti_coords_list-units': 'V',
                   'freq_center': freq_center,
                   'freq_center-units': 'GHz',
                   'freq_range': freq_range,
                   'freq_range-units': 'GHz',
                   'num_steps': num_steps,
                   'num_runs': num_runs,
                   'freq_ind_master_list': freq_ind_master_list,
                   'uwave_power': uwave_power,
                   'uwave_power-units': 'dBm',
                   'readout': readout,
                   'readout-units': 'ns',
                   'sig_counts': sig_counts.astype(int).tolist(),
                   'sig_counts-units': 'counts',
                   'ref_counts': ref_counts.astype(int).tolist(),
                   'ref_counts-units': 'counts'}

        # This will continuously be the same file path so we will overwrite
        # the existing file with the latest version
        file_path = tool_belt.get_file_path(__file__, start_timestamp,
                                            nv_sig['name'], 'incremental')
        tool_belt.save_raw_data(rawData, file_path)
        
        tool_belt.save_figure(raw_fig, file_path)

    # %% Process and plot the data
    
    ret_vals = tool_belt.process_counts(
        sig_counts, ref_counts, 1, readout, norm_style
    )
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
    
    # Fits
    fit_fig, _, fit_func, popt, _ = pulsed_resonance.create_fit_figure(
        freq_center, freq_range, num_steps, norm_avg_sig, norm_avg_sig_ste
    )


    if len(popt) == 3:
        low_freq = round(popt[2],4)
        high_freq = None
        print('Single resonance: ',low_freq,'GHz') 
    elif len(popt) == 6:
        low_freq = round(popt[2],4)
        high_freq = round(popt[5],4)
        print('Low resonance: ',low_freq,'GHz') 
        print('High resonance: ',high_freq,'GHz')
        print('Slitting = ',(high_freq-low_freq)*1000,'MHz')

    # %% Clean up and save the data

    tool_belt.reset_cfm(cxn)

    timestamp = tool_belt.get_time_stamp()

    data = {'timestamp': timestamp,
               'nv_sig': nv_sig,
               # 'nv_sig-units': tool_belt.get_nv_sig_units(),
               'opti_coords_list': opti_coords_list,
               'opti_coords_list-units': 'V',
               'freq_center': freq_center,
               'freq_center-units': 'GHz',
               'freq_range': freq_range,
               'freq_range-units': 'GHz',
               'num_steps': num_steps,
               'num_runs': num_runs,
               'freq_ind_master_list': freq_ind_master_list,
               'uwave_power': uwave_power,
               'uwave_power-units': 'dBm',
               'readout': readout,
               'readout-units': 'ns',
               'sig_counts': sig_counts.astype(int).tolist(),
               'sig_counts-units': 'counts',
               'ref_counts': ref_counts.astype(int).tolist(),
               'ref_counts-units': 'counts',
               'norm_avg_sig': norm_avg_sig.astype(float).tolist(),
               'norm_avg_sig-units': 'arb',
#               'norm_avg_sig_ste': norm_avg_sig_ste.astype(float).tolist(),
#               'norm_avg_sig_ste-units': 'arb',
               }

    nv_name = nv_sig['name']
    
    file_path = tool_belt.get_file_path(__file__, timestamp, nv_name)
    data_file_name = file_path.stem
    tool_belt.save_figure(raw_fig, file_path)

    tool_belt.save_raw_data(data, file_path)

    file_path = tool_belt.get_file_path(__file__, timestamp, nv_name + "-fit")
    tool_belt.save_figure(fit_fig, file_path)
    
    if close_plot:
        plt.close()

    return None, None


def replot(file):

    data = tool_belt.get_raw_data(file)

    freq_center = data['freq_center']
    freq_range = data['freq_range']
    num_steps = data['num_steps']
    num_runs = data['num_runs']
    ref_counts = data['ref_counts']
    sig_counts = data['sig_counts']
    readout = data['readout']
    norm_style = NormStyle.SINGLE_VALUED #data['nv_sig']['norm_style']

    ret_vals = tool_belt.process_counts(
        sig_counts, ref_counts, 1, readout, norm_style
    )
    (
        sig_counts_avg_kcps,
        ref_counts_avg_kcps,
        norm_avg_sig,
        norm_avg_sig_ste,
    ) = ret_vals

    # Raw data
    kpl.init_kplotlib()
    
    raw_fig, ax_sig_ref, ax_norm = pulsed_resonance.create_raw_data_figure(
        freq_center, freq_range, num_steps
    )
    kpl.plot_line_update(ax_sig_ref, line_ind=0, y=sig_counts_avg_kcps)
    kpl.plot_line_update(ax_sig_ref, line_ind=1, y=ref_counts_avg_kcps)
    kpl.plot_line_update(ax_norm, y=norm_avg_sig)
    
    # Fits
    fit_fig, _, fit_func, popt, _ = pulsed_resonance.create_fit_figure(
        freq_center, freq_range, num_steps, norm_avg_sig, norm_avg_sig_ste
    )
# %%

if __name__ == '__main__':

    file = '2023_01_14-12_59_48-E6test-nv1'
    
    replot(file)
