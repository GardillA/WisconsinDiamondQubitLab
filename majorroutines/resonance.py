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
import numpy
import os
import matplotlib.pyplot as plt


# %% Main


def main(cxn, coords, nd_filter, apd_index, freq_center, freq_range,
         num_steps, num_runs, uwave_power, name='untitled'):

    # %% Get the starting time of the function
    
    timestampStart = tool_belt.get_time_stamp()
    
    # %% Initial calculations and setup
    
    # Set up for the pulser - we can't load the sequence yet until after 
    # optimize runs since optimize loads its own sequence
    readout = 100 * 10**6  # 0.1 s
    readout_sec = readout / (10**9)
    uwave_switch_delay = 100 * 10**6  # 0.1 s to open the gate
    sequence_args = [readout, uwave_switch_delay, apd_index]

    file_name = os.path.basename(__file__)
    file_name_no_ext = os.path.splitext(file_name)[0]

    # Calculate the frequencies we need to set
    half_freq_range = freq_range / 2
    freq_low = freq_center - half_freq_range
    freq_high = freq_center + half_freq_range
    freqs = numpy.linspace(freq_low, freq_high, num_steps)

    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    counts = numpy.empty(num_steps)
    counts[:] = numpy.nan

    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.
    ref_counts = numpy.empty([num_runs, num_steps], dtype=numpy.uint32)
    ref_counts[:] = numpy.nan
    sig_counts = numpy.copy(ref_counts)
    counts_norm = numpy.empty([num_runs, num_steps])

    # %% Collect the data

#    tool_belt.set_xyz(cxn, coords)

    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()

    for run_ind in range(num_runs):
        print(run_ind)

        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break

        optimize.main(cxn, coords, apd_index)

        # Load the APD task with two samples for each frequency step
        ret_vals = cxn.pulse_streamer.stream_load(file_name, sequence_args)
        period = ret_vals[0]
        cxn.apd_counter.load_stream_reader(apd_index, period, 2 * num_steps)

        # Take a sample and increment the frequency
        for step_ind in range(num_steps):

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break

            cxn.microwave_signal_generator.set_freq(freqs[step_ind])

            # If this is the first sample then we have to enable the signal
            if (run_ind == 0) and (step_ind == 0):
                cxn.microwave_signal_generator.set_amp(uwave_power)
                cxn.microwave_signal_generator.uwave_on()

            # Start the timing stream
            cxn.pulse_streamer.stream_start()

            new_counts = cxn.apd_counter.read_stream(apd_index, 2)
            if len(new_counts) != 2:
                raise RuntimeError('There should be exactly 2 samples per freq.')

            ref_counts[run_ind, step_ind] = new_counts[0]
            sig_counts[run_ind, step_ind] = new_counts[1]
            try:
                counts_norm[run_ind, step_ind] = new_counts[1] / new_counts[0]
            except Exception:
                counts_norm[run_ind, step_ind] = 0

    # %% Process and plot the data

    # Find the averages across runs
    ref_counts_avg = numpy.average(ref_counts, axis=0)
    sig_counts_avg = numpy.average(sig_counts, axis=0)
    norm_avg_sig = numpy.average(counts_norm, axis=0)

    # Convert to kilocounts per second
    kcps_uwave_off_avg = (ref_counts_avg / (10**3)) / readout_sec
    kcpsc_uwave_on_avg = (sig_counts_avg / (10**3)) / readout_sec

    # Create an image with 2 plots on one row, with a specified size
    # Then draw the canvas and flush all the previous plots from the canvas
    fig, axes_pack = plt.subplots(1, 2, figsize=(17, 8.5))

    # The first plot will display both the uwave_off and uwave_off counts
    ax = axes_pack[0]
    ax.plot(freqs, kcps_uwave_off_avg, 'r-', label = 'Signal')
    ax.plot(freqs, kcpsc_uwave_on_avg, 'g-', label = 'Reference')
    ax.set_title('Non-normalized Count Rate Versus Frequency')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Count rate (kcps)')
    ax.legend()
    # The second plot will show their subtracted values
    ax = axes_pack[1]
    ax.plot(freqs, norm_avg_sig, 'b-')
    ax.set_title('Normalized Count Rate vs Frequency')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Contrast (arb. units)')

    fig.canvas.draw()
    fig.tight_layout()
    fig.canvas.flush_events()

    # %% Turn off the RF and save the data

    cxn.microwave_signal_generator.uwave_off()

    timestampEnd = tool_belt.get_time_stamp()

    rawData = {'timestanpStart': timestampStart,
               'timestampEnd': timestampEnd,
               'name': name,
               'coords': coords,
               'nd_filter': nd_filter,
               'freq_center': freq_center,
               'freq_range': freq_range,
               'num_steps': num_steps,
               'num_runs': num_runs,
               'uwave_power': uwave_power,
               'readout': readout,
               'uwave_switch_delay': uwave_switch_delay,
               'sig_counts': sig_counts.astype(int).tolist(),
               'ref_counts': ref_counts.astype(int).tolist(),
               'norm_avg_sig': norm_avg_sig.astype(float).tolist()}

    filePath = tool_belt.get_file_path(file_name_no_ext, timestampEnd, name)
    tool_belt.save_figure(fig, filePath)
    tool_belt.save_raw_data(rawData, filePath)
