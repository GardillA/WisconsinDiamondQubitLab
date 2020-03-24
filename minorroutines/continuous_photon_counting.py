# -*- coding: utf-8 -*-
"""
This program allows the initial laser pulse and illumination pulse colors to be
specified between 532, 589, and 638. 

The apd then collects the photons duringreadout and plots the binned counts 
over time. This program is a reworked lifetime_v2.

Useful combination of init and illum lasers:
    init green and illum yellow
    init red and illum yellow
    init red and illum green

Created on Tue Mar 24 12:49:55 2020

@author: agardill
"""

# %% Imports


import utils.tool_belt as tool_belt
import majorroutines.optimize as optimize
import numpy
import os
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
import labrad


# %% Functions


def process_raw_buffer(new_tags, new_channels,
                       current_tags, current_channels,
                       gate_open_channel, gate_close_channel):
    
    # The processing here will be bin_size agnostic
    
    # Tack the new data onto the leftover data (leftovers are necessary if
    # the last read contained a gate open without a matching gate close)
    current_tags.extend(new_tags)
    current_channels.extend(new_channels)
    current_channels_array = numpy.array(current_channels)
    
    # Find gate open clicks
    result = numpy.nonzero(current_channels_array == gate_open_channel)
    gate_open_click_inds = result[0].tolist()

    # Find gate close clicks
    result = numpy.nonzero(current_channels_array == gate_close_channel)
    gate_close_click_inds = result[0].tolist()
    
    new_processed_tags = []
    
    # Loop over the number of closes we have since there are guaranteed to
    # be opens
    num_closed_samples = len(gate_close_click_inds)
    for list_ind in range(num_closed_samples):
        
        gate_open_click_ind = gate_open_click_inds[list_ind]
        gate_close_click_ind = gate_close_click_inds[list_ind]
        
        # Extract all the counts between these two indices as a single sample
        rep = current_tags[gate_open_click_ind+1:
                                    gate_close_click_ind]
        rep = numpy.array(rep, dtype=numpy.int64)
        # Make relative to gate open
        rep -= current_tags[gate_open_click_ind]
        new_processed_tags.extend(rep.astype(int).tolist())
        
    # Clear processed tags
    if len(gate_close_click_inds) > 0:
        leftover_start = gate_close_click_inds[-1]
        del current_tags[0: leftover_start+1]
        del current_channels[0: leftover_start+1]
        
    return new_processed_tags, num_closed_samples


# %% Main


def main(nv_sig, apd_indices, readout_time, init_pulse_duration,
                  illum_pulse_duration, aom_ao_589_pwr, ao_638_pwr, 
                  init_color_ind, illum_color_ind,
                  num_reps, num_runs, num_bins):

    with labrad.connect() as cxn:
        main_with_cxn(cxn, nv_sig, apd_indices, readout_time, init_pulse_duration,
                  illum_pulse_duration, aom_ao_589_pwr, ao_638_pwr, 
                  init_color_ind, illum_color_ind,
                  num_reps, num_runs, num_bins)


def main_with_cxn(cxn, nv_sig, apd_indices, readout_time, init_pulse_duration,
                  illum_pulse_duration, aom_ao_589_pwr, ao_638_pwr, 
                  init_color_ind, illum_color_ind,
                  num_reps, num_runs, num_bins):
    
    if len(apd_indices) > 1:
        msg = 'Currently lifetime only supports single APDs!!'
        raise NotImplementedError(msg)
    
    tool_belt.reset_cfm(cxn)

    # %% Define the times to be used in the sequence
    
    shared_params = tool_belt.get_shared_parameters_dict(cxn)

    # In ns
    readout_time = int(readout_time)
    init_pulse_duration = int(init_pulse_duration)
    illum_pulse_duration = int(illum_pulse_duration)

    aom_delay_time = shared_params['532_aom_delay']
    wait_time = shared_params['post_polarization_wait_dur']
    
    # %% Read the optical power forinit and illum pulse light
    
    
    if init_color_ind == 532:
        init_optical_power_pd = tool_belt.opt_power_via_photodiode(532)
        init_optical_power_mW = \
            tool_belt.calc_optical_power_mW(532, init_optical_power_pd)
    elif init_color_ind == 589:
        init_optical_power_pd = tool_belt.opt_power_via_photodiode(589,
           AO_power_settings = aom_ao_589_pwr, nd_filter = nv_sig['nd_filter'])        
        init_optical_power_mW = \
            tool_belt.calc_optical_power_mW(589, init_optical_power_pd)
    elif init_color_ind == 638:
        init_optical_power_pd = tool_belt.opt_power_via_photodiode(638,
           AO_power_settings = ao_638_pwr)        
        init_optical_power_mW = \
            tool_belt.calc_optical_power_mW(638, init_optical_power_pd)
     
    if illum_color_ind == 532:
        illum_optical_power_pd = tool_belt.opt_power_via_photodiode(532)
        illum_optical_power_mW = \
            tool_belt.calc_optical_power_mW(532, illum_optical_power_pd)        
    if illum_color_ind == 589:       
        illum_optical_power_pd = tool_belt.opt_power_via_photodiode(589,
           AO_power_settings = aom_ao_589_pwr, nd_filter = nv_sig['nd_filter']) 
        illum_optical_power_mW = \
            tool_belt.calc_optical_power_mW(589, illum_optical_power_pd)        
    elif illum_color_ind == 638:
        illum_optical_power_pd = tool_belt.opt_power_via_photodiode(638,
           AO_power_settings = ao_638_pwr)        
        illum_optical_power_mW = \
            tool_belt.calc_optical_power_mW(638, illum_optical_power_pd)

    # %% Analyze the sequence

    # pulls the file of the sequence from serves/timing/sequencelibrary
    file_name = os.path.basename(__file__)
    seq_args = [readout_time, init_pulse_duration, illum_pulse_duration, wait_time, 
                aom_delay_time, aom_ao_589_pwr, ao_638_pwr, apd_indices[0],
                init_color_ind, illum_color_ind]
    seq_args = [int(el) for el in seq_args]
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    ret_vals = cxn.pulse_streamer.stream_load(file_name, seq_args_string)
    seq_time = ret_vals[0]

    # %% Report the expected run time

    seq_time_s = seq_time / (10**9)  # s
    expected_run_time = num_reps * num_runs * seq_time_s  # s
    expected_run_time_m = expected_run_time / 60 # m
    print(' \nExpected run time: {:.1f} minutes. '.format(expected_run_time_m))
#    return

    # %% Bit more setup

    # Record the start time
    startFunctionTime = time.time()
    start_timestamp = tool_belt.get_time_stamp()
    
    opti_coords_list = []

    # %% Collect the data
    
    processed_tags = []

    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()

    for run_ind in range(num_runs):
        
        print(' \nRun index: {}'.format(run_ind))

        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break

        # Optimize
        opti_coords = optimize.main_with_cxn(cxn, nv_sig, apd_indices, 532, disable=True)
        opti_coords_list.append(opti_coords)
        
        
        # Expose the stream
        cxn.apd_tagger.start_tag_stream(apd_indices, apd_indices, False)
    
        # Find the gate channel
        # The order of channel_mapping is APD, APD gate open, APD gate close
        channel_mapping = cxn.apd_tagger.get_channel_mapping()
        gate_open_channel = channel_mapping[1]
        gate_close_channel = channel_mapping[2]
            
        # Stream the sequence
        seq_args = [readout_time, init_pulse_duration, illum_pulse_duration, wait_time, 
                aom_delay_time, aom_ao_589_pwr, ao_638_pwr, apd_indices[0],
                init_color_ind, illum_color_ind]
        seq_args = [int(el) for el in seq_args]
        seq_args_string = tool_belt.encode_seq_args(seq_args)
        
        cxn.pulse_streamer.stream_immediate(file_name, int(num_reps),
                                            seq_args_string)
        
        # Initialize state
        current_tags = []
        current_channels = []
        num_processed_reps = 0

        while num_processed_reps < num_reps:

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break
    
            new_tags, new_channels = cxn.apd_tagger.read_tag_stream()
            new_tags = numpy.array(new_tags, dtype=numpy.int64)
            
            ret_vals = process_raw_buffer(new_tags, new_channels,
                                   current_tags, current_channels,
                                   gate_open_channel, gate_close_channel)
            new_processed_tags, num_new_processed_reps = ret_vals
            # MCC test
            if num_new_processed_reps > 750000:
                print('Processed {} reps out of 10^6 max'.format(num_new_processed_reps))
                print('Tell Matt that the time tagger is too slow!')
            
            num_processed_reps += num_new_processed_reps
            
            processed_tags.extend(new_processed_tags)
            

        cxn.apd_tagger.stop_tag_stream()
        
        # %% Save the data we have incrementally for long measurements

        raw_data = {'start_timestamp': start_timestamp,
                    'nv_sig': nv_sig,
                    'nv_sig-units': tool_belt.get_nv_sig_units(),
                    'init_color_ind': init_color_ind,
                    'init_optical_power_pd': init_optical_power_pd,
                    'init_optical_power_pd-units': 'V',
                    'init_optical_power_mW': init_optical_power_mW,
                    'init_optical_power_mW-units': 'mW',
                    'illum_color_ind': illum_color_ind,
                    'illum_optical_power_pd': illum_optical_power_pd,
                    'illum_optical_power_pd-units': 'V',
                    'illum_optical_power_mW': illum_optical_power_mW,
                    'illum_optical_power_mW-units': 'mW',
                    'readout_time': readout_time,
                    'readout_time-units': 'ns',
                    'init_pulse_duration': init_pulse_duration,
                    'init_pulse_duration-units': 'ns',
                    'illum_pulse_duration': illum_pulse_duration,
                    'illum_pulse_duration-units': 'ns',
                    'num_reps': num_reps,
                    'num_runs': num_runs,
                    'run_ind': run_ind,
                    'opti_coords_list': opti_coords_list,
                    'opti_coords_list-units': 'V',
                    'processed_tags': processed_tags,
                    'processed_tags-units': 'ps',
                    }

        # This will continuously be the same file path so we will overwrite
        # the existing file with the latest version
        file_path = tool_belt.get_file_path(__file__, start_timestamp,
                                            nv_sig['name'], 'incremental')
        tool_belt.save_raw_data(raw_data, file_path)

    # %% Hardware clean up

    tool_belt.reset_cfm(cxn)

    # %% Bin the data
    
    readout_time_ps = 1000*readout_time
    
#    start_readout_time_ps = 1000*start_readout_time
#    end_readout_time_ps = 1000*end_readout_time
    binned_samples, bin_edges = numpy.histogram(processed_tags, num_bins,
                                (0, readout_time_ps))
#    print(binned_samples)
    
    # Compute the centers of the bins
    bin_size = readout_time / num_bins
    bin_center_offset = bin_size / 2
    bin_centers = numpy.linspace(0, readout_time, num_bins) + bin_center_offset
#    print(bin_centers)

    # %% Plot

    fig, ax = plt.subplots(1, 1, figsize=(10, 8.5))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    
    ax.plot(bin_centers, binned_samples, 'r-')
    ax.set_title('Lifetime')
    ax.set_xlabel('Readout time (ns)')
    ax.set_ylabel('Counts')
    ax.set_title('{} initial pulse, {} readout'.format(init_color_ind, illum_color_ind))
       
    params_text = '\n'.join(('Init pulse time: {} us'.format(init_pulse_duration/10**3),
                      'Init power: ' + '%.3f'%(init_optical_power_mW)+ 'mW',
                      'Illum power: ' + '%.3f'%(illum_optical_power_mW)+ 'mW',
                      'bin size: ' + '%.1f'%(bin_size) + 'ns'))

    ax.text(0.55, 0.85, params_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    fig.canvas.draw()
    fig.set_tight_layout(True)
    fig.canvas.flush_events()

    # %% Save the data

    endFunctionTime = time.time()
    time_elapsed = endFunctionTime - startFunctionTime
    timestamp = tool_belt.get_time_stamp()
    
    raw_data = {'timestamp': timestamp,
                'time_elapsed': time_elapsed,
                'nv_sig': nv_sig,
                'nv_sig-units': tool_belt.get_nv_sig_units(),
                'init_color_ind': init_color_ind,
                'init_optical_power_pd': init_optical_power_pd,
                'init_optical_power_pd-units': 'V',
                'init_optical_power_mW': init_optical_power_mW,
                'init_optical_power_mW-units': 'mW',
                'illum_color_ind': illum_color_ind,
                'illum_optical_power_pd': illum_optical_power_pd,
                'illum_optical_power_pd-units': 'V',
                'illum_optical_power_mW': illum_optical_power_mW,
                'illum_optical_power_mW-units': 'mW',
                'readout_time': readout_time,
                'readout_time-units': 'ns',
                'init_pulse_duration': init_pulse_duration,
                'init_pulse_duration-units': 'ns',
                'illum_pulse_duration': illum_pulse_duration,
                'illum_pulse_duration-units': 'ns',
                'num_bins': num_bins,
                'num_reps': num_reps,
                'num_runs': num_runs,
                'opti_coords_list': opti_coords_list,
                'opti_coords_list-units': 'V',
                'binned_samples': binned_samples.tolist(),
                'bin_centers': bin_centers.tolist(),
                'processed_tags': processed_tags,
                'processed_tags-units': 'ps',
                }

    file_path = tool_belt.get_file_path(__file__, timestamp, nv_sig['name'])
    tool_belt.save_figure(fig, file_path)
    tool_belt.save_raw_data(raw_data, file_path)

# %%

#def decayExp(t, amplitude, decay):
#    return amplitude * numpy.exp(- t / decay)
#
#def triple_decay(t, a1, d1, a2, d2, a3, d3):
#    return decayExp(t, a1, d1) + decayExp(t, a2, d2) + decayExp(t, a3, d3)
#
## %% Fitting the data
#
#def t1_exponential_decay(open_file_name):
#
#    directory = 'E:/Shared Drives/Kolkowitz Lab Group/nvdata/lifetime_v2/'
#
#    # Open the specified file
#    with open(directory + open_file_name + '.txt') as json_file:
#
#        # Load the data from the file
#        data = json.load(json_file)
#        countsT1_array = numpy.array(data["binned_samples"])
#        readout_time = data["readout_time"]
#        num_bins = data["num_bins"]
#        
#    bin_size = readout_time / num_bins
#    bin_center_offset = bin_size / 2
#    bin_centers = numpy.linspace(0, readout_time, num_bins) + bin_center_offset
#
#    amplitude = 500
#    decay = 100 # us
#    init_params = [amplitude, decay]
#    
#    init_params = [500, 10, 500, 100, 500, 500]
#    
#
##    popt,pcov = curve_fit(decayExp, timeArray, countsT1,
##                              p0=init_params)
##    popt,pcov = curve_fit(triple_decay, bin_centers, countsT1_array,
##                              p0=init_params)
#
#    linspaceTime = numpy.linspace(0, readout_time, num=1000)
#
#
#    fig_fit, ax= plt.subplots(1, 1, figsize=(10, 8))
#    ax.plot(bin_centers, countsT1_array,'bo',label='data')
##    ax.plot(linspaceTime, triple_decay(linspaceTime,*popt),'r-',label='fit')
#    ax.set_xlabel('Wait Time (ns)')
#    ax.set_ylabel('Counts (arb.)')
#    ax.set_title('Lifetime')
#    ax.legend()
#
##    text = "\n".join((r'$A_0 e^{-t / d}$',
##                      r'$A_0 = $' + '%.1f'%(popt[0]),
##                      r'$d = $' + "%.1f"%(decay_time) + " us"))
##    text = "\n".join((r'$A_1 e^{-t / d_1} + A_2 e^{-t / d_2} + A_3 e^{-t / d_3}$',
##                      r'$A_1 = $' + '%.1f'%(popt[0]),
##                      r'$d_1 = $' + "%.1f"%(popt[1]) + " us",
##                      r'$A_2 = $' + '%.1f'%(popt[2]),
##                      r'$d_2 = $' + "%.1f"%(popt[3]) + " us",
##                      r'$A_3 = $' + '%.1f'%(popt[4]),
##                      r'$d_3 = $' + "%.1f"%(popt[5]) + " us"))
#
#
##    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
##    ax.text(0.65, 0.75, text, transform=ax.transAxes, fontsize=12,
##                            verticalalignment="top", bbox=props)
#    ax.set_yscale("log", nonposy='clip')
#    
#    fig_fit.canvas.draw()
#    fig_fit.canvas.flush_events()
#
#    file_path = directory + open_file_name
#    tool_belt.save_figure(fig_fit, file_path+'-triple_fit_semilog')
##    fig_fit.savefig(open_file_name + '-replot.svg')

# %%
    

if __name__ == '__main__':
    file_name = '2019_11/2019_11_26-17_22_11-undoped_Y2O3-633_bandpass'
    
#    t1_exponential_decay(file_name)