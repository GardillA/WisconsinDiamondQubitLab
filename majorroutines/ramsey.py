# -*- coding: utf-8 -*-
"""
Ramsey measruement.

This routine puts polarizes the nv state into 0, then applies a pi/2 pulse to 
put the state into a superposition between the 0 and + or - 1 state. The state
then evolves for a time, tau, of free precesion, and then a second pi/s pulse
is applied. The amount of population in 0 is read out by collecting the 
fluorescence during a readout.

Created on Wed Apr 24 15:01:04 2019

@author: agardill
"""

# %% Imports


import utils.tool_belt as tool_belt
import majorroutines.optimize as optimize
import numpy
import os
import time
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
from random import shuffle

#import json
#from scipy import asarray as ar,exp

# %% Main

def main(cxn, nv_sig, nd_filter, apd_indices,
         uwave_freq, uwave_power, uwave_pi_half_pulse, precession_time_range,
         num_steps, num_reps, num_runs, 
         name='untitled'):
    
    tool_belt.reset_cfm(cxn)
    
#    print(coords)
    
    # %% Defiene the times to be used in the sequence

    # Define some times (in ns)
    # time to intially polarize the nv
    polarization_time = 3 * 10**3
    # time of illumination during which signal readout occurs
    signal_time = 3 * 10**3
    # time of illumination during which reference readout occurs
    reference_time = 3 * 10**3
    # time between polarization and experiment without illumination
    pre_uwave_exp_wait_time = 1 * 10**3
    # time between the end of the experiment and signal without illumination
    post_uwave_exp_wait_time = 1 * 10**3
    # time between signal and reference without illumination
    sig_to_ref_wait_time = pre_uwave_exp_wait_time + post_uwave_exp_wait_time
    # the amount of time the AOM delays behind the gate and rf
    aom_delay_time = 750
    # the amount of time the rf delays behind the AOM and rf
    rf_delay_time = 40
    # the length of time the gate will be open to count photons
    gate_time = 320
    
    # Convert pi_pulse to integer
    uwave_pi_half_pulse = round(uwave_pi_half_pulse)
    
    seq_file_name = 't1_double_quantum.py'
    

    # %% Create the array of relaxation times
    
    # Array of times to sweep through
    # Must be ints since the pulse streamer only works with int64s
    
    min_precession_time = int(precession_time_range[0])
    max_precession_time = int(precession_time_range[1])
    
    taus = numpy.linspace(min_precession_time, max_precession_time,
                          num=num_steps, dtype=numpy.int32)
     
    # %% Fix the length of the sequence to account for odd amount of elements
     
    # Our sequence pairs the longest time with the shortest time, and steps 
    # toward the middle. This means we only step through half of the length
    # of the time array. 
    
    # That is a problem if the number of elements is odd. To fix this, we add 
    # one to the length of the array. When this number is halfed and turned 
    # into an integer, it will step through the middle element.
    
    if len(taus) % 2 == 0:
        half_length_taus = int( len(taus) / 2 )
    elif len(taus) % 2 == 1:
        half_length_taus = int( (len(taus) + 1) / 2 )
        
    # Then we must use this half length to calculate the list of integers to be
    # shuffled for each run
    
    tau_ind_list = list(range(0, half_length_taus))
        
    # %% Create data structure to save the counts
    
    # We create an array of NaNs that we'll fill
    # incrementally for the signal and reference counts. 
    # NaNs are ignored by matplotlib, which is why they're useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.
    
    sig_counts = numpy.empty([num_runs, num_steps], dtype=numpy.uint32)
    sig_counts[:] = numpy.nan
    ref_counts = numpy.copy(sig_counts)
    
    # %% Make some lists and variables to save at the end
    
    opti_coords_list = []
    
    
    # %% Analyze the sequence
    
    # We can simply reuse t1_double_quantum for this if we pass a pi/2 pulse
    # instead of a pi pulse and use the same states for init/readout
    seq_args = [min_precession_time, polarization_time, signal_time, reference_time, 
                sig_to_ref_wait_time, pre_uwave_exp_wait_time, 
                post_uwave_exp_wait_time, aom_delay_time, rf_delay_time, 
                gate_time, uwave_pi_half_pulse, 0,
                max_precession_time, apd_indices[0], 1, 1]
    ret_vals = cxn.pulse_streamer.stream_load(seq_file_name, seq_args, 1)
    seq_time = ret_vals[0]
#    print(sequence_args)
#    print(seq_time)
    
    # %% Let the user know how long this will take
    
    seq_time_s = seq_time / (10**9)  # s
    expected_run_time = num_steps * num_reps * num_runs * seq_time_s / 2  # s
    expected_run_time_m = expected_run_time / 60 # s
    
    print(' \nExpected run time: {:.1f} minutes. '.format(expected_run_time_m))
    
    # %% Get the starting time of the function, to be used to calculate run time

    startFunctionTime = time.time()
    
    # %% Collect the data
    
    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()
    
    for run_ind in range(num_runs):

        print(' \nRun index: {}'.format(run_ind))
        
        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break
        
        # Optimize
        opti_coords = optimize.main(cxn, nv_sig, nd_filter, apd_indices)
        opti_coords_list.append(opti_coords)
        
        # Set up the microwaves - just use the Tektronix
        cxn.signal_generator_tsg4104a.set_freq(uwave_freq)
        cxn.signal_generator_tsg4104a.set_amp(uwave_power)
        cxn.signal_generator_tsg4104a.uwave_on()
            
        # Load the APD
        cxn.apd_tagger.start_tag_stream(apd_indices)  
        
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

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break
            
            print(' \nFirst relaxation time: {}'.format(taus[tau_ind_first]))
            print('Second relaxation time: {}'.format(taus[tau_ind_second])) 
            
            seq_args = [taus[tau_ind_first], polarization_time, signal_time, reference_time, 
                            sig_to_ref_wait_time, pre_uwave_exp_wait_time, 
                            post_uwave_exp_wait_time, aom_delay_time, rf_delay_time, 
                            gate_time, uwave_pi_half_pulse, 0,
                            taus[tau_ind_second], apd_indices[0], 1, 1]
            
            cxn.pulse_streamer.stream_immediate(seq_file_name, num_reps,
                                                seq_args, 1)   
            
            # Each sample is of the form [*(<sig_shrt>, <ref_shrt>, <sig_long>, <ref_long>)]
            # So we can sum on the values for similar index modulus 4 to
            # parse the returned list into what we want.
            new_counts = cxn.apd_tagger.read_counter_separate_gates(1)
            sample_counts = new_counts[0]
            
            count = sum(sample_counts[0::4])
            sig_counts[run_ind, tau_ind_first] = count
            print('First signal = ' + str(count))
            
            count = sum(sample_counts[1::4])
            ref_counts[run_ind, tau_ind_first] = count  
            print('First Reference = ' + str(count))
            
            count = sum(sample_counts[2::4])
            sig_counts[run_ind, tau_ind_second] = count
            print('Second Signal = ' + str(count))

            count = sum(sample_counts[3::4])
            ref_counts[run_ind, tau_ind_second] = count
            print('Second Reference = ' + str(count))
            
        cxn.apd_tagger.stop_tag_stream()

    # %% Hardware clean up
    
    tool_belt.reset_cfm(cxn)
    
    # %% Average the counts over the iterations

    avg_sig_counts = numpy.average(sig_counts, axis=0)
    avg_ref_counts = numpy.average(ref_counts, axis=0)
    
    # %% Calculate the ramsey data, signal / reference over different 
    # relaxation times

    # Replace x/0=inf with 0
    try:
        norm_avg_sig = avg_sig_counts / avg_ref_counts
    except RuntimeWarning as e:
        print(e)
        inf_mask = numpy.isinf(norm_avg_sig)
        # Assign to 0 based on the passed conditional array
        norm_avg_sig[inf_mask] = 0
    
    # %% Plot the signal

    raw_fig, axes_pack = plt.subplots(1, 2, figsize=(17, 8.5))

    ax = axes_pack[0]
    ax.plot(taus / 10**3, avg_sig_counts, 'r-', label = 'signal')
    ax.plot(taus / 10**3, avg_ref_counts, 'g-', label = 'reference')
    ax.set_xlabel('Precession time (us)')
    ax.set_ylabel('Counts')
    ax.legend()

    ax = axes_pack[1]
    ax.plot(taus / 10**3, norm_avg_sig, 'b-')
    ax.set_title('Ramsey Measurement')
    ax.set_xlabel('Precession time (us)')
    ax.set_ylabel('Contrast (arb. units)')

    raw_fig.canvas.draw()
    # fig.set_tight_layout(True)
    raw_fig.canvas.flush_events()
    
    # %% Save the data

    endFunctionTime = time.time()

    timeElapsed = endFunctionTime - startFunctionTime

    timestamp = tool_belt.get_time_stamp()
    
    raw_data = {'timestamp': timestamp,
            'timeElapsed': timeElapsed,
            'name': name,
            'nv_sig': nv_sig,
            'nv_sig-units': tool_belt.get_nv_sig_units(),
            'nv_sig-format': tool_belt.get_nv_sig_format(),
            'opti_coords_list': opti_coords_list,
            'opti_coords_list-units': 'V',
            'nd_filter': nd_filter,
            'gate_time': gate_time,
            'gate_time-units': 'ns',
            'uwave_freq': uwave_freq,
            'uwave_freq-units': 'GHz',
            'uwave_power': uwave_power,
            'uwave_power-units': 'dBm',
            'uwave_pi_half_pulse': uwave_pi_half_pulse,
            'uwave_pi_half_pulse-units': 'ns',
            'precession_time_range': precession_time_range,
            'precession_time_range-units': 'ns',
            'num_steps': num_steps,
            'num_reps': num_reps,
            'num_runs': num_runs,
            'sig_counts': sig_counts.astype(int).tolist(),
            'sig_counts-units': 'counts',
            'ref_counts': ref_counts.astype(int).tolist(),
            'ref_counts-units': 'counts',
            'norm_avg_sig': norm_avg_sig.astype(float).tolist(),
            'norm_avg_sig-units': 'arb'}
    
    file_path = tool_belt.get_file_path(__file__, timestamp, name)
    tool_belt.save_figure(raw_fig, file_path)
    tool_belt.save_raw_data(raw_data, file_path)
    