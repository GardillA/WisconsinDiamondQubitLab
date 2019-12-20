# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 10:45:09 2019

This file is used to determine the cutoff for photon count nuer for individual
measurmenets between the charge states of the NV.

Collect the photon counts under yellow illumination, after reionizing NV into 
NV- with green light. A second collection occurs after ionizing NV to NV0 with 
red light.

@author: yanfeili
"""
import utils.tool_belt as tool_belt
import majorroutines.optimize as optimize
import numpy
import matplotlib.pyplot as plt
import labrad
import photonstatistics as ps

def get_Probability_distribution(aList):

    def get_unique_value(aList):
        unique_value_list = []
        for i in range(0,len(aList)):
            if aList[i] not in unique_value_list:
                unique_value_list.append(aList[i])
        return unique_value_list
    unique_value = get_unique_value(aList)
    relative_frequency = []
    for i in range(0,len(unique_value)):
        relative_frequency.append(aList.count(unique_value[i])/ (len(aList)))

    return unique_value, relative_frequency

#def create_figure():
#    fig, ax = plt.subplots(1, 1, figsize=(10, 8.5))
#    ax.set_xlabel('number of photons (n)')
#    ax.set_ylabel('P(n)')
#
#    return fig

#%% Main
# Connect to labrad in this file, as opposed to control panel
def main(nv_sig, apd_indices, aom_ao_589_pwr, ao_638_pwr, readout_time, 
         ionization_time, num_runs, num_reps):

    with labrad.connect() as cxn:
        main_with_cxn(cxn, nv_sig, apd_indices, aom_ao_589_pwr, ao_638_pwr, 
                      readout_time, ionization_time, num_runs, num_reps)

def main_with_cxn(cxn, nv_sig, apd_indices, aom_ao_589_pwr, ao_638_pwr, 
                  readout_time, ionization_time,num_runs, num_reps):

    tool_belt.reset_cfm(cxn)

# %% Initial Calculation and setup
#    apd_indices = [0]

    shared_params = tool_belt.get_shared_parameters_dict(cxn)
    
    #delay of aoms and laser
    aom_delay = shared_params['532_aom_delay']
    
    illumination_time = readout_time + 10**3
    reionization_time = 10**6
    
    readout_power = aom_ao_589_pwr

    # Set up our data structure, list
    # we repeatively collect photons for tR

    counts = []
    ref_counts = []
    sig_counts=[]
    opti_coords_list = []

#%% Estimate the lenth of the sequance

    seq_args = [readout_time, reionization_time, illumination_time, 
                ionization_time, aom_delay, apd_indices[0], 
                aom_ao_589_pwr, ao_638_pwr]
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    ret_vals = cxn.pulse_streamer.stream_load('determine_n_thresh_with_638.py', seq_args_string)

    seq_time = ret_vals[0]

    seq_time_s = seq_time / (10**9)  # s
    expected_run_time = num_reps * num_runs * seq_time_s + (0.5 * num_runs)  #s
    expected_run_time_m = expected_run_time / 60 # m

    # Ask to continue and timeout if no response in 2 seconds?

    print(' \nExpected run time: {:.1f} minutes. '.format(expected_run_time_m))

#    return

#%% Collect data
    tool_belt.init_safe_stop()


    for run_ind in range(num_runs):

        # Optimize
        opti_coords = optimize.main_with_cxn(cxn, nv_sig, apd_indices, 532)
        opti_coords_list.append(opti_coords)

        print('Run index: {}'. format(run_ind))

        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break

        # Load the APD
        cxn.apd_tagger.start_tag_stream(apd_indices)

        seq_args = [readout_time, reionization_time, illumination_time, 
                    ionization_time, aom_delay, apd_indices[0], 
                    aom_ao_589_pwr, ao_638_pwr]

        seq_args_string = tool_belt.encode_seq_args(seq_args)
        cxn.pulse_streamer.stream_immediate('determine_n_thresh_with_638.py', num_reps, seq_args_string)

        # Get the counts
        new_counts = cxn.apd_tagger.read_counter_simple(num_reps)

        counts.extend(new_counts)

    cxn.apd_tagger.stop_tag_stream()

    sig_counts = counts[0:len(counts):2]
    ref_counts = counts[1:len(counts):2]

#%% plot the data and the fit

    unique_value1, relative_frequency1 = get_Probability_distribution(list(sig_counts))
    unique_value2, relative_frequency2 = get_Probability_distribution(list(ref_counts))

    fig, ax = plt.subplots(1, 1, figsize=(10, 8.5))
    
    g00,g01,y01,y00 = ps.get_curve_fit_NV0(readout_time,readout_power,unique_value1, relative_frequency1)
    gm0,gm1,ym1,ym0 = ps.get_curve_fit_NVm(readout_time,readout_power,unique_value2, relative_frequency2)
    
    ax.plot(unique_value1, relative_frequency1, 'ro')
    ax.plot(unique_value2, relative_frequency2, 'bo')
    ax.set_xlabel('number of photons (n)')
    ax.set_ylabel('P(n)')
    
    text = '\n'.join(('Reionization time (532 nm)' + '%.3f'%(reionization_time/10**3) + 'us',
                      'Illumination time (589 nm)' + '%.3f'%(illumination_time/10**3) + 'us',
                      'Ionization time (638 nm)' + '%.3f'%(ionization_time/10**3) + 'us',
                      'Readout time' + '%.3f'%(readout_time/10**3)+ 'us'))

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.55, 0.6, text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)

#%% Save data
    timestamp = tool_belt.get_time_stamp()

    # turn the list of unique_values into pure integers, for saving
    unique_value1 = [int(el) for el in unique_value1]
    unique_value2 = [int(el) for el in unique_value2]
    relative_frequency1 = [int(el) for el in relative_frequency1]
    relative_frequency2 = [int(el) for el in relative_frequency2]
    sig_counts = [int(el) for el in sig_counts]
    ref_counts = [int(el) for el in ref_counts]

    raw_data = {'timestamp': timestamp,
            'nv_sig': nv_sig,
            'aom_ao_589_pwr': aom_ao_589_pwr,
            'aom_ao_589_pwr-units':'V',
            'ao_638_pwr': ao_638_pwr,
            'ao_638_pwr-units': 'V',
            'readout_time':readout_time,
            'readout_time_unit':'ns',
            'reionization_time': reionization_time,
            'reionization_time-units': 'ns',
            'illumination_time': illumination_time,
            'illumination_time-units': 'ns',
            'ionization_time': ionization_time,
            'ionization_time-units': 'ns',
            'nv_sig-units': tool_belt.get_nv_sig_units(),
            'num_runs': num_runs,
            'num_reps':num_reps,
            'sig_counts': sig_counts,
            'sig_counts-units': 'counts',
            'unique_valuesNV-': unique_value1,
            'unique_values-units': 'num of photons',
            'relative_frequencyNV-': relative_frequency1,
            'relative_frequency-units': 'occurrences',
            'unique_valuesNV0': unique_value2,
            'unique_values-units': 'num of photons',
            'relative_frequencyNV0': relative_frequency2,
            'relative_frequency-units': 'occurrences'
            }

    file_path = tool_belt.get_file_path(__file__, timestamp, nv_sig['name'])
    tool_belt.save_raw_data(raw_data, file_path)

    tool_belt.save_figure(fig, file_path)
