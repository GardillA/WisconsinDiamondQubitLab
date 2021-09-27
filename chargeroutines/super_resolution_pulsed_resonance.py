# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 10:52:28 2021

@author: agardill
"""

# %% Imports


import utils.tool_belt as tool_belt
import majorroutines.optimize as optimize
import numpy
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import labrad
from utils.tool_belt import States
import majorroutines.pulsed_resonance as pulsed_resonance
from random import shuffle

# %%

def build_voltages(adjusted_nv_coords, adjusted_depletion_coords, num_reps):
    start_x_value = adjusted_nv_coords[0]
    start_y_value = adjusted_nv_coords[1]
    
    dep_x_value = adjusted_depletion_coords[0]
    dep_y_value = adjusted_depletion_coords[1]
    
    # we want the sequence to look like the following: 
        # [[nv_coords], [dep_coords], [nv_coords], [nv_coords], 
        #               [dep_coords], [nv_coords], [nv_coords], 
        #                       .... ]
   

    seq_x = [dep_x_value, start_x_value, start_x_value]
    seq_y = [dep_y_value, start_y_value, start_y_value]
    
    x_voltages = seq_x*num_reps*2
    y_voltages = seq_y*num_reps*2
    
    # and then add on the initial coordinate
    x_voltages = [start_x_value] + x_voltages
    y_voltages = [start_y_value] + y_voltages
    
    return x_voltages, y_voltages

def plot_esr(ref_counts, sig_counts, num_runs, freqs = None, freq_center = None, freq_range = None, num_steps = None):
    
    # if all.freqs() == None:
    #     half_freq_range = freq_range / 2
    #     freq_low = freq_center - half_freq_range
    #     freq_high = freq_center + half_freq_range
    #     freqs = numpy.linspace(freq_low, freq_high, num_steps)
    

    ret_vals = pulsed_resonance.process_counts(ref_counts, sig_counts, num_runs)
    avg_ref_counts, avg_sig_counts, norm_avg_sig, ste_ref_counts, ste_sig_counts, norm_avg_sig_ste = ret_vals

    # Convert to kilocounts per second
    # readout_sec = depletion_time / 1e9
    cts_uwave_off_avg = (avg_ref_counts / (num_reps))# * 1000)) / readout_sec
    cts_uwave_on_avg = (avg_sig_counts / (num_reps))# * 1000)) / readout_sec

    # Create an image with 2 plots on one row, with a specified size
    # Then draw the canvas and flush all the previous plots from the canvas
    fig, axes_pack = plt.subplots(1, 2, figsize=(17, 8.5))

    if len(freqs) == 1:
        marker = 'o'
    else:
        marker = '-'
    # The first plot will display both the uwave_off and uwave_off counts
    ax = axes_pack[0]
    ax.plot(freqs, cts_uwave_off_avg, 'r{}'.format(marker), label = 'Reference')
    ax.plot(freqs, cts_uwave_on_avg, 'g{}'.format(marker), label = 'Signal')
    ax.set_title('Non-normalized Count Rate Versus Frequency')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('NV fluorescence (counts)')
    ax.legend()
    # The second plot will show their subtracted values
    ax = axes_pack[1]
    ax.plot(freqs, norm_avg_sig, 'b{}'.format(marker))
    ax.set_title('Normalized Count Rate vs Frequency')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Contrast (arb. units)')

    fig.canvas.draw()
    fig.tight_layout()
    fig.canvas.flush_events()
    
    return fig, norm_avg_sig, norm_avg_sig_ste 
# %% Main


def main(nv_sig, opti_nv_sig, apd_indices, freq_center, freq_range,
         num_steps, num_reps, num_runs, uwave_power, uwave_pulse_dur,
         state=States.LOW):

    with labrad.connect() as cxn:
        resonance_list = main_with_cxn(cxn, nv_sig,opti_nv_sig, apd_indices, freq_center, freq_range,
                  num_steps, num_reps, num_runs, uwave_power, uwave_pulse_dur,
                  state)
    return resonance_list


def main_with_cxn(cxn, nv_sig, opti_nv_sig,apd_indices, freq_center, freq_range,
              num_steps, num_reps, num_runs, uwave_power, uwave_pulse_dur,
              state=States.LOW):

    # %% Initial calculations and setup

    tool_belt.reset_cfm(cxn)

    # Calculate the frequencies we need to set
    half_freq_range = freq_range / 2
    freq_low = freq_center - half_freq_range
    freq_high = freq_center + half_freq_range
    freqs = numpy.linspace(freq_low, freq_high, num_steps)
    # freqs = numpy.array([freq_center])
    freq_ind_list = list(range(num_steps))
    
    opti_interval = 4 # min
    
    nv_coords = nv_sig['coords']
    depletion_coords = nv_sig['depletion_coords']

    # Set up our data structure, an array of NaNs that we'll fill
    # incrementally. NaNs are ignored by matplotlib, which is why they're
    # useful for us here.
    # We define 2D arrays, with the horizontal dimension for the frequency and
    # the veritical dimension for the index of the run.
    ref_counts = numpy.empty([num_runs, num_steps])
    ref_counts[:] = numpy.nan
    sig_counts = numpy.copy(ref_counts)
    
    # imaging_laser_key = 'imaging_laser'
    # imaging_laser_name = nv_sig[imaging_laser_key]
    # imaging_laser_power = tool_belt.set_laser_power(cxn, nv_sig, imaging_laser_key)
    
    init_color = tool_belt.get_registry_entry_no_cxn('wavelength',
                      ['Config', 'Optics', nv_sig['initialize_laser']])
    depletion_color = tool_belt.get_registry_entry_no_cxn('wavelength',
                      ['Config', 'Optics', nv_sig['CPG_laser']])
    
    
    # Set the charge readout (assumed to be yellow here) to the correct filter
    if 'charge_readout_laser_filter' in nv_sig:
        tool_belt.set_filter(cxn, nv_sig, 'charge_readout_laser')
    
    readout_time = nv_sig['charge_readout_dur']
    init_time = nv_sig['initialize_dur']
    depletion_time = nv_sig['CPG_laser_dur']
    readout_power = nv_sig['charge_readout_laser_power']
    ionization_time = nv_sig['nv0_ionization_dur']
    shelf_time = nv_sig['spin_shelf_dur']
    shelf_power = nv_sig['spin_shelf_laser_power']
    
    
    green_laser_name = nv_sig['imaging_laser']
    red_laser_name = nv_sig['nv0_ionization_laser']
    yellow_laser_name = nv_sig['charge_readout_laser']
    sig_gen_name = tool_belt.get_signal_generator_name_no_cxn(state)    
            
    seq_args = [readout_time, init_time, depletion_time, ionization_time, uwave_pulse_dur, shelf_time,
            uwave_pulse_dur,init_color, depletion_color, 
            green_laser_name, yellow_laser_name, red_laser_name, sig_gen_name, 
             apd_indices[0], readout_power, shelf_power ]
    # print(seq_args)
    # return
    seq_args_string = tool_belt.encode_seq_args(seq_args)

    drift_list = []

    # %% Get the starting time of the function

    start_timestamp = tool_belt.get_time_stamp()

    # %% Collect the data

    # Start 'Press enter to stop...'
    tool_belt.init_safe_stop()
    
    start_time = time.time()
    start_function_time = start_time
    
    print(depletion_coords)
    
    for run_ind in range(num_runs):
        print('Run index: {}'. format(run_ind))

        # Break out of the while if the user says stop
        if tool_belt.safe_stop():
            break

        # Optimize and save the coords we found
        optimize.main_with_cxn(cxn, opti_nv_sig, apd_indices)
        drift = tool_belt.get_drift()
        drift_list.append(drift)
        
        adjusted_nv_coords = numpy.array(nv_coords) + drift
        adjusted_depletion_coords = numpy.array(depletion_coords) + drift
        
        
        # Set up the microwaves and laser. Then load the pulse streamer 
        # (must happen after optimize and iq_switch since run their
        # own sequences)
        sig_gen_cxn = tool_belt.get_signal_generator_cxn(cxn, state)
        sig_gen_cxn.set_amp(uwave_power)
        ret_vals = cxn.pulse_streamer.stream_load('super_resolution_rabi.py', seq_args_string)
        
        period = ret_vals[0]
        
        period_s = period/10**9
        period_s_total = (period_s*num_reps*num_steps + 1)
        period_m_total = period_s_total/60
        print('Expected time for this run: {:.1f} m'.format(period_m_total))
        # return

        # Shuffle the freqs we step thru
        shuffle(freq_ind_list)
        
        # Take a sample and increment the frequency
        for step_ind in range(num_steps):
            
            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break
            
            freq_ind = freq_ind_list[step_ind]
            print(freqs[freq_ind])
            
            time_current = time.time()
            if time_current - start_time > opti_interval * 60:
                optimize.main_with_cxn(cxn, opti_nv_sig, apd_indices)
                drift = tool_belt.get_drift()
                drift_list.append(drift)
                
                adjusted_nv_coords = numpy.array(nv_coords) + drift
                adjusted_depletion_coords = numpy.array(depletion_coords) + drift
                
                start_time = time_current
                
            tool_belt.set_xyz(cxn, adjusted_nv_coords)
            
            # Build the list to step through the coords on readout NV and targets
            x_voltages, y_voltages = build_voltages(adjusted_nv_coords, 
                                                  adjusted_depletion_coords, num_reps)
            
            # print(freqs[freq_ind])
            sig_gen_cxn.set_freq(freqs[freq_ind])
            sig_gen_cxn.uwave_on()

            # Start the tagger stream
            cxn.apd_tagger.start_tag_stream(apd_indices)
            cxn.pulse_streamer.stream_load('super_resolution_rabi.py', seq_args_string)
            
            # Load the galvo
            xy_server = tool_belt.get_xy_server(cxn) 
            xy_server.load_arb_scan_xy(x_voltages, y_voltages, int(period))
            
            # Clear the tagger buffer of any excess counts
            # cxn.apd_tagger.clear_buffer()
            
            # Start the timing stream
            cxn.pulse_streamer.stream_start(int(num_reps))

            num_samples = num_reps * 6
            num_read_so_far = 0
            total_samples_list = []
            while num_read_so_far < num_samples:
        
                if tool_belt.safe_stop():
                    break
        
                # Read the samples and update the image
                new_samples = cxn.apd_tagger.read_counter_simple()
                num_new_samples = len(new_samples)
        
                if num_new_samples > 0:
                    for el in new_samples:
                        total_samples_list.append(int(el))
                    num_read_so_far += num_new_samples
            # print(total_samples_list)
            sig_gate_counts = total_samples_list[2::6]
            # sig_gate_counts = total_samples_list[1::6]
            sig_counts[run_ind, freq_ind] = sum(sig_gate_counts)

            ref_gate_counts = total_samples_list[5::6]
            # ref_gate_counts = total_samples_list[4::6]
            ref_counts[run_ind, freq_ind] = sum(ref_gate_counts)

        cxn.apd_tagger.stop_tag_stream()

        # %% Save the data we have incrementally for long measurements

        rawData = {'start_timestamp': start_timestamp,
                   'nv_sig': nv_sig,
                   'nv_sig-units': tool_belt.get_nv_sig_units(),
                   'freq_center': freq_center,
                   'freq_center-units': 'GHz',
                   'freq_range': freq_range,
                   'freq_range-units': 'GHz',
                   'uwave_pulse_dur': uwave_pulse_dur,
                   'uwave_pulse_dur-units': 'ns',
                   'state': state.name,
                   'num_steps': num_steps,
                   'run_ind': run_ind,
                   'uwave_power': uwave_power,
                   'uwave_power-units': 'dBm',
                   'freqs': freqs.tolist(),
                   'drift_list': drift_list,
                   'opti_interval': opti_interval,
                   'sig_counts': sig_counts.astype(int).tolist(),
                   'sig_counts-units': 'counts',
                   'ref_counts': ref_counts.astype(int).tolist(),
                   'ref_counts-units': 'counts'}

        # This will continuously be the same file path so we will overwrite
        # the existing file with the latest version
        file_path = tool_belt.get_file_path(__file__, start_timestamp,
                                            nv_sig['name'], 'incremental')
        tool_belt.save_raw_data(rawData, file_path)


    # %% Process and plot the data
    end_function_time = time.time()
    time_elapsed = end_function_time - start_function_time
    # print(time_elapsed)
    
    fig, norm_avg_sig,norm_avg_sig_ste  = plot_esr(ref_counts, sig_counts, num_runs, freqs)


    # %% Fit the data

    # fit_func, popt, pcov = fit_resonance(freq_range, freq_center, num_steps,
    #                                       norm_avg_sig, norm_avg_sig_ste)
    # if (fit_func is not None) and (popt is not None):
    #     fit_fig = create_fit_figure(freq_range, freq_center, num_steps,
    #                                 norm_avg_sig, fit_func, popt)
    # else:
    #     fit_fig = None

    # %% Clean up and save the data

    tool_belt.reset_cfm(cxn)

    timestamp = tool_belt.get_time_stamp()

    rawData = {'timestamp': timestamp,
               'time_elapsed': time_elapsed,
                'nv_sig': nv_sig,
                'nv_sig-units': tool_belt.get_nv_sig_units(),
                'freq_center': freq_center,
                'freq_center-units': 'GHz',
                'freq_range': freq_range,
                'freq_range-units': 'GHz',
                'uwave_pulse_dur': uwave_pulse_dur,
                'uwave_pulse_dur-units': 'ns',
                'state': state.name,
                'num_steps': num_steps,
                'num_reps': num_reps,
                'num_runs': num_runs,
                'uwave_power': uwave_power,
                'uwave_power-units': 'dBm',
                'freqs': freqs.tolist(),
                'drift_list': drift_list,
                'opti_interval': opti_interval,
                'sig_counts': sig_counts.astype(int).tolist(),
                'sig_counts-units': 'counts',
                'ref_counts': ref_counts.astype(int).tolist(),
                'ref_counts-units': 'counts',
                'norm_avg_sig': norm_avg_sig.astype(float).tolist(),
                'norm_avg_sig-units': 'arb',
                'norm_avg_sig_ste': norm_avg_sig_ste.astype(float).tolist(),
                'norm_avg_sig_ste-units': 'arb'}

    name = nv_sig['name']
    filePath = tool_belt.get_file_path(__file__, timestamp, name)
    tool_belt.save_figure(fig, filePath)
    tool_belt.save_raw_data(rawData, filePath)
    # filePath = tool_belt.get_file_path(__file__, timestamp, name + '-fit')
    # if fit_fig is not None:
    #     tool_belt.save_figure(fit_fig, filePath)

    # # %% Return

    # if fit_func == single_gaussian_dip:
    #     print('Single resonance at {:.4f} GHz'.format(popt[2]))
    #     print('\n')
    #     return popt[2], None
    # elif fit_func == double_gaussian_dip:
    #     print('Resonances at {:.4f} GHz and {:.4f} GHz'.format(popt[2], popt[5]))
    #     print('Splitting of {:d} MHz'.format(int((popt[5] - popt[2]) * 1000)))
    #     print('\n')
    #     return popt[2], popt[5]
    # else:
    #     print('No resonances found')
    #     print('\n')
    #     return None, None
    

# %%

if __name__ == '__main__':

    apd_indices = [0]
    sample_name = 'johnson'    
    
    green_laser = "cobolt_515"
    yellow_laser = 'laserglow_589'
    red_laser = 'cobolt_638'
    
    green_power = 7
    red_power = 120
    nd_yellow = "nd_0.5"
    
    opti_nv_sig = {
        "coords": [-0.033, -0.021, 4.85],
        "name": "{}-nv1_2021_09_07".format(sample_name,),
        "disable_opt": False,
        "expected_count_rate": 15,
        "imaging_laser":green_laser,
        "imaging_laser_power": green_power,
        "imaging_readout_dur": 1e7,
        "collection_filter": "630_lp",
        "magnet_angle": None,
    }  # 14.5 max
    
    
    nv_sig = {
        "coords": [-0.020, 0.282, 4.85],
        "depletion_coords": [-0.032, 0.275, 4.85], # A
        #"depletion_coords": [-0.035, 0.276, 4.85], # B
        
        "name": "{}-dnv7_2021_09_23".format(sample_name,),
        "disable_opt": False,
        "expected_count_rate": 65,
            'imaging_laser': green_laser, 'imaging_laser_power': green_power,
            'imaging_readout_dur': 1E7,
            
            "initialize_laser": green_laser,
            "initialize_laser_power": green_power,
            "initialize_dur": 1e4,
            
            "CPG_laser": red_laser,
            'CPG_laser_power': red_power,
            "CPG_laser_dur": 1.5e4,
        
            'nv0_ionization_laser': red_laser, 'nv0_ionization_laser_power': red_power,
            'nv0_ionization_dur':300,
            
            'spin_shelf_laser': yellow_laser, 'spin_shelf_laser_filter': nd_yellow, 
            'spin_shelf_laser_power': 0.4, 'spin_shelf_dur':0,
            
            'charge_readout_laser': yellow_laser, 'charge_readout_laser_filter': nd_yellow, 
            'charge_readout_laser_power': 0.15, 'charge_readout_dur':50e6,
            
            'collection_filter': '630_lp', 'magnet_angle': 60,
            
            # "resonance_LOW": 2.8351, "rabi_LOW": 169.0, 'uwave_power_LOW': 15.5,  # 15.5 max
            "resonance_LOW": 2.806, "rabi_LOW":169.0, 'uwave_power_LOW': 15.5,  # 15.5 max
            'resonance_HIGH': 2.9445, 'rabi_HIGH': 191.9, 'uwave_power_HIGH': 14.5}   # 14.5 max  
    
    
    freq_range = 0.075
    
    uwave_power = nv_sig['uwave_power_LOW']
    # uwave_pulse_dur =  nv_sig['rabi_LOW'] / 2
    num_steps = 31
    num_reps = int(10**3)
    num_runs = 6
    
    A = [ -0.0318, 0.275, 4.85]
    B = [-0.035, 0.276, 4.85]
    do_plot = False
    
    try:
        # nv_sig['depletion_coords'] = A 
        # nv_sig['resonance_LOW'] = 2.8375
        # nv_sig['rabi_LOW']  = 127.4
        # main(nv_sig, nv_sig, apd_indices, nv_sig['resonance_LOW'], freq_range,
        #           num_steps, num_reps, num_runs, uwave_power, nv_sig['rabi_LOW']/2)
         
        # nv_sig['depletion_coords'] =  B
        # nv_sig['resonance_LOW'] = 2.8375
        # nv_sig['rabi_LOW']  = 89.2
        # main(nv_sig, nv_sig, apd_indices, nv_sig['resonance_LOW'], freq_range,
        #         num_steps, num_reps, num_runs, uwave_power, nv_sig['rabi_LOW']/2)
         
        
        
        if do_plot:
            folder = 'pc_rabi/branch_master/super_resolution_pulsed_resonance/2021_09'
            # ++++++ Add data sets
            file_list_a = [
                '2021_09_25-19_51_47-johnson-dnv7_2021_09_23',
                '2021_09_26-03_10_22-johnson-dnv7_2021_09_23',
                '2021_09_26-14_49_45-johnson-dnv7_2021_09_23',
                '2021_09_26-16_38_47-johnson-dnv7_2021_09_23',
                '2021_09_26-19_39_45-johnson-dnv7_2021_09_23',
                '2021_09_27-03_21_01-johnson-dnv7_2021_09_23',
                ]
            
            file_list_b = [
                '2021_09_25-21_39_13-johnson-dnv7_2021_09_23',
                '2021_09_26-07_39_00-johnson-dnv7_2021_09_23',
                '2021_09_26-15_43_30-johnson-dnv7_2021_09_23',
                '2021_09_26-17_32_32-johnson-dnv7_2021_09_23',
                '2021_09_26-21_27_12-johnson-dnv7_2021_09_23',
                '2021_09_27-08_43_31-johnson-dnv7_2021_09_23',
                ]
            
            total_sig_counts_A = []
            total_ref_counts_A = []
            num_runs_a = 0
            for file in file_list_a:
                data = tool_belt.get_raw_data(file, folder)
        
                freqs = data['freqs']
                num_steps = data['num_steps']
                num_runs = data['num_runs']
                sig_counts = data['sig_counts']
                ref_counts = data['ref_counts']
                total_sig_counts_A = total_sig_counts_A + sig_counts
                total_ref_counts_A = total_ref_counts_A + ref_counts
                num_runs_a = num_runs_a + num_runs
                
            
            fig, norm_avg_sig_a, norm_avg_sig_ste  = plot_esr(total_ref_counts_A, total_sig_counts_A, num_runs_a, freqs, None, None, num_steps)
             
            timestamp = tool_belt.get_time_stamp()
            freq_center = data['freq_center']
            freq_range = data['freq_range']
            uwave_pulse_dur = data['uwave_pulse_dur']
            state = data['state']
            opti_interval = data['opti_interval']
            rawData = {'timestamp': timestamp,
                    'nv_sig': nv_sig,
                    'nv_sig-units': tool_belt.get_nv_sig_units(),
                    'freq_center': freq_center,
                    'freq_center-units': 'GHz',
                    'freq_range': freq_range,
                    'freq_range-units': 'GHz',
                    'uwave_pulse_dur': uwave_pulse_dur,
                    'uwave_pulse_dur-units': 'ns',
                    'state': state,
                    'num_steps': num_steps,
                    'num_reps': num_reps,
                    'num_runs': num_runs_a,
                    'uwave_power': uwave_power,
                    'uwave_power-units': 'dBm',
                    'freqs': freqs,
                    'opti_interval': opti_interval,
                    'sig_counts': total_sig_counts_A,
                    'sig_counts-units': 'counts',
                    'ref_counts': total_ref_counts_A,
                    'ref_counts-units': 'counts',
                    'norm_avg_sig': norm_avg_sig_a.astype(float).tolist(),
                    'norm_avg_sig-units': 'arb',
                    'norm_avg_sig_ste': norm_avg_sig_ste.astype(float).tolist(),
                    'norm_avg_sig_ste-units': 'arb'}
            
            name = nv_sig['name']
            filePath = tool_belt.get_file_path(__file__, timestamp, name)
            tool_belt.save_figure(fig, filePath + '-compilation')
            tool_belt.save_raw_data(rawData, filePath + '-compilation')
             
            total_sig_counts_B = []
            total_ref_counts_B = []
            num_runs_b = 0
            for file in file_list_b:
                data = tool_belt.get_raw_data(file, folder)
        
                freqs = data['freqs']
                num_steps = data['num_steps']
                num_runs = data['num_runs']
                sig_counts = data['sig_counts']
                ref_counts = data['ref_counts']
                total_sig_counts_B = total_sig_counts_B + sig_counts
                total_ref_counts_B = total_ref_counts_B + ref_counts
                num_runs_b = num_runs_b + num_runs
                
            fig, norm_avg_sig_b, norm_avg_sig_ste  =  plot_esr(total_ref_counts_B, total_sig_counts_B, num_runs_b, freqs, None, None, num_steps)
            
            time.sleep(1)
            timestamp = tool_belt.get_time_stamp()
            freq_center = data['freq_center']
            freq_range = data['freq_range']
            uwave_pulse_dur = data['uwave_pulse_dur']
            state = data['state']
            opti_interval = data['opti_interval']
            rawData = {'timestamp': timestamp,
                    'nv_sig': nv_sig,
                    'nv_sig-units': tool_belt.get_nv_sig_units(),
                    'freq_center': freq_center,
                    'freq_center-units': 'GHz',
                    'freq_range': freq_range,
                    'freq_range-units': 'GHz',
                    'uwave_pulse_dur': uwave_pulse_dur,
                    'uwave_pulse_dur-units': 'ns',
                    'state': state,
                    'num_steps': num_steps,
                    'num_reps': num_reps,
                    'num_runs': num_runs_b,
                    'uwave_power': uwave_power,
                    'uwave_power-units': 'dBm',
                    'freqs': freqs,
                    'opti_interval': opti_interval,
                    'sig_counts': total_sig_counts_B,
                    'sig_counts-units': 'counts',
                    'ref_counts': total_ref_counts_B,
                    'ref_counts-units': 'counts',
                    'norm_avg_sig': norm_avg_sig_b.astype(float).tolist(),
                    'norm_avg_sig-units': 'arb',
                    'norm_avg_sig_ste': norm_avg_sig_ste.astype(float).tolist(),
                    'norm_avg_sig_ste-units': 'arb'}
            
            name = nv_sig['name']
            filePath = tool_belt.get_file_path(__file__, timestamp, name)
            tool_belt.save_figure(fig, filePath + '-compilation')
            tool_belt.save_raw_data(rawData, filePath + '-compilation')
            
            fig, ax = plt.subplots()
            ax.plot(freqs, norm_avg_sig_a, 'b-', label = "Point A")
            ax.plot(freqs, norm_avg_sig_b, 'r-', label = "Point B")
            ax.set_title('Normalized Count Rate vs Frequency')
            ax.set_xlabel('Frequency (GHz)')
            ax.set_ylabel('Contrast (arb. units)')
            ax.legend()
            
            
            # ++++ COMPARE +++++
            file_list = ['2021_09_23-03_06_43-johnson-dnv0_2021_09_09',
                          '2021_09_23-08_13_46-johnson-dnv0_2021_09_09']
            label_list = ['Point A', 'Point B']
                
            # fig, ax = plt.subplots(figsize=(8.5, 8.5))
            # for f in [0]:#range(len(file_list)):
            #     file = file_list[f]
            #     data = tool_belt.get_raw_data(file, folder)
        
            #     freq_center = data['freq_center']
            #     freq_range = data['freq_range']
            #     num_steps = data['num_steps']
            #     num_runs = data['num_runs']
            #     norm_avg_sig = data['norm_avg_sig']
                
                
                
            #     freqs = pulsed_resonance.calculate_freqs(freq_range, freq_center, num_steps)
            
            #     smooth_freqs, rel_counts = pulsed_resonance.simulate(2.8351, 
            #                                             freq_range, 0.04,
            #                                             169.0, 169.0/2)
            
            #     ax.plot(freqs, norm_avg_sig, 'b-', label='data')#label_list[f])
            #     ax.plot(smooth_freqs, rel_counts,'r--',  label='simulation')#label_list[f])
            #     ax.set_xlabel('Frequency (GHz)')
            #     ax.set_ylabel('Contrast (arb. units)')
            #     ax.legend(loc='lower right')
            
            # # +++++++ REPLOT ++++++++
            # file = '2021_09_23-03_06_43-johnson-dnv0_2021_09_09'
            # # file = '2021_09_23-08_13_46-johnson-dnv0_2021_09_09'
            # data = tool_belt.get_raw_data(file, folder)
            # ref_counts = data['ref_counts']
            # sig_counts = data['sig_counts']
            # num_runs = data['num_runs']
            # freq_center = data['freq_center']
            # freq_range = data['freq_range']
            # num_steps = data['num_steps']
            
            # plot_esr(ref_counts, sig_counts, num_runs, None, freq_center, freq_range, num_steps)
        
    
    finally:
        # Reset our hardware - this should be done in each routine, but
        # let's double check here
        # tool_belt.reset_cfm()
        # Kill safe stop
        if tool_belt.check_safe_stop_alive():
            print('\n\nRoutine complete. Press enter to exit.')
            tool_belt.poll_safe_stop()
