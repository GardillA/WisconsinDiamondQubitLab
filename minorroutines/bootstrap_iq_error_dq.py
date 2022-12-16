# -*- coding: utf-8 -*-
"""
Dynamical decoupling CPMG



Created on Fri Aug 5 2022

@author: agardill
"""

# %% Imports


import utils.tool_belt as tool_belt
import majorroutines.optimize as optimize
import copy
from numpy import pi
import numpy
import time
import matplotlib.pyplot as plt
from random import shuffle
import labrad
from utils.tool_belt import States
from utils.tool_belt import NormStyle
from scipy.optimize import curve_fit


# %% Constants


im = 0 + 1j
inv_sqrt_2 = 1 / numpy.sqrt(2)
gmuB = 2.8e-3  # gyromagnetic ratio in GHz / G



# %% Main


def main(
    nv_sig,
    phase_range,
    num_steps,
    num_reps,
    num_runs,
    iq_state=States.HIGH,
):

    with labrad.connect() as cxn:
        angle = main_with_cxn(
            cxn,
            nv_sig,
            phase_range,
            num_steps,
            num_reps,
            num_runs,
            iq_state,
        )
        return angle


def main_with_cxn(
    cxn,
    nv_sig,
    phase_range,
    num_steps,
    num_reps,
    num_runs,
    iq_state=States.HIGH,
):

    counter_server = tool_belt.get_server_counter(cxn)
    pulsegen_server = tool_belt.get_server_pulse_gen(cxn)
    arbwavegen_server = tool_belt.get_server_arb_wave_gen(cxn)
    
    tool_belt.reset_cfm(cxn)

    # %% Sequence setup

    laser_key = "spin_laser"
    laser_name = nv_sig[laser_key]
    tool_belt.set_filter(cxn, nv_sig, laser_key)
    laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)
    polarization_time = nv_sig["spin_pol_dur"]
    gate_time = nv_sig["spin_readout_dur"]
    norm_style = nv_sig['norm_style']
    pi_pulse_reps = 1

    seq_file_name = "dynamical_decoupling_dq.py"
    
    rabi_period_low = nv_sig["rabi_{}".format(States.LOW.name)]
    uwave_freq_low = nv_sig["resonance_{}".format(States.LOW.name)]
    uwave_power_low = nv_sig["uwave_power_{}".format(States.LOW.name)]
    uwave_pi_pulse_low = tool_belt.get_pi_pulse_dur(rabi_period_low)
    uwave_pi_on_2_pulse_low = tool_belt.get_pi_on_2_pulse_dur(rabi_period_low)
    rabi_period_high = nv_sig["rabi_{}".format(States.HIGH.name)]
    uwave_freq_high = nv_sig["resonance_{}".format(States.HIGH.name)]
    uwave_power_high = nv_sig["uwave_power_{}".format(States.HIGH.name)]
    uwave_pi_pulse_high = tool_belt.get_pi_pulse_dur(rabi_period_high)
    uwave_pi_on_2_pulse_high = tool_belt.get_pi_on_2_pulse_dur(rabi_period_high)
    
    if iq_state.value == States.LOW.value:
        state_activ = States.LOW
        state_proxy = States.HIGH
        
    elif iq_state.value == States.HIGH.value:
        state_activ = States.HIGH
        state_proxy = States.LOW
        
    
    ### Create the array of the phases to test

    phis = numpy.linspace(
        phase_range[0],
        phase_range[1],
        num=num_steps,
        )
    
    phis_deg = phis*180/pi
    # print(phis_deg)
    # return

    if len(phis) % 2 == 0:
        half_length_phis = int(len(phis) / 2)
    elif len(phis) % 2 == 1:
        half_length_phis = int((len(phis) + 1) / 2)

    # Then we must use this half length to calculate the list of integers to be
    # shuffled for each run

    phi_ind_list = list(range(0, half_length_phis))

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
    phi_index_master_list = [[] for i in range(num_runs)]

    # %% Analyze the sequence
    
    num_reps = int(num_reps)
    
    seq_args = [
        100,
        polarization_time,
        gate_time,
        uwave_pi_pulse_low,
        uwave_pi_on_2_pulse_low,
        uwave_pi_pulse_high,
        uwave_pi_on_2_pulse_high,
        100,
        pi_pulse_reps,
        state_activ.value,
        state_proxy.value,
        laser_name, 
        laser_power
    ]
    seq_args_string = tool_belt.encode_seq_args(seq_args)
    ret_vals = pulsegen_server.stream_load(seq_file_name, seq_args_string)
    seq_time = ret_vals[0]
    # print(seq_args)
    # return
        # print(seq_time)

    # %% Let the user know how long this will take

    seq_time_s = seq_time / (10 ** 9)  # to seconds
    expected_run_time_s = (
        (num_steps / 2) * num_reps * num_runs * seq_time_s
    )  # s
    expected_run_time_m = expected_run_time_s / 60  # to minutes

    print(" \nExpected run time: {:.1f} minutes. ".format(expected_run_time_m))
    # return
    
    # create figure
    raw_fig, axes_pack = plt.subplots(1, 2, figsize=(17, 8.5))
    
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
        sig_gen_low_cxn = tool_belt.get_server_sig_gen(cxn, States.LOW)
        sig_gen_low_cxn.set_freq(uwave_freq_low)
        sig_gen_low_cxn.set_amp(uwave_power_low)
        try:
            sig_gen_low_cxn.load_iq()
        except Exception:
            pass
        sig_gen_low_cxn.uwave_on()
        sig_gen_high_cxn = tool_belt.get_server_sig_gen(cxn, States.HIGH)
        sig_gen_high_cxn.set_freq(uwave_freq_high)
        sig_gen_high_cxn.set_amp(uwave_power_high)
        try:
            sig_gen_high_cxn.load_iq()
        except Exception:
            pass
        sig_gen_high_cxn.uwave_on()

        # Set up the laser
        tool_belt.set_filter(cxn, nv_sig, laser_key)
        laser_power = tool_belt.set_laser_power(cxn, nv_sig, laser_key)

        # Load the APD
        counter_server.start_tag_stream()

        # Shuffle the list of tau indices so that it steps thru them randomly
        shuffle(phi_ind_list)

        for phi_ind in phi_ind_list:

            # 'Flip a coin' to determine which tau (long/shrt) is used first
            rand_boolean = numpy.random.randint(0, high=2)

            if rand_boolean == 1:
                phi_ind_first = phi_ind
                phi_ind_second = -phi_ind - 1
            elif rand_boolean == 0:
                phi_ind_first = -phi_ind - 1
                phi_ind_second = phi_ind

            # add the tau indexxes used to a list to save at the end
            phi_index_master_list[run_ind].append(phi_ind_first)
            phi_index_master_list[run_ind].append(phi_ind_second)

            # Break out of the while if the user says stop
            if tool_belt.safe_stop():
                break

            print(" \nFirst phase: {}".format(phis[phi_ind_first]*180/pi))
            print("Second relaxation time: {}".format(phis[phi_ind_second]*180/pi))

            arbwavegen_server.load_arb_phases([0, phis[phi_ind_first], 0, 0, phis[phi_ind_second], 0])
            # arbwavegen_server.load_cpmg(1)
            
            # Clear the tagger buffer of any excess counts
            # counter_server.clear_buffer()
            pulsegen_server.stream_immediate(
                seq_file_name, num_reps, seq_args_string
            )

            # Each sample is of the form [*(<sig_shrt>, <ref_shrt>, <sig_long>, <ref_long>)]
            # So we can sum on the values for similar index modulus 4 to
            # parse the returned list into what we want.
            new_counts = counter_server.read_counter_separate_gates(1)
            sample_counts = new_counts[0]
            # print(new_counts)

            count = sum(sample_counts[0::4])
            sig_counts[run_ind, phi_ind_first] = count

            count = sum(sample_counts[1::4])
            ref_counts[run_ind, phi_ind_first] = count

            count = sum(sample_counts[2::4])
            sig_counts[run_ind, phi_ind_second] = count

            count = sum(sample_counts[3::4])
            ref_counts[run_ind, phi_ind_second] = count
            
            arbwavegen_server.reset()

        counter_server.stop_tag_stream()

        # %% incremental plotting
      
        # Average the counts over the iterations
        inc_sig_counts = sig_counts[: run_ind + 1]
        inc_ref_counts = ref_counts[: run_ind + 1]
        ret_vals = tool_belt.process_counts(
            inc_sig_counts, inc_ref_counts, num_reps, gate_time, norm_style
        )
        (
            sig_counts_avg_kcps,
            ref_counts_avg_kcps,
            norm_avg_sig,
            norm_avg_sig_ste,
        ) = ret_vals
        
        
        
        ax = axes_pack[0]
        ax.cla()
        ax.plot(phis_deg, sig_counts_avg_kcps, "r-", label="signal")
        ax.plot(phis_deg, ref_counts_avg_kcps, "g-", label="reference")
        ax.set_xlabel(r"Relative phase, $\phi$ (degrees)")
        ax.set_ylabel("kcps")
        ax.legend()
        
        ax = axes_pack[1]
        ax.cla()
        ax.plot(phis_deg, norm_avg_sig, "b-")
        # ax.set_title("CPMG-{} Measurement".format(pi_pulse_reps))
        ax.set_xlabel(r"Relative phase, $\phi$ (degrees)")
        ax.set_ylabel("Contrast (arb. units)")
        
        text_popt = 'Run # {}/{}'.format(run_ind+1,num_runs)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.8, 0.9, text_popt,transform=ax.transAxes,
                verticalalignment='top', bbox=props)
        
        raw_fig.canvas.draw()
        raw_fig.set_tight_layout(True)
        raw_fig.canvas.flush_events()
        
        # %% Save the data we have incrementally for long T1s

        raw_data = {
            "start_timestamp": start_timestamp,
            "nv_sig": nv_sig,
            "nv_sig-units": tool_belt.get_nv_sig_units(cxn),
            "phase_range": phase_range,
            "phase_range-units": "radians",
            "iq_state": iq_state.name,
            "num_steps": num_steps,
            "num_reps": num_reps,
            "run_ind": run_ind,	
            "phis": phis.tolist(),	
            "phis_deg":phis_deg.tolist(),	
            "phi_index_master_list": phi_index_master_list,
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

    # %% Hardware clean up

    tool_belt.reset_cfm(cxn)

    # %% Plot the data

    ret_vals = tool_belt.process_counts(sig_counts, ref_counts, num_reps, gate_time, norm_style)
    (
        sig_counts_avg_kcps,
        ref_counts_avg_kcps,
        norm_avg_sig,
        norm_avg_sig_ste,
    ) = ret_vals
    
    ax = axes_pack[0]
    ax.cla()
    ax.plot(phis_deg, sig_counts_avg_kcps, "r-", label="signal")
    ax.plot(phis_deg, ref_counts_avg_kcps, "g-", label="reference")
    ax.set_xlabel(r"Relative phase, $\phi$ (degrees)")
    ax.set_ylabel("kcps")
    ax.legend()

    ax = axes_pack[1]
    ax.cla()
    ax.plot(phis_deg, norm_avg_sig, "b-")
    ax.set_title(r"$\pi / 2_x$ - $\pi_{\phi}$ - $\pi / 2_x$ (DQ basis)")
    ax.set_xlabel(r"Relative phase, $\phi$ (degrees)")
    ax.set_ylabel("Contrast (arb. units)")

    raw_fig.canvas.draw()
    raw_fig.set_tight_layout(True)
    raw_fig.canvas.flush_events()

    # %% Save the data

    endFunctionTime = time.time()

    timeElapsed = endFunctionTime - startFunctionTime

    timestamp = tool_belt.get_time_stamp()

    raw_data = {
        "timestamp": timestamp,
        "timeElapsed": timeElapsed,
        "nv_sig": nv_sig,
        "nv_sig-units": tool_belt.get_nv_sig_units(cxn),
        "phase_range": phase_range,
        "phase_range-units": "radians",
        "iq_state": iq_state.name,
        "num_steps": num_steps,
        "num_reps": num_reps,
        "run_ind": run_ind,	
        "phis": phis.tolist(),	
        "phis_deg":phis_deg.tolist(),	
        "phi_index_master_list": phi_index_master_list,
        "opti_coords_list": opti_coords_list,
        "opti_coords_list-units": "V",
        "sig_counts": sig_counts.astype(int).tolist(),
        "sig_counts-units": "counts",
        "ref_counts": ref_counts.astype(int).tolist(),
        "ref_counts-units": "counts",
        "norm_avg_sig": norm_avg_sig.astype(float).tolist(),
        "norm_avg_sig-units": "arb",
    }

    nv_name = nv_sig["name"]
    file_path = tool_belt.get_file_path(__file__, timestamp, nv_name)
    tool_belt.save_figure(raw_fig, file_path)
    tool_belt.save_raw_data(raw_data, file_path)

    # Fit and save figs


    return 


    
# %% Run the file


if __name__ == "__main__":

    nd_yellow = "nd_0"
    green_power =8000
    nd_green = 'nd_1.1'
    red_power = 120
    sample_name = "siena"
    # sample_name = "hopper"
    green_laser = "integrated_520"
    yellow_laser = "laserglow_589"
    red_laser = "cobolt_638"
    
    
    sig_base = {
        "disable_opt":False,
        "ramp_voltages": False,
        # "correction_collar": 0.12,
        'expected_count_rate': None,
    
        "spin_laser":green_laser,
        "spin_laser_power": green_power,
        "spin_laser_filter": nd_green,
        "spin_readout_dur": 300,
        "spin_pol_dur": 10000.0,
    
        "imaging_laser":green_laser,
        "imaging_laser_power": green_power,
        "imaging_laser_filter": nd_green,
        "imaging_readout_dur": 1e7,
        
        # "imaging_laser":yellow_laser,
        # "imaging_laser_power": 0.2,
        # "imaging_laser_filter": "nd_1.0",
        # "imaging_readout_dur": 5e7,
    
        "initialize_laser": green_laser,
        "initialize_laser_power": green_power,
        "initialize_laser_dur":  1e3,
        "CPG_laser": green_laser,
        "CPG_laser_power":red_power,
        "CPG_laser_dur": int(1e6),
    
        "nv-_prep_laser": green_laser,
        "nv-_prep_laser-power": None,
        "nv-_prep_laser_dur": 1e3,
        "nv0_prep_laser": red_laser,
        "nv0_prep_laser-power": None,
        "nv0_prep_laser_dur": 1e3,
        
        "nv-_reionization_laser": green_laser,
        "nv-_reionization_laser_power": green_power,
        "nv-_reionization_dur": 1e3,
        
        "nv0_ionization_laser": red_laser,
        "nv0_ionization_laser_power": None,
        "nv0_ionization_dur": 300,
        
        "spin_shelf_laser": yellow_laser,
        "spin_shelf_laser_power": None,
        "spin_shelf_dur": 0,
        
        "charge_readout_laser": yellow_laser,
        "charge_readout_laser_power": 0.15, 
        "charge_readout_laser_filter": nd_yellow,
        "charge_readout_dur": 200e6, 
    
        "collection_filter": "715_sp+630_lp", # NV band only
        "uwave_power_LOW": 12.12,  
        "uwave_power_HIGH": 10,
        
    } 
    
    
    
    nv_sig_1 = copy.deepcopy(sig_base) # 
    nv_sig_1["coords"] = [-0.184, 0.092, 4.05]
    nv_sig_1["name"] = "{}-nv1_2022_10_27".format(sample_name,)
    # nv_sig_1["norm_style"]= NormStyle.POINT_TO_POINT
    nv_sig_1["norm_style"]= NormStyle.SINGLE_VALUED
    nv_sig_1[ "green_power_mW"] = 1.0
    nv_sig_1["expected_count_rate"] = 19
    nv_sig_1[ "spin_readout_dur"] = 300
    nv_sig_1['magnet_angle'] = 151.7
    
    nv_sig_1["resonance_LOW"]= 2.7805
    nv_sig_1["rabi_LOW"]=136.09 # +/- 0.54
    nv_sig_1["resonance_HIGH"]=2.9600
    nv_sig_1["rabi_HIGH"]= 176.3

        

    phase_range = [0, 2*pi]
    num_steps = 25
    num_reps = 2e4
    num_runs = 10
    main(
          nv_sig_1,
          phase_range,
          num_steps,
          num_reps,
          num_runs,
          iq_state=States.HIGH,
      )
    
    
    # file_name = "2022_12_16-14_32_57-siena-nv1_2022_10_27"
    # data = tool_belt.get_raw_data(file_name, 'pc_rabi/branch_master/bootstrap_iq_error_dq/2022_12')
    # norm_avg_sig = data['norm_avg_sig']
    # # norm_avg_sig_ste = data['norm_avg_sig_ste']
    # phis_deg = data['phis_deg']
    
    # tau_lin = numpy.linspace(phis_deg[0], phis_deg[-1], 1000)
    
    # fig, ax = plt.subplots()
    # ax.plot(phis_deg, norm_avg_sig,  "bo")
    
    # fit_func = lambda t, amp, offset, freq:tool_belt.sin_1_at_0_phase(t, amp, offset, freq, pi)
    # init_params = [0.08, 0.9, 5/360]
    # popt, pcov = curve_fit(
    #     fit_func,
    #     phis_deg,
    #     norm_avg_sig,
    #     p0=init_params,
    # )
    # print('Offset = {}'.format(popt[1]))
    # print('Amplitude = {}'.format(popt[0]))
    # ax.plot(
    #         tau_lin,
    #         fit_func(tau_lin, *popt),
    #         "r-",
    #         label="fit",
    #     ) 
    
    # import utils.kplotlib as kpl
    # kpl.init_kplotlib()
    
    # base_text = "Offset = {:.3f} \nAmp = {:.3f}"
    # size = kpl.Size.SMALL
    # text = base_text.format(popt[1], popt[0])
    # kpl.anchored_text(ax, text, kpl.Loc.UPPER_RIGHT, size=size)
    # ax.set_xlabel(r"Relative phase, $\phi$ (degrees)")
    # ax.set_ylabel("Contrast (arb. units)")
    
    
