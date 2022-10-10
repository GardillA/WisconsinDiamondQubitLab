# -*- coding: utf-8 -*-
"""
This file contains functions to control the CFM. Just change the function call
in the main section at the bottom of this file and run the file. Shared or
frequently changed parameters are in the __main__ body and relatively static
parameters are in the function definitions.

Created on Sun Nov 25 14:00:28 2018

@author: mccambria
"""


# %% Imports


import labrad
import numpy
import numpy as np
import time
import copy
import utils.tool_belt as tool_belt
import majorroutines.image_sample as image_sample
import majorroutines.image_sample_temperature as image_sample_temperature
import majorroutines.map_rabi_contrast_NIR as map_rabi_contrast_NIR
import majorroutines.ensemble_image_sample_NIR_differential as ensemble_image_sample_NIR_differential
import majorroutines.ensemble_image_sample_NIR_differential_faster as ensemble_image_sample_NIR_differential_faster
import majorroutines.optimize as optimize
import majorroutines.stationary_count as stationary_count
import majorroutines.resonance as resonance
import majorroutines.pulsed_resonance as pulsed_resonance
import majorroutines.four_point_esr as four_point_esr
import majorroutines.optimize_magnet_angle as optimize_magnet_angle
import majorroutines.rabi as rabi
import majorroutines.discrete_rabi as discrete_rabi
import majorroutines.g2_measurement as g2_measurement
import majorroutines.t1_dq_main as t1_dq_main
import majorroutines.ramsey as ramsey
import majorroutines.spin_echo as spin_echo
import majorroutines.lifetime as lifetime
import majorroutines.lifetime_v2 as lifetime_v2
import chargeroutines.determine_charge_readout_params as determine_charge_readout_params
import chargeroutines.determine_charge_readout_params_moving_target as determine_charge_readout_params_moving_target
import chargeroutines.determine_charge_readout_params_1Dscan_target as determine_charge_readout_params_1Dscan_target
import minorroutines.determine_standard_readout_params as determine_standard_readout_params
import chargeroutines.scc_pulsed_resonance as scc_pulsed_resonance
import debug.test_major_routines as test_major_routines
from utils.tool_belt import States
import time


# %% Major Routines


def do_image_sample(nv_sig, apd_indices, 
                    nv_minus_initialization=False,cbarmin=None,cbarmax=None):

    scan_range = 0.2
    num_steps = 60
    
    # scan_range = 0.5
    # num_steps = 90

    # scan_range = 0.3
    # num_steps = 25

    # scan_range = 1.0
    # num_steps = 120

    # scan_range = 5.0
    # scan_range = 3.0
    # scan_range = 1.5
    # scan_range = 1.2
    # scan_range = 0.75
    # scan_range = 0.3
    # scan_range = 0.6
    # scan_range = 0.15
    # scan_range = 0.1
    # scan_range = 0.075
    #    scan_range = 0.025

    # num_steps = 300
    # num_steps = 200
    # num_steps = 150
    #    num_steps = 135
    # num_steps = 120
    # num_steps = 90
    # num_steps = 30
    # num_steps = 50
    # num_steps = 20

    # For now we only support square scans so pass scan_range twice
    image_sample.main(
        nv_sig,
        scan_range,
        scan_range,
        num_steps,
        apd_indices,
        nv_minus_initialization=nv_minus_initialization,
        cmin=cbarmin,
        cmax=cbarmax,
    )


def do_image_sample_zoom(nv_sig, apd_indices):

    scan_range = 0.2
    num_steps = 35

    image_sample.main(
        nv_sig,
        scan_range,
        scan_range,
        num_steps,
        apd_indices,
    )


def do_image_sample_temperature(nv_sig, apd_indices):
    
    # Sample every 0.02 V for diffraction limited spacing
    image_range = 0.2
    # num_steps = 5
    num_steps = 10
    
    nir_laser_voltage = 1.3
    
    # esr_num_reps = 3e4
    # esr_num_runs = 800
    esr_num_reps = 1e5
    esr_num_runs = 4
    
    image_sample_temperature.main(
        nv_sig,
        image_range,
        num_steps,
        apd_indices,
        nir_laser_voltage,
        esr_num_reps,
        esr_num_runs,
    )
    

def do_ensemble_image_sample(nv_sig, apd_indices, num_steps,  technique, scan_range=0.2, nv_minus_initialization=False,
                             cbarmin_percdiff=None,cbarmax_percdiff=None,cbarmin_counts=None,cbarmax_counts=None,cbarmin_diffcounts=None,cbarmax_diffcounts=None):
    
    list_of_techniquies = ['NIR_temperature_diff','NIR_counts_diff','NIR_counts_diff_faster'] 
    if technique not in list_of_techniquies:
        print('Technique must be in:', list_of_techniquies)
    
    if technique=='NIR_temperature_diff':
        # image_range = scan_range
        # num_steps = 10
        nir_laser_voltage = 1.3
        esr_num_reps = 1e5
        esr_num_runs = 4
        image_sample_temperature.main(
            nv_sig,
            scan_range,
            num_steps,
            apd_indices,
            nir_laser_voltage,
            esr_num_reps,
            esr_num_runs,
        )
    
    if technique=='NIR_counts_diff':
        # scan_range = 0.5
        # num_steps = 90
        # scan_range = .6
        # num_steps = 3
        nir_laser_voltage = 1.3
        ensemble_image_sample_NIR_differential.main(
            nv_sig,
            scan_range,
            num_steps,
            apd_indices,
            nir_laser_voltage,
            nv_minus_initialization
            )
        
    if technique=='NIR_counts_diff_faster':
        # scan_range = 0.5
        # num_steps = 90
        # scan_range = .6
        # num_steps = 3
        nir_laser_voltage = 1.3
        sleeptime=1
        ensemble_image_sample_NIR_differential_faster.main(
            nv_sig,
            scan_range,
            num_steps,
            apd_indices, 
            nir_laser_voltage,
            sleeptime,
            nv_minus_initialization,
            cbarmin_percdiff,cbarmax_percdiff,cbarmin_counts,cbarmax_counts,cbarmin_diffcounts,cbarmax_diffcounts)
        


def do_optimize(nv_sig, apd_indices):

    optimize.main(
        nv_sig,
        apd_indices,
        set_to_opti_coords=False,
        save_data=True,
        plot_data=True,
    )


def do_stationary_count(
    nv_sig,
    apd_indices,
    disable_opt=None,
    nv_minus_initialization=False,
    nv_zero_initialization=False,
):

    run_time = 3 * 60 * 10 ** 9  # ns

    stationary_count.main(
        nv_sig,
        run_time,
        apd_indices,
        disable_opt=disable_opt,
        nv_minus_initialization=nv_minus_initialization,
        nv_zero_initialization=nv_zero_initialization,
    )


def do_g2_measurement(nv_sig, apd_a_index, apd_b_index):

    run_time = 60 * 10  # s
    diff_window = 200  # ns

    g2_measurement.main(
        nv_sig, run_time, diff_window, apd_a_index, apd_b_index
    )


def do_resonance(nv_sig, apd_indices, freq_center=2.87, freq_range=0.2):

    num_steps = 51
    num_runs = 20
    uwave_power = -5.0

    resonance.main(
        nv_sig,
        apd_indices,
        freq_center,
        freq_range,
        num_steps,
        num_runs,
        uwave_power,
        state=States.HIGH,
    )


def do_resonance_state(nv_sig, apd_indices, state):

    freq_center = nv_sig["resonance_{}".format(state.name)]
    uwave_power = -5.0

    #    freq_range = 0.200
    #    num_steps = 51
    #    num_runs = 2

    # Zoom
    freq_range = 0.05
    num_steps = 51
    num_runs = 10

    resonance.main(
        nv_sig,
        apd_indices,
        freq_center,
        freq_range,
        num_steps,
        num_runs,
        uwave_power,
    )


def do_four_point_esr(nv_sig, apd_indices, state):

    detuning=0.004
    d_omega=0.002
    num_reps = 1e5
    num_runs = 4

    ret_vals = four_point_esr.main(
        nv_sig,
        apd_indices,
        num_reps,
        num_runs,
        state,
        detuning,
        d_omega,
        ret_file_name=True,
    )
    
    # print(resonance, res_err)
    return ret_vals


def do_determine_standard_readout_params(nv_sig, apd_indices):
    
    num_reps = 1e5
    max_readouts = [50e3]
    state = States.LOW
    
    determine_standard_readout_params.main(nv_sig, apd_indices, num_reps, 
                                           max_readouts, state=state)


def do_pulsed_resonance(nv_sig, apd_indices, freq_center=2.87, freq_range=0.2):

    num_steps = 51
    # num_reps = 1e5
    num_reps = 4e3
    num_runs = 4
    uwave_power = 16.5
    uwave_pulse_dur = 120

    pulsed_resonance.main(
        nv_sig,
        apd_indices,
        freq_center,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
        uwave_power,
        uwave_pulse_dur,
    )


def do_pulsed_resonance_state(nv_sig, apd_indices, state):

    freq_range = 0.040
    num_steps = 51
    # num_reps = 1e5
    # num_runs = 10
    # num_reps = 5e4
    # num_runs = 20
    num_reps = 4e3
    # num_runs = 16
    num_runs = 4

    # Zoom
    # freq_range = 0.035
    # # freq_range = 0.120
    # num_steps = 51
    # num_reps = 8000
    # num_runs = 3

    composite = False

    res, _ = pulsed_resonance.state(
        nv_sig,
        apd_indices,
        state,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
        composite,
    )
    nv_sig["resonance_{}".format(state.name)] = res
    return res


def do_scc_pulsed_resonance(nv_sig, apd_indices, state=States.LOW):

    opti_nv_sig = nv_sig
    state = "LOW"
    freq_center = nv_sig["resonance_{}".format(state)]
    uwave_power = nv_sig["uwave_power_{}".format(state)]
    uwave_pulse_dur = tool_belt.get_pi_pulse_dur(
        nv_sig["rabi_{}".format(state)]
    )
    freq_range = 0.040
    # num_steps = 21
    num_steps = 1
    num_reps = int(1e3)
    # num_runs = 80
    num_runs = 5

    # for red_dur in numpy.linspace(75, 300, 10):
    #     nv_sig['nv0_ionization_dur'] = red_dur
    #     scc_pulsed_resonance.main(nv_sig, opti_nv_sig, apd_indices,
    #                               freq_center, freq_range,
    #                               num_steps, num_reps, num_runs,
    #                               uwave_power, uwave_pulse_dur)

    scc_pulsed_resonance.main(
        nv_sig,
        opti_nv_sig,
        apd_indices,
        freq_center,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
        uwave_power,
        uwave_pulse_dur,
    )


def do_determine_charge_readout_params(nv_sig, apd_indices,nbins=None,nreps=500):

    # readout_durs = [10*10**3, 50*10**3, 100*10**3, 500*10**3,
    #                 1*10**6, 2*10**6, 3*10**6, 4*10**6, 5*10**6,
    #                 6*10**6, 7*10**6, 8*10**6, 9*10**6, 1*10**7,
    #                 2*10**7, 3*10**7, 4*10**7, 5*10**7]
    # readout_durs = numpy.linspace(10e6, 50e6, 5)
    # readout_durs = [10e6, 25e6, 50e6, 100e6, 200e6, 400e6, 700e6, 1e9, 2e9]
    # readout_durs = [10e6, 25e6, 50e6, 100e6, 200e6, 400e6, 1e9]
    readout_durs = [100e6]
    # readout_durs = [5e6, 10e6, 20e6, 40e6, 100e6, 500e6, 1e9]
    # readout_durs = numpy.linspace(700e6, 1e9, 7)
    # readout_durs = [50e6, 100e6, 200e6, 400e6, 1e9]
    # readout_durs = [2e9]
    readout_durs = [int(el) for el in readout_durs]
    max_readout_dur = max(readout_durs)

    # readout_powers = np.linspace(0.6, 0.8, 5)
    # readout_powers = np.arange(0.6, 1.05, 0.05)
    # readout_powers = np.arange(0.68, 1.04, 0.04)
    # readout_powers = np.linspace(0.9, 1.0, 3)
    # readout_powers = [0.25,0.5,0.75,1.0]
    readout_powers = [1.0]
    # readout_powers = [0.75]
    readout_powers = [round(val, 3) for val in readout_powers]

    # num_reps = 2000
    # num_reps = 1000
    num_reps = nreps

    determine_charge_readout_params.determine_readout_dur_power(
        nv_sig,
        nv_sig,
        apd_indices,
        num_reps,
        max_readout_dur=max_readout_dur,
        bins=nbins,
        readout_powers=readout_powers,
        plot_readout_durs=readout_durs,
    )

 
def do_determine_charge_readout_params_1Dscan_target(nv_sig, apd_indices, x_steps, y_steps, nbins=None,nreps=500,NIRtest=True):
    #this only works for a 1d scan. one of x_steps or y_steps should be a list or array with several points. 
    #The other one should be a list or array with one point, the one to stay at.
    
    readout_dur = [int(100e6)]
    readout_power = [round(1.0,3)]
    
    determine_charge_readout_params_1Dscan_target.main(
        nv_sig,
        nv_sig,
        np.asarray(x_steps), 
        np.asarray(y_steps), 
        apd_indices,
        nreps,
        max_readout_dur=max(readout_dur),
        bins=nbins,
        readout_powers=readout_power,
        plot_readout_durs=readout_dur,
        fit_threshold_full_model= False,
        NIR_test = NIRtest,
    )
    
            

def do_optimize_magnet_angle(nv_sig, apd_indices):

    angle_range = [0, 150]
    num_angle_steps = 6
    freq_center = 2.87
    freq_range = 0.200
    num_freq_steps = 51
    # num_freq_runs = 30
    num_freq_runs = 15

    # Pulsed
    uwave_power = 16.5
    uwave_pulse_dur = 85
    num_freq_reps = 5000

    # CW
    # uwave_power = -5.0
    # uwave_pulse_dur = None
    # num_freq_reps = None

    optimize_magnet_angle.main(
        nv_sig,
        apd_indices,
        angle_range,
        num_angle_steps,
        freq_center,
        freq_range,
        num_freq_steps,
        num_freq_reps,
        num_freq_runs,
        uwave_power,
        uwave_pulse_dur,
    )


def do_rabi(nv_sig, apd_indices, state, uwave_time_range=[0, 200]):

    num_steps = 51
    # num_reps = 2e4
    # num_runs = 20
    # num_reps = 1e5x
    # # num_runs = 5
    # num_runs = 10
    num_reps = 4e3
    # num_runs = 16
    num_runs = 4

    period = rabi.main(
        nv_sig,
        apd_indices,
        uwave_time_range,
        state,
        num_steps,
        num_reps,
        num_runs,
    )
    nv_sig["rabi_{}".format(state.name)] = period


def do_discrete_rabi(nv_sig, apd_indices, state, max_num_pi_pulses=4):

    # num_reps = 2e4
    num_reps = 1e3
    # num_reps = 100
    # num_runs = 2
    num_runs = 4

    discrete_rabi.main(
        nv_sig, apd_indices, state, max_num_pi_pulses, num_reps, num_runs
    )

    # for iq_delay in numpy.linspace(670, 680, 11):
    #     discrete_rabi.main(nv_sig, apd_indices,
    #                         state, max_num_pi_pulses, num_reps, num_runs, iq_delay)


def paper_figure1_data(nv_sig, apd_indices):
    # T1 experiment parameters, formatted:
    # [[init state, read state], relaxation_time_range, num_steps, num_reps]
    num_runs = 800
    num_reps = 2000
    num_steps = 12
    min_tau = int(500e3)
    max_tau = int(15e6)
    tau_linspace = numpy.linspace(min_tau, max_tau, num_steps)
    num_steps = num_steps - 2
    max_tau = tau_linspace[-3]
    t1_exp_array = numpy.array(
        [
            [
                [States.LOW, States.LOW],
                [min_tau, max_tau],
                num_steps,
                num_reps,
                num_runs,
            ],
            # [[States.LOW, States.HIGH], [min_tau, max_tau], num_steps, num_reps, num_runs],
            # [[States.ZERO, States.ZERO], [min_tau, max_tau], num_steps, num_reps, num_runs],
        ],
        dtype=object,
    )

    t1_dq_main.main(nv_sig, apd_indices, t1_exp_array, num_runs)


def do_t1_dq_scc(nv_sig, apd_indices):
    # T1 experiment parameters, formatted:
    # [[init state, read state], relaxation_time_range, num_steps, num_reps]
    num_runs = 400
    num_reps = 20
    num_steps = 6
    min_tau = 500e3
    max_tau_omega = int(12e9)
    max_tau_gamma = int(5e9)
    # max_tau_omega = int(5.3e9)
    # max_tau_gamma = int(3e9)
    t1_exp_array = numpy.array(
        [
            [
                [States.ZERO, States.HIGH],
                [min_tau, max_tau_omega],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.ZERO],
                [min_tau, max_tau_omega],
                num_steps,
                num_reps,
                num_runs,
            ],
            # [[States.ZERO, States.HIGH], [min_tau, max_tau_omega//3], num_steps, num_reps, num_runs],
            # [[States.ZERO, States.ZERO], [min_tau, max_tau_omega//3], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.LOW], [min_tau, max_tau_gamma], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.HIGH], [min_tau, max_tau_gamma], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.LOW], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.HIGH], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            [
                [States.LOW, States.HIGH],
                [min_tau, max_tau_gamma],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.LOW],
                [min_tau, max_tau_gamma],
                num_steps,
                num_reps,
                num_runs,
            ],
            # [[States.LOW, States.HIGH], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            # [[States.LOW, States.LOW], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
        ],
        dtype=object,
    )

    t1_dq_main.main(
        nv_sig, apd_indices, t1_exp_array, num_runs, scc_readout=True
    )


def do_t1_dq(nv_sig, apd_indices):
    # T1 experiment parameters, formatted:
    # [[init state, read state], relaxation_time_range, num_steps, num_reps]
    num_runs = 500
    num_reps = 1000
    num_steps = 12
    min_tau = 10e3
    max_tau_omega = int(18e6)
    max_tau_gamma = int(8.5e6)
    t1_exp_array = numpy.array(
        [
            [
                [States.ZERO, States.HIGH],
                [min_tau, max_tau_omega],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.ZERO],
                [min_tau, max_tau_omega],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.HIGH],
                [min_tau, max_tau_omega // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.ZERO],
                [min_tau, max_tau_omega // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
            # [[States.HIGH, States.LOW], [min_tau, max_tau_gamma], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.HIGH], [min_tau, max_tau_gamma], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.LOW], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.HIGH], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            [
                [States.LOW, States.HIGH],
                [min_tau, max_tau_gamma],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.LOW],
                [min_tau, max_tau_gamma],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.HIGH],
                [min_tau, max_tau_gamma // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.LOW],
                [min_tau, max_tau_gamma // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
        ],
        dtype=object,
    )

    t1_dq_main.main(nv_sig, apd_indices, t1_exp_array, num_runs)


def do_t1_dq_knill(nv_sig, apd_indices):
    # T1 experiment parameters, formatted:
    # [[init state, read state], relaxation_time_range, num_steps, num_reps]
    num_runs = 100
    num_reps = 1000
    num_steps = 12
    min_tau = 10e3
    max_tau_omega = int(18e6)
    max_tau_gamma = int(8.5e6)
    t1_exp_array = numpy.array(
        [
            [
                [States.ZERO, States.HIGH],
                [min_tau, max_tau_omega],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.ZERO],
                [min_tau, max_tau_omega],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.HIGH],
                [min_tau, max_tau_omega // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.ZERO, States.ZERO],
                [min_tau, max_tau_omega // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
            # [[States.HIGH, States.LOW], [min_tau, max_tau_gamma], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.HIGH], [min_tau, max_tau_gamma], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.LOW], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            # [[States.HIGH, States.HIGH], [min_tau, max_tau_gamma//3], num_steps, num_reps, num_runs],
            [
                [States.LOW, States.HIGH],
                [min_tau, max_tau_gamma],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.LOW],
                [min_tau, max_tau_gamma],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.HIGH],
                [min_tau, max_tau_gamma // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
            [
                [States.LOW, States.LOW],
                [min_tau, max_tau_gamma // 3],
                num_steps,
                num_reps,
                num_runs,
            ],
        ],
        dtype=object,
    )

    t1_dq_main.main(
        nv_sig, apd_indices, t1_exp_array, num_runs, composite_pulses=True
    )


def do_ramsey(nv_sig, apd_indices):

    detuning = 2.5  # MHz
    precession_time_range = [0, 4 * 10 ** 3]
    num_steps = 151
    num_reps = 3 * 10 ** 5
    num_runs = 1

    ramsey.main(
        nv_sig,
        apd_indices,
        detuning,
        precession_time_range,
        num_steps,
        num_reps,
        num_runs,
    )


def do_spin_echo(nv_sig, apd_indices):

    # T2* in nanodiamond NVs is just a couple us at 300 K
    # In bulk it's more like 100 us at 300 K
    max_time = 120  # us
    num_steps = max_time  # 1 point per us
    precession_time_range = [1e3, max_time * 10 ** 3]
    num_reps = 4e3
    # num_runs = 5
    num_runs = 20

    #    num_steps = 151
    #    precession_time_range = [0, 10*10**3]
    #    num_reps = int(10.0 * 10**4)
    #    num_runs = 6

    state = States.LOW

    angle = spin_echo.main(
        nv_sig,
        apd_indices,
        precession_time_range,
        num_steps,
        num_reps,
        num_runs,
        state,
    )
    return angle


def do_spin_echo_battery(nv_sig, apd_indices):

    do_pulsed_resonance_state(nv_sig, apd_indices, States.LOW)
    do_pulsed_resonance_state(nv_sig, apd_indices, States.HIGH)
    do_rabi(nv_sig, apd_indices, States.LOW, uwave_time_range=[0, 400])
    do_rabi(nv_sig, apd_indices, States.HIGH, uwave_time_range=[0, 400])
    angle = do_spin_echo(nv_sig, apd_indices)
    return angle



def do_map_rabi_contrast_diff_NIR(nv_sig, apd_indices, state, uwave_time_range, image_range, num_steps, rabi_num_reps, rabi_num_runs):
    
    nir_laser_voltage = 1.3
    
    map_rabi_contrast_NIR.main(nv_sig,
    image_range,
    num_steps,
    apd_indices,
    nir_laser_voltage,
    rabi_num_reps,
    rabi_num_runs,
    state,
    uwave_time_range
    )


def do_nir_battery(nv_sig, apd_indices):

    # do_image_sample(nv_sig, apd_indices)
    # do_pulsed_resonance_state(nv_sig, apd_indices, States.LOW)
    # do_pulsed_resonance_state(nv_sig, apd_indices, States.HIGH)
    # do_rabi(nv_sig, apd_indices, States.LOW, uwave_time_range=[0, 500])
    # do_rabi(nv_sig, apd_indices, States.HIGH, uwave_time_range=[0, 500])
    # do_discrete_rabi(nv_sig, apd_indices, States.LOW, 4)
    # do_discrete_rabi(nv_sig, apd_indices, States.HIGH, 4)
    # do_spin_echo(nv_sig, apd_indices)
    # do_determine_charge_readout_params(nv_sig,apd_indices)
    # do_determine_standard_readout_params(nv_sig,apd_indices)

    with labrad.connect() as cxn:
        power_supply = cxn.power_supply_mp710087
        power_supply.output_on()
        power_supply.set_voltage(1.3)
    time.sleep(1)
    
  
    with labrad.connect() as cxn:
        power_supply = cxn.power_supply_mp710087
        power_supply.output_off()
    time.sleep(1)
        

    # do_image_sample(nv_sig, apd_indices)
    # do_pulsed_resonance_state(nv_sig, apd_indices, States.LOW)
    # do_pulsed_resonance_state(nv_sig, apd_indices, States.HIGH)
    # do_rabi(nv_sig, apd_indices, States.LOW, uwave_time_range=[0, 500])
    # do_rabi(nv_sig, apd_indices, States.HIGH, uwave_time_range=[0, 500])
    # do_discrete_rabi(nv_sig, apd_indices, States.LOW, 4)
    # do_discrete_rabi(nv_sig, apd_indices, States.HIGH, 4)
    # nv_sig["spin_pol_dur"] = 1e6
    # do_t1_dq_knill(nv_sig, apd_indices)
    # do_determine_charge_readout_params(nv_sig,apd_indices)
    # do_determine_standard_readout_params(nv_sig,apd_indices)
    


def do_nir_temp_differential(nv_sig, apd_indices):
    
    dD_dT = -74e-6  # GHz / K

    low_res, low_res_err = do_four_point_esr(nv_sig, apd_indices, States.LOW)
    high_res, high_res_err = do_four_point_esr(nv_sig, apd_indices, States.HIGH)
    zfs = (low_res + high_res) / 2
    zfs_err = np.sqrt(low_res_err**2 + high_res_err**2) / 2
    d_temp = zfs / dD_dT
    d_temp_err = zfs_err / abs(dD_dT)

    with labrad.connect() as cxn:
        power_supply = cxn.power_supply_mp710087
        power_supply.output_on()
        power_supply.set_voltage(1.3)
    time.sleep(1)

    nir_low_res, nir_low_res_err = do_four_point_esr(nv_sig, apd_indices, States.LOW)
    nir_high_res, nir_high_res_err = do_four_point_esr(nv_sig, apd_indices, States.HIGH)
    nir_zfs = (nir_low_res + nir_high_res) / 2
    nir_zfs_err = np.sqrt(nir_low_res_err**2 + nir_high_res_err**2) / 2
    nir_d_temp = nir_zfs / dD_dT
    nir_d_temp_err = nir_zfs_err / abs(dD_dT)

    with labrad.connect() as cxn:
        power_supply = cxn.power_supply_mp710087
        power_supply.output_off()
    time.sleep(1)
    
    print(low_res, low_res_err)
    print(high_res, high_res_err)
    print(zfs, zfs_err)
    print(d_temp, d_temp_err)
    
    print(nir_low_res, nir_low_res_err)
    print(nir_high_res, nir_high_res_err)
    print(nir_zfs, nir_zfs_err)
    print(nir_d_temp, nir_d_temp_err)
    
    print((nir_zfs - zfs) / dD_dT)
    print(np.sqrt(nir_zfs_err**2 + zfs_err**2) / abs(dD_dT))


def do_nir_temp_differential2(nv_sig, apd_indices):
    
    tool_belt.init_safe_stop()
    
    sig_files = []
    ref_files = []
    
    while not tool_belt.safe_stop():

        with labrad.connect() as cxn:
            power_supply = cxn.power_supply_mp710087
            power_supply.output_on()
            power_supply.set_voltage(1.3)
        time.sleep(10)
        
        _, _, low_file = do_four_point_esr(nv_sig, apd_indices, States.LOW)
        _, _, high_file = do_four_point_esr(nv_sig, apd_indices, States.HIGH)
        sig_files.append([low_file, high_file])
    
        with labrad.connect() as cxn:
            power_supply = cxn.power_supply_mp710087
            power_supply.output_off()
        time.sleep(10)
    
        _, _, low_file = do_four_point_esr(nv_sig, apd_indices, States.LOW)
        _, _, high_file = do_four_point_esr(nv_sig, apd_indices, States.HIGH)
        ref_files.append([low_file, high_file])
        
    with labrad.connect() as cxn:
        power_supply = cxn.power_supply_mp710087
        power_supply.output_off()
    
    
    timestamp = tool_belt.get_time_stamp()
    rawData = {
        "timestamp": timestamp,
        "nv_sig": nv_sig,
        "nv_sig-units": tool_belt.get_nv_sig_units(),
        "sig_files": sig_files,
        "ref_files": ref_files,
    }
    nv_name = nv_sig["name"]
    filePath = tool_belt.get_file_path(__file__, timestamp, nv_name)
    tool_belt.save_raw_data(rawData, filePath)        


def do_test_major_routines(nv_sig, apd_indices):
    """Run this whenver you make a significant code change. It'll make sure
    you didn't break anything in the major routines.
    """

    test_major_routines.main(nv_sig, apd_indices)
    

# %% Run the file


if __name__ == "__main__":

    # %% Shared parameters

    # apd_indices = [0]
    # apd_indices = [1]
    apd_indices = [0, 1]

    sample_name = "hopper"

    green_laser = "laserglow_532"
    yellow_laser = "laserglow_589"
    red_laser = "cobolt_638"

    nv_sig = {
        'coords': [0.0, 0.0, 0], 'name': '{}-search'.format(sample_name),
        'disable_opt': True, "disable_z_opt": False, 'expected_count_rate': 1300,

        # 'imaging_laser': green_laser, 'imaging_laser_filter': "nd_0", 'imaging_readout_dur': 1e7,
        # 'imaging_laser': green_laser, 'imaging_laser_filter': "nd_0", 'imaging_readout_dur': 1e8,
        "imaging_laser": green_laser, "imaging_laser_filter": "nd_0.5", "imaging_readout_dur": 1e7,
        # 'imaging_laser': green_laser, 'imaging_laser_filter': "nd_0", 'imaging_readout_dur': 1e9,
        # 'imaging_laser': yellow_laser, 'imaging_laser_power': 1.0, 'imaging_readout_dur': 1e8,
        # 'imaging_laser': red_laser, 'imaging_readout_dur': 1e7,
        # 'spin_laser': green_laser, 'spin_laser_filter': 'nd_0.5', 'spin_pol_dur': 1E6, 'spin_readout_dur': 350,
        # 'spin_laser': green_laser, 'spin_laser_filter': 'nd_0.5', 'spin_pol_dur': 2e3, 'spin_readout_dur': 350,
        "spin_laser": green_laser,
        "spin_laser_filter": "nd_0",
        "spin_pol_dur": 10e3,
        "spin_readout_dur": 350,
        # 'spin_laser': green_laser, 'spin_laser_filter': 'nd_0', 'spin_pol_dur': 1E4, 'spin_readout_dur': 300,
        "nv-_reionization_laser": green_laser,
        "nv-_reionization_dur": 1e6,
        "nv-_reionization_laser_filter": "nd_1.0",
        # 'nv-_reionization_laser': green_laser, 'nv-_reionization_dur': 1E5, 'nv-_reionization_laser_filter': 'nd_0.5',
        "nv-_prep_laser": green_laser,
        "nv-_prep_laser_dur": 1e6,
        "nv-_prep_laser_filter": "nd_0",
        # 'nv-_prep_laser': green_laser, 'nv-_prep_laser_dur': 1E4, 'nv-_prep_laser_filter': 'nd_0.5',
        "nv0_ionization_laser": red_laser,
        "nv0_ionization_dur": 75,
        "nv0_prep_laser": red_laser,
        "nv0_prep_laser_dur": 75,
        "spin_shelf_laser": yellow_laser,
        "spin_shelf_dur": 0,
        "spin_shelf_laser_power": 1.0,
        # 'spin_shelf_laser': green_laser, 'spin_shelf_dur': 50,
        "initialize_laser": green_laser,
        "initialize_dur": 1e4,
        "charge_readout_laser": yellow_laser,
        "charge_readout_dur": 100e6,
        "charge_readout_laser_power": 1.0,
        # "charge_readout_laser": yellow_laser, "charge_readout_dur": 10e6, "charge_readout_laser_power": 1.0,

        'collection_filter': None, 'magnet_angle': None,
        'resonance_LOW': 2.8059, 'rabi_LOW': 226.9, 'uwave_power_LOW': 16.5,
        'resonance_HIGH': 2.9363, 'rabi_HIGH': 300, 'uwave_power_HIGH': 16.5,
        }


    # %% Functions to run

    try:
        
        # pass

        tool_belt.init_safe_stop()

        # Increasing x moves the image down, increasing y moves the image left
        # with labrad.connect() as cxn:
        #     cxn.cryo_piezos.write_xy(150,0) 

        

        # tool_belt.set_drift([0.0, 0.0, 0.0])  # Totally reset
        # drift = tool_belt.get_drift()
        # tool_belt.set_drift([0.0, 0.0, drift[2]])  # Keep z
        # tool_belt.set_drift([drift[0], drift[1], 0.0])  # Keep xy
        # tool_belt.set_xyz(labrad.connect(), [0.0,0.0,0.0])
     
        # do_ensemble_image_sample(nv_sig, apd_indices, scan_range=1.3, num_steps=60, 
        #                               technique='NIR_counts_diff_faster',nv_minus_initialization=False,
        #                                 cbarmin_percdiff=-.0,cbarmax_percdiff=.125,
        #                                 cbarmin_diffcounts=0,cbarmax_diffcounts=220,
        #                                 cbarmin_counts=0,cbarmax_counts=2200
        #                                 )
        
        do_image_sample(nv_sig, apd_indices)
        # do_image_sample_zoom(nv_sig, apd_indices)
        # do_image_sample(nv_sig, apd_indices, nv_minus_initialization=True)
        # do_image_sample_zoom(nv_sig, apd_indices, nv_minus_initialization=True)
        # do_optimize(nv_sig, apd_indices)
        # do_stationary_count(nv_sig, apd_indices, disable_opt=True)
        # do_stationary_count(nv_sig, apd_indices, disable_opt=True, nv_minus_initialization=True)
        # do_stationary_count(nv_sig, apd_indices, disable_opt=True, nv_zero_initialization=True)
        # do_resonance(nv_sig, apd_indices, 2.87, 0.200)
        # do_resonance_state(nv_sig , apd_indices, States.LOW)
        # do_resonance_state(nv_sig, apd_indices, States.HIGH)
        # do_pulsed_resonance(nv_sig, apd_indices, 2.87, 0.200)
        # do_pulsed_resonance_state(nv_sig, apd_indices, States.LOW)
        # do_pulsed_resonance_state(nv_sig, apd_indices, States.HIGH)
        # # # do_scc_resonance(nv_sig, apd_indices, States.LOW)
        # # # do_scc_resonance(nv_sig, apd_indices, States.HIGH)
        # # # do_optimize_magnet_angle(nv_sig, apd_indices)
        # # # do_optimize_magnet_angle_fine(nv_sig, apd_indices)
        # # # do_spin_echo_battery(nv_sig, apd_indices)

        # do_rabi(nv_sig, apd_indices, States.LOW, uwave_time_range=[0, 500])
        # do_rabi(nv_sig, apd_indices, States.HIGH, uwave_time_range=[0, 500])
        # do_discrete_rabi(nv_sig, apd_indices, States.LOW, 4)0203c
        # do_discrete_rabi(nv_sig, apd_indices, States.LOW, 4)
        # do_spin_echo(nv_sig, apd_indices)
        # do_g2_measurement(nv_sig, 0, 1)  # 0, (394.6-206.0)/31 = 6.084 ns, 164.3 MHz; 1, (396.8-203.6)/33 = 5.855 ns, 170.8 MHz
        # do_t1_battery(nv_sig, apd_indices)
        # do_t1_interleave_knill(nv_sig, apd_indices)
        # for i in range(4):
        #     do_t1_dq_knill_battery(nv_sig, apd_indices)
        # do_nir_battery(nv_sig, apd_indices)
        # do_determine_standard_readout_params(nv_sig, apd_indices)

        # do_four_point_esr(nv_sig, apd_indices, States.LOW)
        # do_four_point_esr(nv_sig, apd_indices, States.HIGH)
        # do_nir_temp_differential(nv_sig, apd_indices)
        # do_nir_temp_differential2(nv_sig, apd_indices)
        # do_test_four_point_vs_avg(nv_sig, apd_indices)
        # do_image_sample_temperature(nv_sig, apd_indices)
        # do_ensemble_image_sample(nv_sig, apd_indices,scan_range=.4, num_steps=40, technique='NIR_counts_diff_faster',nv_minus_initialization=False)
        # do_pulsed_resonance(nv_sig, apd_indices, 2.87, 0.200)
        # do_pulsed_resonance_state(nv_sig, apd_indices, States.LOW)
        # do_pulsed_resonance_state(nv_sig, apd_indices, States.HIGH)
        # do_rabi(nv_sig, apd_indices, States.LOW, uwave_time_range=[0, 400])
        # do_rabi(nv_sig, apd_indices, States.HIGH, uwave_time_range=[0, 400])
        # do_spin_echo(nv_sig, apd_indices)

        # SCC characterization
        # do_determine_charge_readout_params(nv_sig,apd_indices,nbins=200,nreps=100)
       
        # do_determine_charge_readout_params_1Dscan_target(nv_sig, apd_indices, 
        #                                                   x_steps=[0.0], y_steps=np.arange(-.55,0.55,.015), 
        #                                                   nbins=200,nreps=500,NIRtest=True)
        
        # do_determine_charge_readout_params_1Dscan_target(nv_sig, apd_indices, 
        #                                                   x_steps=np.arange(-.55,0.55,.015), y_steps=[0.0], 
        #                                                   nbins=200,nreps=500,NIRtest=True)
        
        # do_scc_pulsed_resonance(nv_sig, apd_indices)

        # Automatic T1 setup
        # do_stationary_count(nv_sig, apd_indices)
        # do_pulsed_resonance_state(nv_sig, apd_indices, States.LOW)
        # do_pulsed_resonance_state(nv_sig, apd_indices, States.HIGH)
        # do_rabi(nv_sig, apd_indices, States.LOW, uwave_time_range=[0, 1000])
        # do_rabi(nv_sig, apd_indices, States.HIGH, uwave_time_range=[0, 1000])
        # do_discrete_rabi(nv_sig, apd_indices, States.LOW, 6)
        # do_discrete_rabi(nv_sig, apd_indices, States.HIGH, 6)
        # nv_sig["spin_pol_dur"] = 1e6
        # # # # # do_t1_interleave_knill(nv_sig, apd_indices)
        # # do_t1_dq(nv_sig, apd_indices)
        # do_t1_dq_knill(nv_sig, apd_indices)

    except Exception as exc:
        # recipient = "cdfox@wisc.edu"
        recipient = "cambria@wisc.edu"
        tool_belt.send_exception_email(email_to=recipient)
        raise exc

    finally:

        # Reset our hardware - this should be done in each routine, but
        # let's double check here
        tool_belt.reset_cfm()
        tool_belt.reset_safe_stop()

        # Leave green on
        # with labrad.connect() as cxn:
        #     cxn.pulse_streamer.constant([3], 0.0, 0.0)

