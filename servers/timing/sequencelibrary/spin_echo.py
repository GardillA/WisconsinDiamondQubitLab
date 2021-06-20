# -*- coding: utf-8 -*-
"""
Created on Sat May  4 08:34:08 2019

@author: Aedan
"""

from pulsestreamer import Sequence
from pulsestreamer import OutputState
import numpy
import utils.tool_belt as tool_belt
from utils.tool_belt import States

LOW = 0
HIGH = 1


def get_seq(pulser_wiring, args):

    # %% Parse wiring and args

    # The first 11 args are ns durations and we need them as int64s
    durations = []
    for ind in range(13):
        durations.append(numpy.int64(args[ind]))

    # Unpack the durations
    tau_shrt, polarization_time, signal_time, reference_time,  \
            sig_to_ref_wait_time, pre_uwave_exp_wait_time,  \
            post_uwave_exp_wait_time, aom_delay_time, rf_delay_time,  \
            gate_time, pi_pulse, pi_on_2_pulse, tau_long = durations

    # Get the APD indices
    apd_index = args[13]

    # Signify which signal generator to use
    sig_gen_name = args[14]
    
    # Laser specs
    laser_name = args[15]
    laser_power = args[16]

    pulser_do_apd_gate = pulser_wiring['do_apd_{}_gate'.format(apd_index)]

    sig_gen_gate_chan_name = 'do_{}_gate'.format(sig_gen_name)
    pulser_do_sig_gen_gate = pulser_wiring[sig_gen_gate_chan_name]

    # %% Write the microwave sequence to be used.

    # In t1, the sequence is just a pi pulse, wait for a relaxation time, then
    # then a second pi pulse

    # I define both the time of this experiment, which is useful for the AOM
    # and gate sequences to dictate the time for them to be LOW
    # And I define the actual uwave experiement to be plugged into the rf
    # sequence. I hope that this formatting works.

    # With future protocols--ramsey, spin echo, etc--it will be easy to use
    # this format of sequence building and just change this section of the file

    uwave_experiment_shrt = pi_on_2_pulse + tau_shrt + pi_pulse + \
                            tau_shrt + pi_on_2_pulse
    uwave_experiment_long = pi_on_2_pulse + tau_long + pi_pulse + \
                            tau_long + pi_on_2_pulse

    # %% Couple calculated values

    prep_time = aom_delay_time + rf_delay_time + \
        polarization_time + pre_uwave_exp_wait_time + \
        uwave_experiment_shrt + post_uwave_exp_wait_time

    up_to_long_gates = prep_time + signal_time + sig_to_ref_wait_time + \
        reference_time + pre_uwave_exp_wait_time + \
        uwave_experiment_long + post_uwave_exp_wait_time

    # %% Calclate total period. This is fixed for each tau index

    # The period is independent of the particular tau, but it must be long
    # enough to accomodate the longest tau
    period = aom_delay_time + rf_delay_time + polarization_time + \
        pre_uwave_exp_wait_time + uwave_experiment_shrt + post_uwave_exp_wait_time + \
        signal_time + sig_to_ref_wait_time + reference_time + pre_uwave_exp_wait_time + \
        uwave_experiment_long + post_uwave_exp_wait_time + \
        signal_time + sig_to_ref_wait_time + reference_time

    # %% Define the sequence

    seq = Sequence()

    # APD gating
    pre_duration = prep_time
    short_sig_to_short_ref = signal_time + sig_to_ref_wait_time - gate_time
    short_ref_to_long_sig = up_to_long_gates - (prep_time + signal_time + sig_to_ref_wait_time + gate_time)
    long_sig_to_long_ref = signal_time + sig_to_ref_wait_time - gate_time
    post_duration = reference_time - gate_time
    train = [(pre_duration, LOW),
             (gate_time, HIGH),
             (short_sig_to_short_ref, LOW),
             (gate_time, HIGH),
             (short_ref_to_long_sig, LOW),
             (gate_time, HIGH),
             (long_sig_to_long_ref, LOW),
             (gate_time, HIGH),
             (post_duration, LOW)]
    seq.setDigital(pulser_do_apd_gate, train)

    # Laser
    if laser_power == -1:
        laser_high = HIGH
        laser_low = LOW
    else:
        laser_high = laser_power
        laser_low = 0.0
    train = [(rf_delay_time + polarization_time, laser_high),
             (pre_uwave_exp_wait_time + uwave_experiment_shrt + post_uwave_exp_wait_time, laser_low),
             (signal_time, laser_high),
             (sig_to_ref_wait_time, laser_low),
             (reference_time, laser_high),
             (pre_uwave_exp_wait_time + uwave_experiment_long + post_uwave_exp_wait_time, laser_low),
             (signal_time, laser_high),
             (sig_to_ref_wait_time, laser_low),
             (reference_time + aom_delay_time, laser_high)]
    if laser_power == -1:
        pulser_laser_mod = pulser_wiring['do_{}_dm'.format(laser_name)]
        seq.setDigital(pulser_laser_mod, train)
    else:
        pulser_laser_mod = pulser_wiring['ao_{}_am'.format(laser_name)]
        seq.setAnalog(pulser_laser_mod, train)

    # Pulse the microwave for tau
    pre_duration = aom_delay_time + polarization_time + pre_uwave_exp_wait_time
    mid_duration = post_uwave_exp_wait_time + signal_time + sig_to_ref_wait_time + \
        reference_time + pre_uwave_exp_wait_time
    post_duration = post_uwave_exp_wait_time + signal_time + \
        sig_to_ref_wait_time + reference_time + rf_delay_time

    train = [(pre_duration, LOW)]
    train.extend([(pi_on_2_pulse, HIGH), (tau_shrt, LOW)])
    train.extend([(pi_pulse, HIGH)])
    train.extend([(tau_shrt, LOW), (pi_on_2_pulse, HIGH)])
    train.extend([(mid_duration, LOW)])
    train.extend([(pi_on_2_pulse, HIGH), (tau_long, LOW)])
    train.extend([(pi_pulse, HIGH)])
    train.extend([(tau_long, LOW), (pi_on_2_pulse, HIGH)])
    train.extend([(post_duration, LOW)])
    seq.setDigital(pulser_do_sig_gen_gate, train)

    final_digital = [pulser_wiring['do_sample_clock']]
    final = OutputState(final_digital, 0.0, 0.0)
    return seq, final, [period]

if __name__ == '__main__':
    wiring = {'do_sample_clock': 0,
              'do_apd_0_gate': 4,
              'do_532_aom': 1,
              'do_signal_generator_tsg4104a_gate': 2,
              'do_uwave_gate_1': 3,
              }

            # tau_shrt, polarization_time, signal_time, reference_time
    args = [500, 3000, 3000, 3000,
            # sig_to_ref_wait_time, pre_uwave_exp_wait_time
            2000, 1000,
            # post_uwave_exp_wait_time, aom_delay_time, rf_delay_time
            1000, 0, 0, 
            # gate_time, pi_pulse, pi_on_2_pulse, tau_long
            320, 1000, 500, 2500,
            # apd_index, state_value 
            0, States.LOW.value]
    seq, final, ret_vals = get_seq(wiring, args)
    seq.plot()
