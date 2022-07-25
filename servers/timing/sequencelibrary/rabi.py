# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:39:27 2019

@author: mccambria
"""

from pulsestreamer import Sequence
from pulsestreamer import OutputState
import numpy
import utils.tool_belt as tool_belt
from utils.tool_belt import States

LOW = 0
HIGH = 1


def get_seq(pulse_streamer, config, args):

    # %% Parse wiring and args

    # The first 9 args are ns durations and we need them as int64s
    durations = []
    for ind in range(4):
        durations.append(numpy.int64(args[ind]))
        
    # Unpack the durations
    tau, polarization_time, gate_time, max_tau = durations

    # Get the APD indices
    apd_index = args[4]

    # Signify which signal generator to use
    state = args[5]
    state = States(state)
    sig_gen_name = config['Microwaves']['sig_gen_{}'.format(state.name)]
    
    # Laser specs
    laser_name = args[6]
    laser_power = args[7]

    # Get what we need out of the wiring dictionary
    pulser_wiring = config['Wiring']['PulseStreamer']
    key = 'do_apd_{}_gate'.format(apd_index)
    pulser_do_apd_gate = pulser_wiring[key]
    sig_gen_gate_chan_name = 'do_{}_gate'.format(sig_gen_name)
    pulser_do_sig_gen_gate = pulser_wiring[sig_gen_gate_chan_name]

    # %% Couple calculated values
    
    aom_delay_time = config['Optics'][laser_name]['delay']
    uwave_delay_time = config['Microwaves'][sig_gen_name]['delay']
    signal_wait_time = config['CommonDurations']['uwave_buffer']
    background_wait_time = signal_wait_time
    reference_wait_time = 2 * signal_wait_time
    reference_time = signal_wait_time

    prep_time = polarization_time + signal_wait_time + \
        tau + signal_wait_time
    end_rest_time = max_tau - tau

    # The period is independent of the particular tau, but it must be long
    # enough to accomodate the longest tau
    # period = aom_delay_time + polarization_time + reference_wait_time + \
    #     reference_wait_time + polarization_time + reference_wait_time + \
    #     reference_time + max_tau
    # print(period)

    # %% Define the sequence

    seq = Sequence()

    # APD gating - first high is for signal, second high is for reference
    pre_duration = aom_delay_time + prep_time
#    mid_duration = period - (pre_duration + (2 * gate_time) + post_duration)
    mid_duration = polarization_time + reference_wait_time - gate_time
    train = [(pre_duration, LOW),
             (gate_time, HIGH),
             (mid_duration, LOW),
             (gate_time, HIGH), 
             (0, LOW)]
    seq.setDigital(pulser_do_apd_gate, train)
    period = 0
    for el in train:
        period += el[0]
    # print(period)

    # # Ungate (high) the APD channel for the background
    # gateBackgroundTrain = [( AOMDelay + preparationTime + polarizationTime + referenceWaitTime + referenceTime + backgroundWaitTime, low),
    #                       (gateTime, high), (endRestTime - gateTime, low)]
    # pulserSequence.setDigital(pulserDODaqGate0, gateBackgroundTrain)

    # Pulse the laser with the AOM for polarization and readout
    train = [(polarization_time, HIGH),
             (signal_wait_time + tau + signal_wait_time, LOW),
             (polarization_time, HIGH),
             (reference_wait_time, LOW),
             (gate_time, HIGH),
             (aom_delay_time, LOW)]
    tool_belt.process_laser_seq(pulse_streamer, seq, config,
                                laser_name, laser_power, train)
    # total_dur = 0
    # for el in train:
    #     total_dur += el[0]
    # print(total_dur)

    # Pulse the microwave for tau
    pre_duration = aom_delay_time + polarization_time + signal_wait_time - uwave_delay_time
    post_duration = signal_wait_time + polarization_time + \
        reference_wait_time + gate_time + uwave_delay_time
    train = [(pre_duration, LOW), (tau, HIGH), (post_duration, LOW)]
    seq.setDigital(pulser_do_sig_gen_gate, train)
    # total_dur = 0
    # for el in train:
    #     total_dur += el[0]
    # print(total_dur)

    final_digital = [pulser_wiring['do_sample_clock']]
    final = OutputState(final_digital, 0.0, 0.0)
    return seq, final, [period]


if __name__ == '__main__':
    config = tool_belt.get_config_dict()
    # tool_belt.set_delays_to_zero(config)
    args = [10, 1000.0, 350, 200, 1, 1, 'integrated_520', None]
    seq = get_seq(None, config, args)[0]
    seq.plot()
