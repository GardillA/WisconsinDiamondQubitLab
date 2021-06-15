# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 16:19:44 2019

@author: mccambria
"""

from pulsestreamer import Sequence
from pulsestreamer import OutputState
import numpy
import utils.tool_belt as tool_belt
from utils.tool_belt import States

LOW = 0
HIGH = 1


def get_seq(pulser_wiring, args):

    # Unpack the args
    readout, am_589_power,  uwave_switch_delay, apd_index, state_value, color_ind = args

    readout = numpy.int64(readout)
    readout = numpy.int64(readout)
    uwave_switch_delay = numpy.int64(uwave_switch_delay)
    wait_time = numpy.int64(300)
    period = readout + wait_time + uwave_switch_delay + readout + wait_time

    # Get what we need out of the wiring dictionary
    pulser_do_daq_clock = pulser_wiring['do_sample_clock']
    pulser_do_apd_gate = pulser_wiring['do_apd_{}_gate'.format(apd_index)]

    if color_ind == 532:
        pulser_do_aom = pulser_wiring['do_532_aom']
    elif color_ind == 638:
        pulser_do_aom = pulser_wiring['do_638_laser']
    elif color_ind == 589:
        pulser_ao_aom = pulser_wiring['ao_589_aom']

    sig_gen_name = tool_belt.get_signal_generator_name(States(state_value))
    sig_gen_gate_chan_name = 'do_{}_gate'.format(sig_gen_name)
    pulser_do_sig_gen_gate = pulser_wiring[sig_gen_gate_chan_name]

    seq = Sequence()

    # Collect two samples
#    train = [(readout + clock_pulse, LOW),
#             (clock_pulse, HIGH),
#             (clock_pulse, LOW),
#             (uwave_switch_delay + readout + clock_pulse, LOW),
#             (clock_pulse, HIGH),
#             (clock_pulse, LOW)]
#    seq.setDigital(pulser_do_daq_clock, train)

    # Ungate the APD channel for the readouts
    train = [(readout, HIGH), (wait_time, LOW),
             (uwave_switch_delay, LOW),
             (readout, HIGH), (wait_time, LOW)]
    seq.setDigital(pulser_do_apd_gate, train)

    # Uwave should be on for the first measurement and off for the second
    train = [(readout, LOW), (wait_time, LOW),
             (uwave_switch_delay, HIGH),
             (readout, HIGH), (wait_time, LOW)]
    seq.setDigital(pulser_do_sig_gen_gate, train)

    # The laser should always be on
    if color_ind == 532 or color_ind == 638:
        train = [(period, HIGH)]
        seq.setDigital(pulser_do_aom, train)
    if color_ind == 589:
        train = [(period, am_589_power)]
        seq.setAnalog(pulser_ao_aom, train)

    final_digital = [pulser_do_daq_clock]
    final = OutputState(final_digital, 0.0, 0.0)
    return seq, final, [period]


if __name__ == '__main__':
    wiring = {'do_sample_clock': 0,
              'do_apd_0_gate': 1,
              'do_532_aom': 2,
              'do_signal_generator_sg394_gate': 7,
              'do_signal_generator_bnc835_gate': 4}
    args = [100000000, 1000000, 0, 1]
    seq, final, ret_vals = get_seq(wiring, args)
    seq.plot()
