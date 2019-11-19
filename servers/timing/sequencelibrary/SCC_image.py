# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 21:24:36 2019

@author: mccambria
"""

from pulsestreamer import Sequence
from pulsestreamer import OutputState
import numpy

LOW = 0
HIGH = 1


def get_seq(pulser_wiring, args):

    # Unpack the args
    delay, readout, illumination, reionization, power, apd_index = args

    # Get what we need out of the wiring dictionary
    pulser_do_daq_clock = pulser_wiring['do_sample_clock']
    pulser_do_daq_gate = pulser_wiring['do_apd_{}_gate'.format(apd_index)]
    pulser_do_aom = pulser_wiring['do_532_aom']
    pulser_ao_aom = pulser_wiring['ao_589_aom']

    # Convert the 32 bit ints into 64 bit ints
    delay = numpy.int64(delay)
    readout = numpy.int64(readout)
    illumination = numpy.int64(illumination)
    reionization = numpy.int64(reionization)
    period = numpy.int64(delay + illumination + reionization + readout + 1300)

    seq = Sequence()

    # The clock signal will be high for 100 ns with buffers of 100 ns on
    # either side. During the buffers, everything should be low. The buffers
    # account for any timing jitters/delays and ensure that everything we
    # expect to be on one side of the clock signal is indeed on that side.
    train = [(delay + reionization + 1000 + readout + 100, LOW), (100, HIGH), (100, LOW)]
    seq.setDigital(pulser_do_daq_clock, train)

    train = [(delay + reionization + 100, LOW), (readout, HIGH), (300, LOW)]
    seq.setDigital(pulser_do_daq_gate, train)
    
    # pulse the green laser to reionize to NV-
    train = [(reionization, HIGH), (1000 + illumination, LOW)]
    seq.setDigital(pulser_do_aom, train)
    
    # pulse the yellow laser to readout the counts
    train = [(reionization + 1000, 0.0), (illumination, power), (100, 0.0)]
    seq.setAnalog(pulser_ao_aom, train)  

    
    final_digital = [3]
    final = OutputState(final_digital, 0.0, 0.0)


    return seq, final, [period]


if __name__ == '__main__':
    wiring = {'do_sample_clock': 0,
              'do_apd_0_gate': 1,
              'do_532_aom': 2,
              'ao_589_aom': 0}
    args = [0, 500, 500, 250, 1.0, 0]
    seq, ret_vals, _ = get_seq(wiring, args)
    seq.plot()
