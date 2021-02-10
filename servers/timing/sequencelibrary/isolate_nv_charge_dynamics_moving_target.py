# -*- coding: utf-8 -*-
"""
Created on Sat Mar  24 08:34:08 2020

Thsi file is for use with the isolate_nv_charge_dynamics_moving_target' routine.

This sequence has three pulses, seperated by wait times that allow time for
the galvo to move. We also have two clock pulses instructing the galvo to move, 
followed by a clock pulse at the end of the sequence to signifiy the counts to read.

@author: Aedan
"""

from pulsestreamer import Sequence
from pulsestreamer import OutputState
import numpy

LOW = 0
HIGH = 1

def get_seq(pulser_wiring, args):

    # %% Parse wiring and args

    # The first 2 args are ns durations and we need them as int64s
    durations = []
    for ind in range(7):
        durations.append(numpy.int64(args[ind]))

    # Unpack the durations
    initialization_time, pulse_time, readout_time, \
        delay_532, delay_589, delay_638, \
        galvo_time = durations
                
    aom_ao_589_pwr = args[7]

    # Get the APD index
    apd_index = args[8]
    
    init_color = args[9]
    pulse_color = args[10]
    read_color = args[11]
    ungate_moving_pulse = args[12]

    # Get what we need out of the wiring dictionary
    pulser_do_apd_gate = pulser_wiring['do_apd_{}_gate'.format(apd_index)]
    pulser_do_clock = pulser_wiring['do_sample_clock']
    pulser_do_532_aom = pulser_wiring['do_532_aom']
    pulser_ao_589_aom = pulser_wiring['ao_589_aom']
    pulser_do_638_aom = pulser_wiring['do_638_laser']
    
    total_laser_delay = delay_532 + delay_589 + delay_638

    # %% Calclate total period.
    
#    We're going to readout the entire sequence, including the clock pulse
    period = total_laser_delay + initialization_time + pulse_time + readout_time \
        + 2 * galvo_time + 3* 100
    
    # %% Define the sequence

    seq = Sequence()
    
    # APD 
#    train = [(period - readout_time - 100, LOW), (readout_time, HIGH), (100, LOW)]
#    train = [(readout_time, HIGH), (100, LOW)]

    
    if ungate_moving_pulse == True:
        train = [(total_laser_delay, LOW), (initialization_time, LOW),
             (100 + galvo_time, LOW), (pulse_time, LOW),
             (100 + galvo_time, LOW), (readout_time, HIGH),
             (100, LOW)]
    else:
        train = [(total_laser_delay, LOW), (initialization_time, HIGH),
             (100 + galvo_time, LOW), (pulse_time, HIGH),
             (100 + galvo_time, LOW), (readout_time, HIGH),
             (100, LOW)]
    seq.setDigital(pulser_do_apd_gate, train)
    
    # clock 
    # I needed to add 100 ns between the redout and the clock pulse, otherwise 
    # the tagger misses some of the gate open/close clicks
    train = [(total_laser_delay + initialization_time+ 100, LOW),(100, HIGH),
             (galvo_time + pulse_time, LOW), (100, HIGH), 
             (galvo_time + readout_time, LOW), (100, HIGH),
             (100, LOW)] 
#    train = [(period + 100, LOW), (100, HIGH), (100, LOW)]
    seq.setDigital(pulser_do_clock, train)
    
#    train = [(period, HIGH)]
#    seq.setDigital(pulser_do_532_aom, train)
    
    # start each laser sequence
    train_532 = [(total_laser_delay - delay_532 , LOW)]
    train_589 = [(total_laser_delay - delay_589, LOW)]
    train_638 = [(total_laser_delay - delay_638, LOW)]
   
    galvo_delay_train = [(100 + galvo_time, LOW)]
    
    # add the initialization pulse segment
    init_train_on = [(initialization_time, HIGH)]
    init_train_off = [(initialization_time, LOW)]
    if init_color == 532:
        train_532.extend(init_train_on)
        train_589.extend(init_train_off)
        train_638.extend(init_train_off)
    if init_color == 589:
        init_train_on = [(initialization_time, aom_ao_589_pwr)]
        train_532.extend(init_train_off)
        train_589.extend(init_train_on)
        train_638.extend(init_train_off)
    if init_color == 638:
        train_532.extend(init_train_off)
        train_589.extend(init_train_off)
        train_638.extend(init_train_on)
        
    train_532.extend(galvo_delay_train)
    train_589.extend(galvo_delay_train)
    train_638.extend(galvo_delay_train)
    
    # add the pulse pulse segment
    pulse_train_on = [(pulse_time, HIGH)]
    pulse_train_off = [(pulse_time, LOW)]
    if pulse_color == 532:
        train_532.extend(pulse_train_on)
        train_589.extend(pulse_train_off)
        train_638.extend(pulse_train_off)
    if pulse_color == 589:
        pulse_train_on = [(pulse_time, aom_ao_589_pwr)]
        train_532.extend(pulse_train_off)
        train_589.extend(pulse_train_on)
        train_638.extend(pulse_train_off)
    if pulse_color == 638:
        train_532.extend(pulse_train_off)
        train_589.extend(pulse_train_off)
        train_638.extend(pulse_train_on)
        
    train_532.extend(galvo_delay_train)
    train_589.extend(galvo_delay_train)
    train_638.extend(galvo_delay_train)
    
    # add the readout pulse segment
    read_train_on = [(readout_time, HIGH)]
    read_train_off = [(readout_time, LOW)]
    if read_color == 532:
        train_532.extend(read_train_on)
        train_589.extend(read_train_off)
        train_638.extend(read_train_off)
    if read_color == 589:
        read_train_on = [(readout_time, aom_ao_589_pwr)]
        train_532.extend(read_train_off)
        train_589.extend(read_train_on)
        train_638.extend(read_train_off)
    if read_color == 638:
        train_532.extend(pulse_train_off)
        train_589.extend(pulse_train_off)
        train_638.extend(read_train_on)
        
    train_532.extend([(100, LOW)])
    train_589.extend([(100, LOW)])
    train_638.extend([(100, LOW)])
        
        

    seq.setDigital(pulser_do_532_aom, train_532)
    seq.setAnalog(pulser_ao_589_aom, train_589)
    seq.setDigital(pulser_do_638_aom, train_638)    
        
    final_digital = []
    final = OutputState(final_digital, 0.0, 0.0)
    return seq, final, [period]

if __name__ == '__main__':
    wiring = {'ao_638_laser': 0,
              'ao_589_aom':1,
              'do_sample_clock': 0,
              'do_apd_0_gate': 4,
              'do_532_aom': 1,
              'do_638_laser': 7
              }

    seq_args = [1000, 1500, 3000, 0, 0, 0, 500, 0.5, 0, 532, 532, 638]
#    seq_args = [1000, 100000, 200000, 140, 1080, 90, 2000000, 0.7, 0, 638, 532, 589]

    seq, final, ret_vals = get_seq(wiring, seq_args)
    seq.plot()
