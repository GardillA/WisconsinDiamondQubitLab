# -*- coding: utf-8 -*-
"""
Input server for the LASER COMPONENTS COUNT-100C APD. Communicates via the DAQ.

Created on Tue Apr  9 08:52:34 2019

@author: mccambria

### BEGIN NODE INFO
[info]
name = Apd
version = 1.0
description =

[startup]
cmdline = %PYTHON% %FILE%
timeout = 20

[shutdown]
message = 987654321
timeout = 5
### END NODE INFO
"""

from labrad.server import LabradServer
from labrad.server import setting
from twisted.internet.defer import ensureDeferred
import numpy
import nidaqmx
import nidaqmx.stream_readers as stream_readers
from nidaqmx.constants import TriggerType
from nidaqmx.constants import Level
from nidaqmx.constants import AcquisitionType


class Apd(LabradServer):
    name = 'Apd'

    def initServer(self):
        config = ensureDeferred(self.get_config())
        config.addCallback(self.on_get_config)
        self.tasks = []
        self.stream_reader_state = []

    async def get_config(self):
        p = self.client.registry.packet()
        p.cd('Config')
        p.get('daq0_name')
        p.cd('Wiring')
        p.get('daq_di_pulser_clock')
        # The DAQ supports 4 counts max
        for index in range(4):
            p.get('daq_ctr_apd{}'.format(index))
            p.get('daq_ci_apd{}'.format(index))
            p.get('daq_di_pulser_apdgate{}'.format(index))
        result = await p.send()
        return result['get']

    def on_get_config(self, config):
        # The counters share a clock, but everything else is distinct
        self.dev_name = config[0]
        self.daq_di_pulser_clock = config[1]
        self.daq_ctr_apd = []
        self.daq_ci_apd = []
        self.daq_di_pulser_apdgate = []
        # Determine how many counters to set up
        # We assume that all elements past the first empty are also empty
        try:
            first_empty = config.index('')
            counter_config = config[1: first_empty]
        except ValueError:
            # If there are no empties then just take the rest of config
            counter_config = config[1: len(config)]
        num_counters = len(counter_config) / 3
        if not num_counters.is_integer():
            raise ValueError('Number of counters in registry '
                             'is not an integer')
        # Loop through the possible counters
        for index in range(int(num_counters)):
            config_index = 3 * index
            self.daq_ctr_apd[index] = config[config_index+1]
            self.daq_ci_apd[index] = config[config_index+2]
            self.daq_di_pulser_apdgate[index] = config[config_index+3]

    @setting(0, apd_index='i', period='i', total_num_to_read='i')
    def load_stream_reader(self, c, apd_index, period, total_num_to_read):

        task = nidaqmx.Task('Apd-load_stream_reader_' + str(apd_index))
        self.tasks[apd_index] = task

        chan_name = self.dev_name + '/ctr' + self.daq_ctr_apd[apd_index]
        chan = task.ci_channels.add_ci_count_edges_chan(chan_name)
        chan.ci_count_edges_term = 'PFI' + self.daq_ci_apd[apd_index]

        # Set up the input stream
        input_stream = nidaqmx.task.InStream(task)
        reader = stream_readers.CounterReader(input_stream)
        # Just collect whatever data is available when we read
        reader.verify_array_shape = False

        # Set up the gate ('pause trigger')
        # Pause when low - i.e. read only when high
        task.triggers.pause_trigger.trig_type = TriggerType.DIGITAL_LEVEL
        task.triggers.pause_trigger.dig_lvl_when = Level.LOW
        gate_chan_name = 'PFI' + self.daq_di_pulser_apdgate[apd_index]
        task.triggers.pause_trigger.dig_lvl_src = gate_chan_name

        # Configure the sample to advance on the rising edge of the PFI input.
        # The frequency specified is just the max expected rate in this case.
        # We'll stop once we've run all the samples.
        clock_chan_name = 'PFI' + self.daq_di_pulser_clock
        freq = float(1/(period*(10**-9)))  # freq in seconds as a float
        task.timing.cfg_samp_clk_timing(freq, source=clock_chan_name,
                                        sample_mode=AcquisitionType.CONTINUOUS)

        # Initialize the state dictionary for this stream
        self.stream_reader_state[apd_index] = {}
        state_dict = self.stream_reader_state[apd_index]
        state_dict['reader'] = reader
        state_dict['num_read_so_far'] = 0
        state_dict['total_num_to_read'] = total_num_to_read
        # Something funny is happening if we get more
        # than 1000 samples in one read
        state_dict['buffer_size'] = min(total_num_to_read, 1000)
        state_dict['last_value'] = 0  # Last cumulative value we read

        # Start the task. It will start counting immediately so we'll have to
        # discard the first sample.
        task.start()

    @setting(1, apd_index='i', returns='*i')
    def read_stream(self, c, apd_index):

        # Unpack the state dictionary
        state_dict = self.stream_reader_state[apd_index]
        reader = state_dict['reader']
        num_read_so_far = state_dict['num_read_so_far']
        total_num_to_read = state_dict['total_num_to_read']
        buffer_size = state_dict['buffer_size']

        # The counter task begins counting as soon as the task starts.
        # The AO channel writes its first samples only on the first clock
        # signal after the task starts. This means there's one
        # sample from the counter stream that we don't want to record.
        # We do need it for a calculation below, however.
        if num_read_so_far == 0:
            state_dict['last_value'] = reader.read_one_sample_uint32()

        # Initialize the read sample array to its maximum possible size.
        new_samples_cum = numpy.zeros(buffer_size,
                                      dtype=numpy.uint32)

        # Read the samples currently in the DAQ memory.
        num_new_samples = reader.read_many_sample_uint32(new_samples_cum)
        if num_new_samples >= buffer_size:
            raise Warning('The DAQ buffer contained more samples than '
                          'expected. Validate your parameters and '
                          'increase bufferSize if necessary.')

        # Check if we collected more samples than we need, which may happen
        # if the pulser runs longer than necessary. If so, just to throw out
        # excess samples.
        if num_read_so_far + num_new_samples > total_num_to_read:
            num_new_samples = total_num_to_read - num_read_so_far
        new_samples_cum = new_samples_cum[0: num_new_samples]

        # The DAQ counter reader returns cumulative counts, which is not what
        # we want. So we have to calculate the difference between samples
        # n and n-1 in order to get the actual count for the nth sample.
        new_samples_diff = numpy.zeros(num_new_samples)
        for index in range(num_new_samples):
            if index == 0:
                last_value = state_dict['last_value']
            else:
                last_value = new_samples_cum[index-1]

            new_samples_diff[index] = new_samples_cum[index] - last_value

        state_dict['last_value'] = new_samples_cum[num_new_samples-1]

        # Update the current count
        state_dict['num_read_so_far'] = num_read_so_far + num_new_samples

        return new_samples_diff

    @setting(2, apd_index='i')
    def close_task(self, c, apd_index):
        try:
            task = self.tasks[apd_index]
            task.close()
        except Exception:
            pass


__server__ = Apd()

if __name__ == '__main__':
    from labrad import util
    util.runServer(__server__)
