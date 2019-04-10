# -*- coding: utf-8 -*-
"""
Output server for the Thorlabs GVS212 galvanometer. Controlled by the DAQ.

Created on Mon Apr  8 19:50:12 2019

@author: mccambria

### BEGIN NODE INFO
[info]
name = Galvo
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
import nidaqmx
import nidaqmx.stream_writers as stream_writers
from nidaqmx.constants import AcquisitionType
import numpy


class Galvo(LabradServer):
    name = 'Galvo'

    def initServer(self):
        config = ensureDeferred(self.get_config())
        config.addCallback(self.on_get_config)

    async def get_config(self):
        p = self.client.registry.packet()
        p.cd('Config')
        p.get('daq0_name')
        p.cd('Wiring')
        p.get('daq_ao_galvo_x')
        p.get('daq_ao_galvo_y')
        p.get('daq_di_pulser_clock')
        result = await p.send()
        return result['get']

    def on_get_config(self, config):
        self.daq_name = config[0]
        self.daq_ao_galvo_x = config[1]
        self.daq_ao_galvo_y = config[2]
        self.daq_di_pulser_clock = config[3]

    @setting(0, xVoltage='v[]', yVoltage='v[]')
    def write(self, c, xVoltage, yVoltage):
        with nidaqmx.Task() as task:
            # Set up the output channels
            chanName = self.daq_name + '/AO' + self.daq_ao_galvo_x
            task.ao_channels.add_ao_voltage_chan(chanName,
                                                 min_val=-10.0, max_val=10.0)
            chanName = self.daq_name + '/AO' + self.daq_ao_galvo_y
            task.ao_channels.add_ao_voltage_chan(chanName,
                                                 min_val=-10.0, max_val=10.0)
            task.write([xVoltage, yVoltage])

    @setting(1, returns='*2v[]')
    def read(self, c):
        with nidaqmx.Task() as task:
            # Set up the internal channels
            chan_name = self.daq_name + '/_ao' + \
                self.daq_ao_galvo_x + '_vs_aognd'
            task.ai_channels.add_ai_voltage_chan(chan_name,
                                                 min_val=-10.0, max_val=10.0)
            chan_name = self.daq_name + '/_ao' + \
                self.daq_ao_galvo_y + '_vs_aognd'
            task.ai_channels.add_ai_voltage_chan(chan_name,
                                                 min_val=-10.0, max_val=10.0)
            voltages = task.read()

        return voltages

    def load_stream_writer(self, task_name, voltages, period):
        # Close the existing task and create a new one
        if self.task is not None:
            self.task.close()
        task = nidaqmx.Task(task_name)
        self.task = task

        # Clear other existing stream state attributes
        self.writer = None
        self.stream_voltages = None
        self.buffer_pos = None

        # Set up the output channels
        chan_name = self.daq_name + '/AO' + self.daq_ao_galvo_x
        task.ao_channels.add_ao_voltage_chan(chan_name,
                                             min_val=-10.0, max_val=10.0)
        chan_name = self.daq_name + '/AO' + self.daq_ao_galvo_y
        task.ao_channels.add_ao_voltage_chan(chan_name,
                                             min_val=-10.0, max_val=10.0)

        # Set up the output stream
        output_stream = nidaqmx.task.OutStream(task)
        writer = stream_writers.AnalogMultiChannelWriter(output_stream)

        # Configure the sample to advance on the rising edge of the PFI input.
        # The frequency specified is just the max expected rate in this case.
        # We'll stop once we've run all the samples.
        chan_name = 'PFI' + self.daq_di_pulser_clock
        freq = float(1/(period*(10**-9)))  # freq in seconds as a float
        task.timing.cfg_samp_clk_timing(freq, source=chan_name,
                                        sample_mode=AcquisitionType.CONTINUOUS)

        # Start the task before writing so that the channel will sit on
        # the last value when the task stops. The first sample won't actually
        # be written until the first clock signal.
        task.start()

        # We'll write incrementally if there are more than 4000 samples
        # per channel since the DAQ buffer supports 8191 samples max
        if voltages.shape[1] > 4000:
            # Refill the buffer every 3000 samples
            task.register_every_n_samples_transferred_from_buffer_event(3000, self.fill_buffer)
            buffer_voltages = voltages[:, 0:4000]
            # Set up the stream state attributes
            self.stream_writer = writer
            self.stream_voltages = voltages
            self.stream_buffer_pos = 4000
        else:
            buffer_voltages = voltages
        writer.write_many_sample(buffer_voltages)

    def fill_buffer(self):
        # Check if there are more than 3000 samples left to write
        voltages = self.stream_voltages
        buffer_pos = self.stream_buffer_pos
        if voltages.shape[1] - buffer_pos > 3000:
            next_buffer_pos = buffer_pos + 3000
            buffer_voltages = voltages[:, buffer_pos:next_buffer_pos]
            self.stream_buffer_pos = next_buffer_pos
        else:
            buffer_voltages = voltages[:, buffer_pos:]
        cont_buffer_voltages = numpy.ascontiguousarray(buffer_voltages)
        self.stream_writer.write_many_sample(cont_buffer_voltages)

    @setting(2, x_center='v[]', y_center='v[]',
             x_range='v[]', y_range='v[]', num_steps='i', period='i',
             returns='iv[]v[]iv[]v[]v[]')
    def load_scan(self, c, x_center, y_center,
                  x_range, y_range, num_steps, period):

        # Force the scan to have square pixels by only applying num_steps
        # to the shorter axis
        if x_range <= y_range:
            pixel_size = x_range / num_steps
            x_num_steps = num_steps
            y_num_steps = y_range // pixel_size
            y_range = pixel_size * y_num_steps
        else:
            pixel_size = y_range / num_steps
            y_num_steps = num_steps
            x_num_steps = x_range // pixel_size
            x_range = pixel_size * x_num_steps

        # Calculate x and y offsets
        x_offset = x_center - (x_range / 2)
        y_offset = y_center - (y_range / 2)

        # Set up vectors for the number of samples in each direction
        # [0, 1, 2, ... length - 1]
        x_steps = numpy.arange(x_num_steps)
        y_steps = numpy.arange(y_num_steps)

        # Apply scale and offset to get the voltages we'll apply to the galvo
        # Note that the polar/azimuthal angles, not the actual x/y positions
        # are linear in these voltages. For a small range, however, we don't
        # really care.
        x_voltages_1d = (pixel_size * x_steps) + x_offset
        y_voltages_1d = (pixel_size * y_steps) + y_offset

        # Winding cartesian product
        # The x values are repeated and the y values are mirrored and tiled
        # The comments below shows what happens for [1, 2, 3], [4, 5, 6]

        # [1, 2, 3] => [1, 2, 3, 3, 2, 1]
        x_inter = numpy.concatenate((x_voltages_1d,
                                     numpy.flipud(x_voltages_1d)))
        # [1, 2, 3, 3, 2, 1] => [1, 2, 3, 3, 2, 1, 1, 2, 3]
        if y_num_steps % 2 == 0:  # Even x size
            x_voltages = numpy.tile(x_inter, int(y_num_steps/2))
        else:  # Odd x size
            x_voltages = numpy.tile(x_inter, int(numpy.floor(y_num_steps/2)))
            x_voltages = numpy.concatenate((x_voltages, x_voltages_1d))

        # [4, 5, 6] => [4, 4, 4, 5, 5, 5, 6, 6, 6]
        y_voltages = numpy.repeat(y_voltages_1d, x_num_steps)

        voltages = numpy.vstack((x_voltages, y_voltages))
        self.load_stream_writer('Galvo-set_up_sweep', voltages, period)

        x_low = x_voltages_1d[0]
        x_high = x_voltages_1d[len(x_voltages_1d) - 1]
        y_low = y_voltages_1d[0]
        y_high = y_voltages_1d[len(y_voltages_1d) - 1]

        return x_num_steps, x_low, x_high, y_num_steps, y_low, y_high, pixel_size

    @setting(3)
    def close_task(self, c):
        task = self.task
        if task is not None:
            task.close()


__server__ = Galvo()

if __name__ == '__main__':
    from labrad import util
    util.runServer(__server__)
