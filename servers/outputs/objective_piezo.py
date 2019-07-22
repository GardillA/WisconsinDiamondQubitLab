# -*- coding: utf-8 -*-
"""
Output server for the PI E709 objective piezo.

Created on Thu Apr  4 15:58:30 2019

@author: mccambria

### BEGIN NODE INFO
[info]
name = objective_piezo
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
from pipython import GCSDevice
from pipython import pitools 
import time
import nidaqmx
import logging
import numpy
import nidaqmx.stream_writers as stream_writers


class ObjectivePiezo(LabradServer):
    name = 'objective_piezo'
    logging.basicConfig(level=logging.DEBUG, 
                format='%(asctime)s %(levelname)-8s %(message)s',
                datefmt='%y-%m-%d_%H-%M-%S',
                filename='E:/Shared drives/Kolkowitz Lab Group/nvdata/labrad_logging/{}.log'.format(name))

    def initServer(self):
        config = ensureDeferred(self.get_config())
        config.addCallback(self.on_get_config)

    async def get_config(self):
        p = self.client.registry.packet()
        p.cd('Config')
        p.get('objective_piezo_model')
        p.get('gcs_dll_path')
        p.get('objective_piezo_serial')
        p.cd(['Wiring', 'Daq'])
        p.get('ao_objective_piezo')
        p.get('di_clock')
        result = await p.send()
        return result['get']

    def on_get_config(self, config):
        # Load the generic device
        self.piezo = GCSDevice(devname=config[0], gcsdll=config[1])
        # Connect the specific device with the serial number
        self.piezo.ConnectUSB(config[2])
        # Just one axis for this device
        self.axis = self.piezo.axes[0]
        self.daq_ao_objective_piezo = config[3]
        self.daq_di_clock = config[4]
        logging.debug('Init complete')

    def load_stream_writer(self, c, task_name, voltages, period):

        # Close the existing task if there is one
        if self.task is not None:
            self.close_task_internal()

        # Write the initial voltages and stream the rest
        num_voltages = voltages.shape[1]
        self.write(c, voltages[0, 0], voltages[1, 0])
        stream_voltages = voltages[:, 1:num_voltages]
        stream_voltages = numpy.ascontiguousarray(stream_voltages)
        num_stream_voltages = num_voltages - 1

        # Create a new task
        task = nidaqmx.Task(task_name)
        self.task = task

        # Set up the output channels
        task.ao_channels.add_ao_voltage_chan(self.daq_ao_objective_piezo,
                                             min_val=0.0, max_val=10.0)

        # Set up the output stream
        output_stream = nidaqmx.task.OutStream(task)
        writer = stream_writers.AnalogMultiChannelWriter(output_stream)

        # Configure the sample to advance on the rising edge of the PFI input.
        # The frequency specified is just the max expected rate in this case.
        # We'll stop once we've run all the samples.
        freq = float(1/(period*(10**-9)))  # freq in seconds as a float
        task.timing.cfg_samp_clk_timing(freq, source=self.daq_di_clock,
                                        samps_per_chan=num_stream_voltages)

        writer.write_many_sample(stream_voltages)

        # Close the task once we've written all the samples
        task.register_done_event(self.close_task_internal)

        task.start()

    def close_task_internal(self, task_handle=None, status=None,
                            callback_data=None):
        task = self.task
        if task is not None:
            task.close()
            self.task = None
        return 0
        
    def waitonoma(self, pidevice, axes=None, timeout=300, predelay=0, postdelay=0):
        """This is ripped from the pitools module and adapted for Python 3
        (ie I changed xrange to range...). It hangs until the device completes
        an absolute open-piezo motion.
        """
        axes = pitools.getaxeslist(pidevice, axes)
        numsamples = 5
        positions = []
        maxtime = time.time() + timeout
        pitools.waitonready(pidevice, timeout, predelay)
        while True:
            positions.append(pidevice.qPOS(axes).values())
            positions = positions[-numsamples:]
            if len(positions) < numsamples:
                continue
            isontarget = True
            for vals in zip(*positions):
                isontarget &= 0.01 > sum([abs(vals[i] - vals[i + 1]) for i in range(len(vals) - 1)])
            if isontarget:
                return
            if time.time() > maxtime:
                pitools.stopall(pidevice)
                raise SystemError('waitonoma() timed out after %.1f seconds' % timeout)
            time.sleep(0.01)
        time.sleep(postdelay)

    @setting(2, voltage='v[]')
    def write(self, c, voltage):
        """Write the specified voltages to the piezo"""

        # Close the stream task if it exists
        # This can happen if we quit out early
        if self.task is not None:
            self.close_task_internal()

        with nidaqmx.Task() as task:
            # Set up the output channels
            task.ao_channels.add_ao_voltage_chan(self.daq_ao_objective_piezo,
                                                 min_val=0.0, max_val=10.0)
            task.write(voltage)

    @setting(1, returns='*v[]')
    def read(self, c):
        """Return the current voltages on the x and y channels."""
        with nidaqmx.Task() as task:
            # Set up the internal channels - to do: actual parsing...
            if self.daq_ao_objective_piezo == 'dev1/AO2':
                chan_name = 'dev1/_ao2_vs_aognd'
            task.ai_channels.add_ai_voltage_chan(chan_name,
                                                 min_val=0.0, max_val=10.0)
            voltage = task.read()
        return voltage
    
    @setting(3, center='v[]', scan_range='v[]',
             num_steps='i', period='i', returns='*v[]')
    def load_z_scan(self, c, center, scan_range, num_steps, period):

        half_scan_range = scan_range / 2
        low = center - half_scan_range
        high = center + half_scan_range
        voltages = numpy.linspace(low, high, num_steps)
        self.load_stream_writer(c, 'ObjectivePiezo-load_z_scan', voltages, period)
        return voltages

    @setting(0, voltage='v[]')
    def write_voltage(self, c, voltage):
        """Write a voltage to the piezo.

        Params:
            voltage: float
                The voltage to write
        """

        self.piezo.SVO(self.axis, False)  # Turn off the feedback servo
        self.piezo.SVA(self.axis, voltage)  # Write the voltage
        # Wait until the device completes the motion
        logging.debug('Began motion at {}'.format(time.time()))
#        self.waitonoma(self.piezo, timeout=5)
        try:
            self.waitonoma(self.piezo, timeout=5)
        except Exception as e:
            logging.error(e)
        logging.debug('Completed motion at {}'.format(time.time()))

    @setting(1, returns='v[]')
    def read_position(self, c):
        """Read the position of the piezo.

        Returns
            float
                The current position of the piezo in microns
        """

        return self.piezo.qPOS()[self.axis]


__server__ = ObjectivePiezo()

if __name__ == '__main__':
    from labrad import util
    util.runServer(__server__)
