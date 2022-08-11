# -*- coding: utf-8 -*-
"""
Output server for the SRS SG394 microwave signal generator.

Created on Wed Apr 10 12:53:38 2019

@author: mccambria

### BEGIN NODE INFO
[info]
name = signal_generator_sg394
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
import pyvisa as visa  # Docs here: https://pyvisa.readthedocs.io/en/master/
import nidaqmx
import nidaqmx.stream_writers as stream_writers
from nidaqmx.constants import AcquisitionType
import socket
import logging


class SignalGeneratorSg394(LabradServer):
    name = 'signal_generator_sg394'
    pc_name = socket.gethostname()

    def initServer(self):
        filename = 'E:/Shared drives/Kolkowitz Lab Group/nvdata/pc_{}/labrad_logging/{}.log'
        filename = filename.format(self.pc_name, self.name)
        logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%y-%m-%d_%H-%M-%S', filename=filename)
        config = ensureDeferred(self.get_config())
        config.addCallback(self.on_get_config)

    async def get_config(self):
        p = self.client.registry.packet()
        p.cd(['', 'Config', 'DeviceIDs'])
        p.get('signal_generator_sg394_visa_address')
        p.cd(['', 'Config', 'Wiring', 'Daq'])
        p.get('di_clock')  
        result = await p.send()
        return result['get']

    def on_get_config(self, config):
        resource_manager = visa.ResourceManager()
        self.sig_gen = resource_manager.open_resource(config[0])
        # Set the VISA read and write termination. This is specific to the
        # instrument - you can find it in the instrument's programming manual
        self.sig_gen.read_termination = '\r\n'
        self.sig_gen.write_termination = '\r\n'
        # Set our channels for FM
        self.daq_di_pulser_clock = config[1]
        # self.daq_ao_sig_gen_mod = config[2]
        self.task = None    # Initialize state variable
        self.reset(None)
        logging.info('Init complete')

    @setting(0)
    def uwave_on(self, c):
        """Turn on the signal. This is like opening an internal gate on
        the signal generator.
        """

        self.sig_gen.write('ENBR 1')

    @setting(1)
    def uwave_off(self, c):
        """Turn off the signal. This is like closing an internal gate on
        the signal generator.
        """

        self.sig_gen.write('ENBR 0')

    @setting(2, freq='v[]')
    def set_freq(self, c, freq):
        """Set the frequency of the signal.

        Params
            freq: float
                The frequency of the signal in GHz
        """

        # Determine how many decimal places we need
        precision = len(str(freq).split('.')[1])
        self.sig_gen.write('FREQ {0:.{1}f} GHZ'.format(freq, precision))

    @setting(3, amp='v[]')
    def set_amp(self, c, amp):
        """Set the amplitude of the signal.

        Params
            amp: float
                The amplitude of the signal in dBm
        """

        # Determine how many decimal places we need
        precision = len(str(amp).split('.')[1])
        cmd = 'AMPR {0:.{1}f} DBM'.format(amp, precision)
        # logging.info(cmd)
        self.sig_gen.write(cmd)

    # def load_stream_writer(self, task_name, voltages, period):
    #     # Close the existing task and create a new one
    #     if self.task is not None:
    #         self.task.close()
    #     task = nidaqmx.Task(task_name)
    #     self.stream_task = task

    #     # Set up the output channels
    #     task.ao_channels.add_ao_voltage_chan(self.daq_ao_sig_gen_mod,
    #                                          min_val=-1.0, max_val=1.0)

    #     # Set up the output stream
    #     output_stream = nidaqmx.task.OutStream(task)
    #     writer = stream_writers.AnalogMultiChannelWriter(output_stream)

    #     # Configure the sample to advance on the rising edge of the PFI input.
    #     # The frequency specified is just the max expected rate in this case.
    #     # We'll stop once we've run all the samples.
    #     freq = float(1/(period*(10**-9)))  # freq in seconds as a float
    #     task.timing.cfg_samp_clk_timing(freq, source=self.daq_di_pulser_clock,
    #                                     sample_mode=AcquisitionType.CONTINUOUS)

    #     # Start the task before writing so that the channel will sit on
    #     # the last value when the task stops. The first sample won't actually
    #     # be written until the first clock signal.
    #     task.start()

    #     writer.write_many_sample(voltages)

    # @setting(4, fm_range='v[]', voltages='*v[]', period='i')
    # def load_fm(self, c, fm_range, voltages, period):
    #     """Set up frequency modulation via an external voltage. This has never
    #     been used or tested and needs work.
    #     """

    #     # Set up the DAQ AO that will control the modulation
    #     self.load_stream_writer('UwaveSigGen-load_fm', voltages, period)
    #     # Simple FM is type 1, subtype 0
    #     self.sig_gen.write('TYPE 1')
    #     self.sig_gen.write('STYP 0')
    #     # Set the range of the modulation
    #     precision = len(str(fm_range).split('.')[1])
    #     self.sig_gen.write('FDEV {0:.{1}f}GHZ'.format(fm_range, precision))
    #     # Set to an external source
    #     self.sig_gen.write('MFNC 5')
    #     # Turn on FM
    #     self.sig_gen.write('MODL 1')


    @setting(5)
    def mod_off(self, c):
        """Turn off the modulation."""

        self.sig_gen.write('MODL 0')
        task = self.task
        if task is not None:
            task.close()

    @setting(7)
    def load_iq(self, c):
        """
        Set up external IQ modulation
        """

        # QAM is type 7
        self.sig_gen.write('TYPE 7')
        # STYP 1 is vector modulation
        # self.sig_gen.write('STYP 1')
        # External mode is modulation function 5
        self.sig_gen.write('QFNC 5')
        # Turn on modulation
        cmd = 'MODL 1'
        self.sig_gen.write(cmd)
        logging.info(cmd)
        

    @setting(6)
    def reset(self, c):
        self.sig_gen.write('FDEV 0')
        self.uwave_off(c)
        self.mod_off(c)
        # Clean up the DAQ task!
        # if self.task is not None:
        #     crash = 1/0
        # Set the DAQ AO to 0
        # with nidaqmx.Task() as task:
        #     # Set up the output channels
        #     task.ao_channels.add_ao_voltage_chan(self.daq_ao_sig_gen_mod,
        #                                          min_val=-1.0, max_val=1.0)
        #     task.write(0.0)
    
    @setting(10,
             returns="i")
    def querry_err(self, c, option):
        
        err = self.sig_gen.write('LERR?')
        return err
    
    @setting(11,
             returns="i")
    def querry_mode(self, c):
        self.sig_gen.write('TYPE 7')
        self.sig_gen.write('QFNC 1')
        # boo = self.sig_gen.write('QFNC?')
        
        cmd = 'QFNC?'
        logging.info(cmd)
        boo = self.sig_gen.write(cmd)
        logging.info(boo)
        return boo

__server__ = SignalGeneratorSg394()

if __name__ == '__main__':
    from labrad import util
    util.runServer(__server__)
