# -*- coding: utf-8 -*-
"""Template for scripts that should be run directly from the files themselves
(as opposed to from the control panel, for example).

Created on Sun Jun 16 11:22:40 2019

@author: mccambria
"""


# %% Imports

import serial
import time

# %% Constants


# %% Functions


# %% Main


def main(dev):
    """When you run the file, we'll call into main, which should contain the
    body of the script.
    """
    
    # This setup sequence is adapted straight from section 2.1 (USB Interface
    # adapter) of the Elliptec communication protocol manual
    time.sleep(0.1)
    dev.flush()
    time.sleep(0.1)
    # Reset - not supported by pyftdi apparently
    
    print('writing...')
    print(dev.write(b'0in'))
    print('written')
    print(dev.read(33))
    
    print(dev.write(b'0ma00000020'))


# %% Run the file


# The __name__ variable will only be '__main__' if you run this file directly.
# This allows a file's functions, classes, etc to be imported without running
# the script that you set up here.
if __name__ == '__main__':

    # Set up your parameters to be passed to main here

    # Run the script
    try:
        dev = serial.Serial('COM5', 9600, serial.EIGHTBITS,
                            serial.PARITY_NONE, serial.STOPBITS_ONE)
        main(dev)
    except Exception as e:
        print(e)
    finally:
        dev.flush()
        del dev
        