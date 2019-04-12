# -*- coding: utf-8 -*-
"""
Scan the galvos over the designated area, collecting counts at each point.
Generate an image of the sample.

Created on Tue Apr  9 15:18:53 2019

@author: Matt
"""

import numpy
import utils.tool_belt as tool_belt
import time
from twisted.logger import Logger
log = Logger()


def populate_img_array(valsToAdd, imgArray, writePos):
    """
    We scan the sample in a winding pattern. This function takes a chunk
    of the 1D list returned by this process and places each value appropriately
    in the 2D image array. This allows for real time imaging of the sample's
    fluorescence.

    Note that this function could probably be much faster. At least in this
    context, we don't care if it's fast. The implementation below was
    written for simplicity.

    Params:
        valsToAdd: numpy.ndarray
            The increment of raw data to add to the image array
        imgArray: numpy.ndarray
            The xDim x yDim array of fluorescence counts
        writePos: tuple(int)
            The last x, y write position on the image array. [] will default
            to the bottom right corner.
        startingPos: SweepStartingPosition
            Sweep starting position of the winding pattern

    Returns:
        numpy.ndarray: The updated imgArray
        tuple(int): The last x, y write position on the image array
    """

    yDim = imgArray.shape[0]
    xDim = imgArray.shape[1]

    if len(writePos) == 0:
        writePos[:] = [xDim, yDim - 1]

    xPos = writePos[0]
    yPos = writePos[1]

    # Figure out what direction we're heading
    headingLeft = ((yDim - 1 - yPos) % 2 == 0)

    for val in valsToAdd:
        if headingLeft:
            # Determine if we're at the left x edge
            if (xPos == 0):
                yPos = yPos - 1
                imgArray[yPos, xPos] = val
                headingLeft = not headingLeft  # Flip directions
            else:
                xPos = xPos - 1
                imgArray[yPos, xPos] = val
        else:
            # Determine if we're at the right x edge
            if (xPos == xDim - 1):
                yPos = yPos - 1
                imgArray[yPos, xPos] = val
                headingLeft = not headingLeft  # Flip directions
            else:
                xPos = xPos + 1
                imgArray[yPos, xPos] = val
    writePos[:] = [xPos, yPos]


def on_click_image(event):
    """
    Click handler for images. Prints the click coordinates to the console.

    Params:
        event: dictionary
            Dictionary containing event details
    """

    try:
        print('{:.3f}, {:.3f}'.format(event.xdata, event.ydata))
    except TypeError:
        # Ignore TypeError if you click in the figure but out of the image
        pass


def main(cxn, name, coords, x_range, y_range,
         num_steps, readout, apd_index, continuous=False):

    # %% Some initial calculations

    x_center, y_center, z_center = coords

    if x_range != y_range:
        raise RuntimeError('x and y resolutions must match for now.')

    # The galvo's small angle step response is 400 us
    # Let's give ourselves a buffer of 500 us (500000 ns)
    delay = int(0.5 * 10**6)

    total_num_samples = num_steps**2

    # %% Load the PulseStreamer

    # We require bookends on samples so stream one extra cycle
    seq_cycles = total_num_samples + 1
    period = cxn.pulse_streamer.stream_load('simple_readout.py', seq_cycles,
                                            [delay, readout, apd_index])

    # %% Set up the galvo

    x_voltages, y_voltages = cxn.galvo.load_sweep_scan(x_center, y_center,
                                                       x_range, y_range,
                                                       num_steps, period)

    x_num_steps = len(x_voltages)
    x_low = x_voltages[0]
    x_high = x_voltages[x_num_steps-1]
    y_num_steps = len(y_voltages)
    y_low = y_voltages[0]
    y_high = y_voltages[y_num_steps-1]

    pixel_size = x_voltages[1] - x_voltages[0]

    # %% Set the piezo

    cxn.objective_piezo.write_voltage(z_center)

    # %% Set up the APD

    cxn.apd_counter.load_stream_reader(apd_index, period, total_num_samples)

    # %% Set up the image display

    # Initialize imgArray and set all values to NaN so that unset values
    # are not interpreted as 0 by matplotlib's colobar
    img_array = numpy.empty((x_num_steps, y_num_steps))
    img_array[:] = numpy.nan
    img_write_pos = []

    # For the image extent, we need to bump out the min/max x/y by half the
    # pixel size in each direction so that the center of each pixel properly
    # lies at its x/y voltages.
    half_pixel_size = pixel_size / 2
    img_extent = [x_high + half_pixel_size, x_low - half_pixel_size,
                  y_low - half_pixel_size, y_high + half_pixel_size]

    fig = tool_belt.create_image_figure(img_array, img_extent,
                                        clickHandler=on_click_image)

    # %% Collect the data

    cxn.pulse_streamer.stream_start()

    timeout_duration = ((period*(10**-9)) * total_num_samples) + 10
    timeout_inst = time.time() + timeout_duration

    num_read_so_far = 0

    tool_belt.init_safe_stop()

    while num_read_so_far < total_num_samples:

        if time.time() > timeout_inst:
            log.failure('Timed out before all samples were collected.')
            break

        if tool_belt.safe_stop():
            break

        # Read the samples and update the image
        new_samples = cxn.apd_counter.read_stream(apd_index)
        num_new_samples = len(new_samples)
        if num_new_samples > 0:
            populate_img_array(new_samples, img_array, img_write_pos)
            tool_belt.update_image_figure(fig, img_array)
            num_read_so_far += num_new_samples

    # %% Save the data

    timeStamp = tool_belt.get_time_stamp()

    rawData = {'timeStamp': timeStamp,
               'name': name,
               'xyz_centers': [x_center, y_center, z_center],
               'x_range': x_range,
               'y_range': y_range,
               'num_steps': num_steps,
               'readout': int(readout),
               'resolution': [x_num_steps, x_num_steps],
               'img_array': img_array.astype(int).tolist()}

    filePath = tool_belt.get_file_path('scan_sample', timeStamp, name)
    tool_belt.save_figure(fig, filePath)
    tool_belt.save_raw_data(rawData, filePath)

    # %% Clean up

    # Stop the pulser
    cxn.pulse_streamer.constant_default()

    # Close tasks
    cxn.galvo.close_task()
    cxn.apd_counter.close_task(apd_index)

    # Return to center
    cxn.galvo.write(x_center, y_center)
