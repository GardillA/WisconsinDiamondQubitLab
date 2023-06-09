# -*- coding: utf-8 -*-
"""
Template matching.

Created on Fri Jan  20 08:30:00 2023

@author: carter fox 
"""

####################### Imports #######################

import cv2
import numpy as np
from matplotlib import pyplot as plt
import json
import utils.tool_belt as tool_belt
import sys
import utils.positioning as positioning
import utils.tool_belt as tool_belt
import utils.common as common

def convert_to_8bit(img, min_val=0, max_val=255):
    img = img.astype(np.float64)
    img -= np.nanmin(img)  # Set the lowest value to 0
    img *= (max_val/np.nanmax(img))
    img = img.astype(np.uint8)
    img += min_val
    return img

def get_img_extent(data):
    x_voltages = data['x_positions_1d']
    y_voltages = data['y_positions_1d']
    x_low = x_voltages[0]
    x_high = x_voltages[-1]
    y_low = y_voltages[0]
    y_high = y_voltages[-1]

    half_pixel_size = (x_voltages[1] - x_voltages[0]) / 2
    img_extent = [x_low - half_pixel_size, x_high + half_pixel_size,
                  y_low - half_pixel_size, y_high + half_pixel_size]
    return img_extent


def get_shift(nv_sig, haystack_file, needle_file, close_plot=False):
    '''
    This routine will take a "needle" file and use it as a template to 
    find it's relative position in a "haystack" file. The two image files have 
    to have the same resolution (V / pixels) to work.
    
    It will return the shift in X and Y of the needle file with respect to the
    central coordinates of the haystack file.
    
    This is intended to be used for a large haystack file that is centered on 
    an NV of interest. The autotracker can then be run with a needle scan
    to find the offset from the original position. The opposite sign of the 
    x and y shifts should be added to the drift.
    
    Ex: this fundtion returns an x_shift = 0.04 and y_shift = -0.05, then the 
    values of -0.04 and 0.05 should be added to the x and y drift, respectively.
    '''
    
    diff_lim_spot_diam = 0.015  # expected st dev of the gaussian in volts
    
    ####################### Process Haystack File #######################
        
    haystack_data = tool_belt.get_raw_data(haystack_file)
    haystack_x_range = haystack_data['x_range']
    haystack_num_steps = haystack_data['num_steps']
    
    x_voltages = haystack_data['x_positions_1d']
    y_voltages = haystack_data['y_positions_1d']
    
    haystack_img_extent = get_img_extent(haystack_data)
    og_center_x_volts = (x_voltages[-1] + x_voltages[0])/2
    og_center_y_volts = (y_voltages[-1] + y_voltages[0])/2
    
    x_range = haystack_data['x_range']
    y_range = haystack_data['y_range']
    min_range = min(x_range, y_range)
    haystack_num_steps = haystack_data['num_steps']
    volts_per_pixel = min_range / haystack_num_steps
    
    haystack_img_array = np.array(haystack_data['img_array'])
    haystack_img_array = convert_to_8bit(haystack_img_array)
    
    # Blur the haystack image wit ha Gaussian Blur
    haystack_img_array = cv2.GaussianBlur(haystack_img_array, (5, 5), 0)
    
    ####################### Initial calculations #######################
    
    # expected st dev of the gaussian in pixels - must be odd
    diff_lim_spot_pixels = int(diff_lim_spot_diam / volts_per_pixel)
    if diff_lim_spot_pixels % 2 == 1:
        diff_lim_spot_pixels += 1
    
    ####################### Process Needle File #######################
        
    needle_data = tool_belt.get_raw_data(needle_file)
    needle_img_array = np.array(needle_data['img_array'])
    needle_x_range = needle_data['x_range']
    needle_num_steps = needle_data['num_steps']
    
    if needle_x_range/needle_num_steps != haystack_x_range/haystack_num_steps:
        print('images must have same scale of volts per pixel!')
        raise RuntimeError
    
    needle_x = needle_data['x_positions_1d']
    needle_y = needle_data['y_positions_1d']
    w, h = len(needle_x),len(needle_y)
    needle_img_extent = get_img_extent(needle_data)
    needle_img_array = convert_to_8bit(needle_img_array)
    
    
    ####################### Run the matching #######################
    
    method = eval('cv2.TM_CCOEFF_NORMED')
    res = cv2.matchTemplate(haystack_img_array, needle_img_array, method)
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    top_left = max_loc
    bottom_right = (top_left[0] + w, top_left[1] + h)
    
    center_x = int( w/2 + top_left[0] )
    center_y = int( h/2 + top_left[1] )
    
    # if x_voltages[1] > x_voltages[0]:
    #     x_voltages.reverse()
    if y_voltages[1] > y_voltages[0]:
        y_voltages.reverse()
    center_x_volts = x_voltages[center_x]# + volts_per_pixel/2
    center_y_volts = y_voltages[center_y]# + volts_per_pixel/2
    # print(center_x_volts,center_y_volts)
    # print(center_x_volts)
    # print(og_center_x_volts)

    # print(manual_coor_shift)
    shift_x_volts = round(center_x_volts - og_center_x_volts,4) 
    shift_y_volts = round(center_y_volts - og_center_y_volts,4)
    
    processed_img_array = np.copy(haystack_img_array)
    processed_img_array = cv2.cvtColor(processed_img_array, cv2.COLOR_GRAY2RGB)
    
    ####################### Plotting #######################
    
    rawData = {'haystack_fname':haystack_file,
               'needle_fname':needle_file,
               'shift_x_volts':shift_x_volts,
               'shift_y_volts':shift_y_volts,
               'haystack_data':haystack_data,
               'needle_data':needle_data,
               'nv_sig':needle_data['nv_sig']
               }
    timestamp = tool_belt.get_time_stamp()
    fname = 'auto_tracker'
    filePath = tool_belt.get_file_path(__file__, timestamp, fname)
    tool_belt.save_raw_data(rawData, filePath)
    
    fig, axes_pack = plt.subplots(1, 3, figsize=(15, 5))
    ax = axes_pack[0]
    ax.set_title('haystack processed')
    ax.set_xlabel('x [V]')
    ax.set_ylabel('y [V]')
    ax.imshow(haystack_img_array, extent=tuple(haystack_img_extent))
    ax = axes_pack[1]
    ax.set_xlabel('x [V]')
    ax.set_ylabel('y [V]')
    ax.imshow(needle_img_array, extent=tuple(needle_img_extent))
    ax.set_title('needle')
    ax = axes_pack[2]
    ax.set_xlabel('x [V]')
    ax.set_ylabel('y [V]')
    ax.imshow(processed_img_array, extent=tuple(haystack_img_extent))
    cv2.rectangle(haystack_img_array,top_left, bottom_right, 255, 1)
    cv2.circle(haystack_img_array, [center_x,center_y], 0, (255, 0, 0), 1)
    plt.imshow(haystack_img_array, extent=tuple(haystack_img_extent))
    ax.set_title('x_shift, y_shift = {} V, {} V'.format(shift_x_volts,shift_y_volts))
    # fig.tight_layout()
    # fig_manager = plt.get_current_fig_manager()
    tool_belt.save_figure(fig, filePath)
        
    if close_plot:
        plt.close()
        
    return shift_x_volts, shift_y_volts

    
if __name__ == '__main__':
    
    ####################### Files #######################
    needle_file = '2023_06_09-10_16_53-WiQD-nv1_XY'
    # needle_file = '2023_06_09-10_12_12-WiQD-nv1_XY'
    haystack_file = '2023_06_08-16_09_16-WiQD-nv1_XY'
    data_needle = tool_belt.get_raw_data(needle_file)
    nv_sig= data_needle['nv_sig']
    
    shift_x_volts, shift_y_volts = get_shift(nv_sig, haystack_file, needle_file,close_plot=False)
    print('x shift = {} V'.format(shift_x_volts))
    print('y shift = {} V'.format(shift_y_volts))