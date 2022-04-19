# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:56:11 2021

@author: agardill
"""

import matplotlib.pyplot as plt 
import numpy
import math
import utils.tool_belt as tool_belt
# import majorroutines.image_sample as image_sample
import csv

min_counts = 0
max_counts = 12

def create_image_figure(imgArray, imgExtent, title):
    f_size = 8
    tick_f_size = 8
    clb_f_size = 8
    fig_w = 4#2
    fig_l =4#2
    
    mu = u"\u03BC" 
    
    fig_tick_l = 3
    fig_tick_w = 0.75
    
    clb_tick_1 = 3
    clb_tick_w = 0.75

    
    # Tell matplotlib to generate a figure with just one plot in it
    fig, ax = plt.subplots(dpi=300)
    fig.set_figwidth(fig_w)
    fig.set_figheight(fig_l)
    
    # ax.set_xticks([-1000,0,1000])
    # ax.set_yticks([-1000,0,1000])
    ax.tick_params(which = 'both', length=fig_tick_l, width=fig_tick_w, colors='k',
                    direction='in',grid_alpha=0.7, labelsize = tick_f_size)
    

        
    # Tell the axes to show an image
    img = ax.imshow(imgArray, cmap='inferno', vmin = 0, vmax = 1,
                    extent=tuple(imgExtent))

    # Add a colorbar
    clb = plt.colorbar(img)
    clb.set_label("NV- population", rotation=270, fontsize = f_size)
    clb.ax.tick_params( length=clb_tick_1, width=clb_tick_w, grid_alpha=0.7, labelsize = clb_f_size)
    # clb.set_ticks([0, 0.2, 0.4, 0.6,0.8,1])
    
    
    # Label axes
    ax.set_xlabel('x (nm)', fontsize = f_size)
    ax.set_ylabel('y (nm)', fontsize = f_size)
    ax.set_title(title)


    # Draw the canvas and flush the events to the backend
    fig.canvas.draw()
    # plt.tight_layout()
    fig.canvas.flush_events()

    return fig


def create_figure(file_name,  sub_folder, threshold, scale):
    data = tool_belt.get_raw_data(file_name, sub_folder)
    nv_sig = data['nv_sig']
    CPG_laser_dur = nv_sig['CPG_laser_dur']
    # readout_counts_avg = numpy.array(data['readout_counts_avg'])
    # readout_counts_array = numpy.array(data['readout_counts_array'])
    num_steps_b = data['num_steps_b']    
    a_voltages_1d = data['a_voltages_1d']
    b_voltages_1d = data['b_voltages_1d']
    img_range_2D= data['img_range_2D']
    # offset_2D = data["offset_2D"]
    # drift_list = data['drift_list_master']
    axes = [0,1]
    
    # convert single shot measurements to NV- population
    raw_counts = numpy.array(data['readout_counts_array'])
    for r in range(len(raw_counts)):
        row = raw_counts[r]
        for c in range(len(row)):
            current_val = raw_counts[r][c]
            if current_val < threshold:
                set_val = 0
            elif current_val >= threshold:
                set_val = 1
            raw_counts[r][c] = set_val
    readout_counts_avg = numpy.average(raw_counts, axis = 1)
    
    half_range_a = img_range_2D[axes[0]]/2
    half_range_b = img_range_2D[axes[1]]/2
    a_low = -half_range_a 
    a_high = half_range_a
    b_low = -half_range_b
    b_high = half_range_b

    pixel_size_a = (a_voltages_1d[1] - a_voltages_1d[0])
    pixel_size_b = (b_voltages_1d[1] - b_voltages_1d[0])

    half_pixel_size_a = pixel_size_a / 2
    half_pixel_size_b = pixel_size_b / 2
    
    img_extent = [(a_low - half_pixel_size_a)*scale,
                  (a_high + half_pixel_size_a)*scale, 
                 
                 (b_low - half_pixel_size_b)*scale, 
                 (b_high + half_pixel_size_b)*scale ]
    um_scaled = True

    split_counts = numpy.split(readout_counts_avg, num_steps_b)
    readout_image_array = numpy.vstack(split_counts)
    r = 0
    
    for i in reversed(range(len(readout_image_array))):
        if r % 2 == 0:
            readout_image_array[i] = list(reversed(readout_image_array[i]))
        r += 1
    
    readout_image_array = numpy.flipud(readout_image_array)
    title = 'SPaCE - {} ms depletion pulse'.format(CPG_laser_dur)
    


    title = '{} ms'.format(CPG_laser_dur/10**6)
    fig = create_image_figure(readout_image_array, img_extent, title)
    # Redraw the canvas and flush the changes to the backend
    fig.canvas.draw()
    fig.canvas.flush_events()
        
    # file_path = 'E:/Shared drives/Kolkowitz Lab Group/nvdata/paper_materials/charge siv/Figures/Fig 1 - E6 and SiV comparison/'
    # file_name = timestamp + '-' + nv_sig['name']
#    tool_belt.save_figure(fig, file_path + file_name +'-charge')
    

    return fig  

    
# %%
if __name__ == '__main__':
    
    folder = 'pc_rabi/branch_CFMIII/SPaCE_digital/2021_12'
    
    # 5 ms zoom
    file_list = [
        '2021_12_27-11_09_25-johnson-nv0_2021_12_22',
        '2021_12_27-11_17_00-johnson-nv0_2021_12_22',
        '2021_12_27-11_18_06-johnson-nv0_2021_12_22',
        '2021_12_27-11_19_28-johnson-nv0_2021_12_22',
        '2021_12_27-11_20_52-johnson-nv0_2021_12_22',
        '2021_12_27-11_21_52-johnson-nv0_2021_12_22',
        '2021_12_27-11_24_00-johnson-nv0_2021_12_22',
        '2021_12_27-11_24_53-johnson-nv0_2021_12_22',
        '2021_12_27-11_25_55-johnson-nv0_2021_12_22',
        '2021_12_27-11_27_01-johnson-nv0_2021_12_22',
        '2021_12_27-11_28_04-johnson-nv0_2021_12_22'
    ]
    threshold = 9
    for file in file_list:
        create_figure(file,  folder, threshold, 0.99e3)
        
        
    folder = 'pc_rabi/branch_CFMIII/SPaCE_digital/2021_11' #don't forget to get rid of vmax/vmin in plot function
    file = '2021_11_17-15_38_14-johnson-nv1_2021_11_17'
    
    # create_figure(file,  folder, 1000)