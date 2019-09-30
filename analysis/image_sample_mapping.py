# -*- coding: utf-8 -*-
"""Saves the data necessary to relocate specific NVs. Also can probably (?)
illustrate a mapping from an NV list to an image_sample.

Created on Mon Jun 10 13:54:07 2019

@author: mccambria
"""


# %% Imports


import json
import majorroutines.image_sample as image_sample
import utils.tool_belt as tool_belt
import matplotlib.pyplot as plt
import os
import labrad
from pathlib import Path


# %% Functions


def illustrate_mapping(file_name):

    data = tool_belt.get_raw_data(__file__, file_name)
    image_sample_file_name = data['image_sample_file_name']
    nv_sig_list = data['nv_sig_list']
        
    fig = image_sample.create_figure(image_sample_file_name)
    axes = fig.get_axes()
    ax = axes[0]
    images = ax.get_images()
    im = images[0]
    im.set_clim(0, 100)
    fig.set_size_inches(8.5, 8.5)

    # Get the expected radius of an NV
    try:
        with labrad.connect() as cxn:
            shared_params = tool_belt.get_shared_parameters_dict(cxn)
        airy_radius_nm = shared_params['airy_radius']
        galvo_nm_per_volt = shared_params['galvo_nm_per_volt']
        airy_radius_volts = airy_radius_nm / galvo_nm_per_volt
    except Exception:
        airy_radius_volts = 0.004

    for ind in range(len(nv_sig_list)):
        coords = nv_sig_list[ind]['coords']
        circle = plt.Circle(tuple(coords[0:2]), 2*airy_radius_volts,
                            ec='g', fill=False, lw=2.0)
        ax.add_patch(circle)
        
    return fig
    
def generate_mapping_files(sample_name, micrometer_coords,
                           image_sample_file_name, nv_sig_list):
    
    raw_data = {
            'sample_name': sample_name,
            'micrometer_coords': micrometer_coords,
            'micrometer_coords-units': 'mm',
            'image_sample_file_name': image_sample_file_name,
            'nv_sig_list': nv_sig_list,
            'nv_sig_list-units': tool_belt.get_nv_sig_units(),
            }
    
    file_name = '{}-mapping'.format(image_sample_file_name)
    file_path = tool_belt.get_file_path(__file__, name=file_name)

    tool_belt.save_raw_data(raw_data, file_path)
    fig = illustrate_mapping(file_name)
    
    tool_belt.save_figure(fig, file_path)


# %% Run the file


if __name__ == '__main__':
    
#    image_sample_file_name = '2019-07-25_18-37-46_ayrton12_search'

    # Ignore this...
#    if False:
    # Circle NVs from an existing mapping
    file_name = '2019-06-10_15-26-39_ayrton12_mapping'
    illustrate_mapping(file_name)
#    else:

#    coords_list = [   [-0.379, 0.280]]
#
#    sample_name = 'ayrton12'
#    micrometer_coords = [3.154]
#    image_sample_file_name = '2019-07-25_18-37-46_ayrton12_search'
#
#    nv_sig_list = []
#    for ind in range(len(coords_list)):
#        coords = coords_list[ind]
#        name = '{}-nv{}_2019_07_25'.format(sample_name, ind)
#        nd_filter = 'nd_1.5'
#        nv_sig = {'coords': coords, 'name': name, nd_filter: nd_filter}
#        nv_sig_list.append(nv_sig)
#    
#    generate_mapping_files(sample_name, micrometer_coords,
#                          image_sample_file_name, nv_sig_list)
