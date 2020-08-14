# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:17:30 2020

Plot spectra data from a json file with the wavelengths and counts

Specifically for plotting data from 3/10 and 3/9 of Er samples without graphene

@author: agardill
"""

import utils.tool_belt as tool_belt
from scipy.optimize import curve_fit
import numpy
import matplotlib.pyplot as plt

# %%

props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

#data_path = 'C:/Users/Public/Documents/Jobin Yvon/SpectraData'

data_path = 'E:/Shared drives/Kolkowitz Lab Group/nvdata/spectra/Brar'

#folder = '5_nm_Er_graphene/2020_03_02'

#folder = '5_nm_Er_NO_graphene_NO_ig'

# %%
def wavelength_range_calc(wavelength_range, wavelength_list):
    wvlngth_range_strt = wavelength_range[0]
    wvlngth_range_end = wavelength_range[1]
    
    wvlngth_list_strt = wavelength_list[0]
    wvlngth_list_end = wavelength_list[-1]
    step_size = wavelength_list[1] - wavelength_list[0]
    
    # Try to calculate the indices for plotting, if they are given
    try: 
        start_index = int( (wvlngth_range_strt - wvlngth_list_strt) /  step_size )
    except Exception:
        pass
    try: 
        end_index = -(int((wvlngth_list_end - wvlngth_range_end) /  step_size) + 1)
    except Exception:
        pass
            
    
    #If both ends of range are specified
    if wvlngth_range_strt and wvlngth_range_end:
        
        #However, if the range is outside the actual data, just plot all data
        # Or if one of the bounds is outside the range, do not include that range
        if wvlngth_range_strt < wvlngth_list_strt and \
                                    wvlngth_range_end > wvlngth_list_end:
            print('Wavelength range outside of data range. Plotting full range of data')
            plot_strt = 0
            plot_end = -1
        elif wvlngth_range_strt < wvlngth_list_strt and \
            wvlngth_range_end <= wvlngth_list_end:
                
            print('Requested lower bound of wavelength range outside of data range. Starting data at lowest wavelength')
            plot_strt = 0
            plot_end = end_index
            
        elif wvlngth_range_strt >= wvlngth_list_strt and \
            wvlngth_range_end > wvlngth_list_end:
                
            print('Requested upper bound of wavelength range outside of data range. Plotting up to highest wavelength')
            plot_strt = start_index
            plot_end = -1
        else:
            #If they are given, then plot the requested range
            plot_strt = start_index
            plot_end = end_index
            
    # If the lower bound is given, but not the upper bound
    elif wvlngth_range_strt and not wvlngth_range_end:
        plot_end = -1
        
        if wvlngth_range_strt < wvlngth_list_strt:
            print('Requested lower bound of wavelength range outside of data range. Plotting full range of data')
            plot_strt = 0
        else:
            plot_strt = start_index
            
    # If the upper bound is given, but not the lower bound
    elif wvlngth_range_end and not wvlngth_range_strt:
        plot_strt = 0
        
        if wvlngth_range_end > wvlngth_list_end:
            print('Requested upper bound of wavelength range outside of data range. Plotting full range of data')
            plot_end = -1
        else:
            plot_end = end_index
    # Lastly, if no range is given
    elif not wvlngth_range_end and not wvlngth_range_strt:
       plot_strt = 0 
       plot_end = -1
        
    return plot_strt, plot_end

# %%
    
def plot_spectra(file,folder, wavelength_range = [None, None], vertical_range = [None, None], plot_title = ''):
    data = tool_belt.get_raw_data(folder, file,
                 nvdata_dir=data_path)
    wavelengths = numpy.array(data['wavelengths'])
    counts = numpy.array(data['counts'])
    
    plot_strt_ind, plot_end_ind  = wavelength_range_calc(wavelength_range, wavelengths)
    
    # subtract off a constant background
    counts_cnst_bkgd = counts - numpy.average(counts[plot_end_ind-8:plot_end_ind])
    
    return wavelengths[plot_strt_ind : plot_end_ind], counts[plot_strt_ind : plot_end_ind]
        
        
        
def plot_spectra_list(parent_folder, file_list, title, label_list, y_range, x_range = None):
    fig, ax= plt.subplots(1, 1, figsize=(10, 8))
    
    for f in range(len(file_list)):
        file = file_list[f]
        wvlngth, counts = plot_spectra(file, parent_folder)
        ax.plot(wvlngth, counts, label =label_list[f])
        
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Counts')
    ax.set_title(title)
    ax.set_ylim(y_range)
    if x_range:
        ax.set_xlim(x_range)
    ax.legend()    
   
    
def august_cap_noncap_plots():
    label_list = ['8/10/2020', '8/11/2020', '8/12/2020', '8/13/2020', '8/14/2020']
    
    # capped
    title = 'Capped 5 nm Er'
    parent_folder = '2020_08_10 5 nm capped'
    file_list = ['2020_08_10-c-550','2020_08_11-c-550','2020_08_12-c-550',
                 '2020_08_13-c-550','2020_08_14-c-550'] 
    plot_spectra_list(parent_folder, file_list, title, label_list, y_range =[500,1500], x_range = [547, 580] )

    file_list = ['2020_08_10-c-670','2020_08_11-c-670','2020_08_12-c-670',
                 '2020_08_13-c-670','2020_08_14-c-670'] 
    plot_spectra_list(parent_folder, file_list, title, label_list, y_range =[580,750], x_range = [644, 692])
    
    # noncapped
    title = 'Noncapped 5 nm Er'
    parent_folder = '2020_08_10 5 nm noncapped'
    file_list = ['2020_08_10-nc-550','2020_08_11-nc-550','2020_08_12-nc-550',
                 '2020_08_13-nc-550','2020_08_14-nc-550'] 
    plot_spectra_list(parent_folder, file_list, title, label_list, y_range =[500,1500], x_range = [547, 580])

    file_list = ['2020_08_10-nc-670','2020_08_11-nc-670','2020_08_12-nc-670',
                 '2020_08_13-nc-670','2020_08_14-nc-670'] 
    plot_spectra_list(parent_folder, file_list, title, label_list, y_range =[580,750], x_range = [644, 692])
    
    
    
    
# %%
    
# capped
folder_c = '2020_08_10 5 nm capped'
#file_c = '2020_08_14-c-550'
file_c = '2020_08_14-c-670'

# noncapped
folder_nc = '2020_08_10 5 nm noncapped'
#file_nc = '2020_08_14-nc-550'
file_nc = '2020_08_14-nc-670'


if __name__ == '__main__':
#    august_cap_noncap_plots()
    
    wvlngth_1, counts_1 = plot_spectra(file_c, folder_c) 

    wvlngth_2, counts_2 = plot_spectra(file_nc, folder_nc)
 


    fig, ax= plt.subplots(1, 1, figsize=(10, 8))
    
#    print(counts_1)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Counts')
    ax.set_title('Capped vs noncapped (8/14/2020)')
    ax.set_ylim([580, 750]) 
    ax.set_xlim([644, 692]) 
    ax.plot(wvlngth_1, numpy.array(counts_1), label ='capped')
    ax.plot(wvlngth_2, numpy.array(counts_2), label = 'noncapped')
    ax.legend()

    
    
    
    