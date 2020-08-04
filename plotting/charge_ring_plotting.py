# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 09:05:07 2020

@author: agardill
"""

import matplotlib.pyplot as plt 
import numpy
import utils.tool_belt as tool_belt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy import integrate

# %%

def gaussian(r, constrast, sigma, center, offset):
    return offset+ constrast * numpy.exp(-((r-center)**2) / (2 * (sigma**2)))

def double_gaussian_dip(freq, low_constrast, low_sigma, low_center, low_offset,
                        high_constrast, high_sigma, high_center, high_offset):
    low_gauss = gaussian(freq, low_constrast, low_sigma, low_center, low_offset)
    high_gauss = gaussian(freq, high_constrast, high_sigma, high_center, high_offset)
    return low_gauss + high_gauss

def exponential(x, a, d, c):
    return c+a*(1-numpy.exp(-x/d))

def exponential_power(x, a, d, c, f):
    return c+a*(1-numpy.exp(-(x/d)**f))

def power_law(x, c, a, p):
    return c - a/x**p

def sqrt(x, a):
    return a*x**0.5


# %%
    
def r_vs_power_plot(nv_sig, ring_radius_list, ring_err_list, power_list, 
                    power_err_list, sub_folder, 
                    img_range, num_steps, green_pulse_time, readout):
    power_fig, ax = plt.subplots(1,1, figsize = (8, 8))
    ax.errorbar(power_list, ring_radius_list, xerr = power_err_list, yerr = ring_err_list, fmt = 'o')
    ax.set_xlabel('Green optical power (mW)')
    ax.set_ylabel('Charge ring radius (um)')
    ax.legend()
            
    timestamp = tool_belt.get_time_stamp()
    
    # save this dile 
    rawData = {'timestamp': timestamp,
               'nv_sig': nv_sig,
               'nv_sig-units': tool_belt.get_nv_sig_units(),
               'image_range': img_range,
               'image_range-units': 'V',
               'num_steps': num_steps,
               'green_pulse_time': green_pulse_time,
               'green_pulse_time-units': 'ns',
               'readout': readout,
               'readout-units': 'ns',
               'power_list': power_list,
               'power_list-units': 'mW',
               'power_err_list': power_err_list.tolist(),
               'power_err_list-units': 'mW',
               'ring_radius_list':ring_radius_list,
               'ring_radius_list-units': 'um',
               'ring_err_list': ring_err_list,
               'ring_err_list-units': 'um'}
    
    filePath = tool_belt.get_file_path("image_sample", timestamp, nv_sig['name'], subfolder = sub_folder)
    tool_belt.save_raw_data(rawData, filePath + '_radius_vs_power')
    
    tool_belt.save_figure(power_fig, filePath + '_radius_vs_power')
    
def r_vs_time_plot(nv_sig, ring_radius_list, ring_err_list, green_time_list, 
                     sub_folder, 
                    img_range, num_steps, green_pulse_time, readout):
    power_fig, ax = plt.subplots(1,1, figsize = (8, 8))
    ax.errorbar(numpy.array(green_time_list)/10**9, ring_radius_list, yerr = ring_err_list, fmt = 'o')
    ax.set_xlabel('Green pulse time (s)')
    ax.set_ylabel('Charge ring radius (um)')
    ax.legend()
            
    timestamp = tool_belt.get_time_stamp()
    
    # save this dile 
    rawData = {'timestamp': timestamp,
               'nv_sig': nv_sig,
               'nv_sig-units': tool_belt.get_nv_sig_units(),
               'image_range': img_range,
               'image_range-units': 'V',
               'num_steps': num_steps,
               'green_pulse_time': green_pulse_time,
               'green_pulse_time-units': 'ns',
               'readout': readout,
               'readout-units': 'ns',
               'green_time_list': green_time_list,
               'green_time_list-units': 'ns',
               'ring_radius_list':ring_radius_list,
               'ring_radius_list-units': 'um',
               'ring_err_list': ring_err_list,
               'ring_err_list-units': 'um'}
    
    filePath = tool_belt.get_file_path("image_sample", timestamp, nv_sig['name'], subfolder = sub_folder)
    tool_belt.save_raw_data(rawData, filePath + '_radius_vs_time')
    
    tool_belt.save_figure(power_fig, filePath + '_radius_vs_time')
# %%

def radial_distrbution_power(folder_name, sub_folder):
    labls = ['Area A', 'Area B', 'Area C']
   # create a file list of the files to analyze
    file_list  = tool_belt.get_file_list(folder_name, '.txt')
    file_list = [ 'A.txt','B.txt', 'C.txt']
    
    # create lists to fill with data
    power_list = []
    radii_array = []
    counts_r_array = []
    fig, ax = plt.subplots(1,1, figsize = (8, 8))
    l = 0
    for file in file_list:
#        try:
            data = tool_belt.get_raw_data(folder_name, file[:-4])
            # Get info from file
            timestamp = data['timestamp']
            nv_sig = data['nv_sig']
            coords = nv_sig['coords']
            x_voltages = data['x_voltages']
            y_voltages = data['y_voltages']
            num_steps = data['num_steps']
            img_range= data['image_range']
            dif_img_array = numpy.array(data['dif_img_array'])
            green_pulse_time = data['green_pulse_time']
            readout = data['readout']
#            opt_volt = data['green_optical_voltage']
            opt_power = data['green_opt_power']
            
            # Initial calculations
            x_coord = coords[0]          
            y_coord = coords[1]
            half_x_range = img_range / 2
            x_high = x_coord + half_x_range

            pixel_size = x_voltages[1] - x_voltages[0]
            half_pixel_size = pixel_size / 2
            
            # List to hold the values of each pixel within the ring
            counts_r = []
            # New 2D array to put the radial values of each pixel
            r_array = numpy.empty((num_steps, num_steps))
            
            # Calculate the radial distance from each point to center
            for i in range(num_steps):
                x_pos = x_voltages[i] - x_coord
                for j in range(num_steps):
                    y_pos = y_voltages[j]  - y_coord
                    r = numpy.sqrt(x_pos**2 + y_pos**2)
                    r_array[i][j] = r
            
            # define bound on each ring radius, which will be one pixel in size
            low_r = 0
            high_r = pixel_size
            
            # step throguh the radial ranges for each ring, add pixel within ring to list
            while high_r <= (x_high + half_pixel_size):
                ring_counts = []
                for i in range(num_steps):
                    for j in range(num_steps): 
                        radius = r_array[i][j]
                        if radius >= low_r and radius < high_r:
                            ring_counts.append(dif_img_array[i][j])
                # average the counts of all counts in a ring
                counts_r.append(numpy.average(ring_counts))
                # advance the radial bounds
                low_r = high_r
                high_r = high_r + pixel_size
            
            # define the radial values as center values of pizels along x, convert to um
            radii = numpy.array(x_voltages[int(num_steps/2):])*35
            # plot
            ax.plot(radii, counts_r, label  = labls[l])#'{} mW green pulse'.format('%.2f'%opt_power))
            power_list.append(opt_power)
            radii_array.append(radii.tolist())
            counts_r_array.append(counts_r)
            l = l+1
            
            integrated_counts = integrate.simps(counts_r, x = radii)
            print(integrated_counts)
            # try to fit the radial distribution to a double gaussian(work in prog)    
#            try:
#                contrast_low = 500
#                sigma_low = 5
#                center_low = -5
#                offset_low = 100
#                contrast_high = 300
#                sigma_high = 5
#                center_high = 20
#                offset_high = 100            
#                guess_params = (contrast_low, sigma_low, center_low, offset_low,
#                                contrast_high, sigma_high, center_high, offset_high)
#                
#                popt, pcov = curve_fit(double_gaussian_dip, radii[1:], counts_r[1:], p0=guess_params)
#                radii_linspace = numpy.linspace(radii[0], radii[-1], 1000)
#                
#                ax.plot(radii_linspace, double_gaussian_dip(radii_linspace, *popt))
#                print('fit succeeded')
#                
#                power_list.append(opt_power)
#                ring_radius_list.append(popt[6])
#                ring_err_list.append(pcov[6][6])
#                
#            except Exception:
#                print('fit failed' )
                
            
#        except Exception:
#            continue
        
    ax.set_xlabel('Radius (um)')
    ax.set_ylabel('Avg counts around ring (kcps)')
    ax.set_title('Varying position in diamond, similar depth\n50 s, 2 mW green pulse')
            #'Varying green pulse power, {} s'.format(green_pulse_time / 10**9))
    ax.legend()
 
    # save data from this file
    rawData = {'timestamp': timestamp,
               'file_list': file_list,
               'nv_sig': nv_sig,
               'nv_sig-units': tool_belt.get_nv_sig_units(),
               'num_steps': num_steps,
               'green_pulse_time': green_pulse_time,
               'green_pulse_time-units': 'ns',
               'readout': readout,
               'readout-units': 'ns',
               'power_list': power_list,
               'power_list-units':'mW',
               'radii_array': radii_array,
               'radii_array-units': 'um',
               'counts_r_array': counts_r_array,
               'counts_r_array-units': 'kcps'}
    
    filePath = tool_belt.get_file_path("image_sample", timestamp, nv_sig['name'], subfolder = sub_folder)
#    print(filePath)
    tool_belt.save_raw_data(rawData, filePath + '_radius')
    
    tool_belt.save_figure(fig, filePath + '_radius')
    
    return
                
# %%

def radial_distrbution_time(folder_name, sub_folder):
    # create a file list of the files to analyze
    file_list  = tool_belt.get_file_list(folder_name, '.txt')
    file_list = ['0.1.txt', '1.txt', '5.txt', '10.txt', '25.txt', '50.txt', '75.txt', '100.txt', '250.txt', '1000.txt' ]
    # create lists to fill with data
    green_time_list = []
    radii_array = []
    counts_r_array = []
    
    fig, ax = plt.subplots(1,1, figsize = (8, 8))
    for file in file_list:
        try:
            data = tool_belt.get_raw_data(folder_name, file[:-4])
            # Get info from file
            timestamp = data['timestamp']
            nv_sig = data['nv_sig']
            coords = nv_sig['coords']
            x_voltages = data['x_voltages']
            y_voltages = data['y_voltages']
            num_steps = data['num_steps']
            img_range= data['image_range']
            dif_img_array = numpy.array(data['dif_img_array'])
            green_pulse_time = data['green_pulse_time']
            readout = data['readout']
            opt_volt = data['green_optical_voltage']
            opt_power = data['green_opt_power']
            
            # Initial calculations
            x_coord = coords[0]          
            y_coord = coords[1]
            half_x_range = img_range / 2
            x_high = x_coord + half_x_range

            pixel_size = x_voltages[1] - x_voltages[0]
            half_pixel_size = pixel_size / 2
            
            # List to hold the values of each pixel within the ring
            counts_r = []
            # New 2D array to put the radial values of each pixel
            r_array = numpy.empty((num_steps, num_steps))
            
            # Calculate the radial distance from each point to center
            for i in range(num_steps):
                x_pos = x_voltages[i] - x_coord
                for j in range(num_steps):
                    y_pos = y_voltages[j]  - y_coord
                    r = numpy.sqrt(x_pos**2 + y_pos**2)
                    r_array[i][j] = r
            
            # define bound on each ring radius, which will be one pixel in size
            low_r = 0
            high_r = pixel_size
            
            # step throguh the radial ranges for each ring, add pixel within ring to list
            while high_r <= (x_high + half_pixel_size):
                ring_counts = []
                for i in range(num_steps):
                    for j in range(num_steps): 
                        radius = r_array[i][j]
                        if radius >= low_r and radius < high_r:
                            ring_counts.append(dif_img_array[i][j])
                # average the counts of all counts in a ring
                counts_r.append(numpy.average(ring_counts))
                # advance the radial bounds
                low_r = high_r
                high_r = high_r + pixel_size
            
            # define the radial values as center values of pizels along x, convert to um
            radii = numpy.array(x_voltages[int(num_steps/2):])*35
            # plot
#            fig, ax = plt.subplots(1,1, figsize = (8, 8))
            ax.plot(radii, counts_r, label  = '{} s green pulse'.format(green_pulse_time/10**9))
            green_time_list.append(green_pulse_time)
            radii_array.append(radii.tolist())
            counts_r_array.append(counts_r)
            
            integrated_counts = integrate.simps(counts_r, x = radii)
            print(integrated_counts)
            
#            # try to fit the radial distribution to a double gaussian(work in prog)    
#            try:
#                contrast_low = 500
#                sigma_low = 5
#                center_low = -5
#                offset_low = 100
#                contrast_high = 300
#                sigma_high = 5
#                center_high = 20
#                offset_high = 100            
#                guess_params = (contrast_low, sigma_low, center_low, offset_low,
#                                contrast_high, sigma_high, center_high, offset_high)
#                
#                popt, pcov = curve_fit(double_gaussian_dip, radii[1:], counts_r[1:], p0=guess_params)
#                radii_linspace = numpy.linspace(radii[0], radii[-1], 1000)
#                
#                ax.plot(radii_linspace, double_gaussian_dip(radii_linspace, *popt))
#                print('fit succeeded')
#                
#                green_time_list.append(green_pulse_time)
#                ring_radius_list.append(popt[6])
#                ring_err_list.append(pcov[6][6])
#                
#            except Exception:
#                print('fit failed' )
                
        except Exception:
            continue
    
    ax.set_xlabel('Radius (um)')
    ax.set_ylabel('Avg counts around ring (kcps)')
    ax.set_title('Varying green pulse time, {} mW'.format('%.2f'%opt_power))
    ax.legend()
 
    # save data from this file
    rawData = {'timestamp': timestamp,
               'file_list': file_list,
               'nv_sig': nv_sig,
               'nv_sig-units': tool_belt.get_nv_sig_units(),
               'num_steps': num_steps,
               'green_optical_voltage': opt_volt,
               'green_optical_voltage-units': 'V',
               'green_opt_power': opt_power,
               'green_opt_power-units': 'mW',
               'readout': readout,
               'readout-units': 'ns',
               'green_time_list': green_time_list,
               'green_time_list-units':'ns',
               'radii_array': radii_array,
               'radii_array-units': 'um',
               'counts_r_array': counts_r_array,
               'counts_r_array-units': 'kcps'}
    
    filePath = tool_belt.get_file_path("image_sample", timestamp, nv_sig['name'], subfolder = sub_folder)
    tool_belt.save_raw_data(rawData, filePath + '_radius')
    
    tool_belt.save_figure(fig, filePath + '_radius')
            
    return
    
# %%

def radial_distrbution_wait_time(folder_name, sub_folder):
    # create a file list of the files to analyze
    file_list  = tool_belt.get_file_list(folder_name, '.txt')
    file_list = [ '0.txt', '1.txt','100.txt', '1000.txt', '10000.txt', '100000.txt']
    # create lists to fill with data
    wait_time_list = []
    radial_counts_list = []
    
    fig, ax = plt.subplots(1,1, figsize = (8, 8))
              
    for file in file_list:
        try:
            data = tool_belt.get_raw_data(folder_name, file[:-4])
            # Get info from file
            timestamp = data['timestamp']
            nv_sig = data['nv_sig']
            coords = nv_sig['coords']
            x_voltages = data['x_voltages']
            y_voltages = data['y_voltages']
            num_steps = data['num_steps']
            img_range= data['image_range']
            dif_img_array = numpy.array(data['dif_img_array'])
            green_pulse_time = data['green_pulse_time']
            wait_time = data['wait_time']
            readout = data['readout']
            opt_volt = data['green_optical_voltage']
            opt_power = data['green_opt_power']
            
            # Initial calculations
            x_coord = coords[0]          
            y_coord = coords[1]
            half_x_range = img_range / 2
            x_high = x_coord + half_x_range

            pixel_size = x_voltages[1] - x_voltages[0]
            half_pixel_size = pixel_size / 2
            
            # List to hold the values of each pixel within the ring
            counts_r = []
            # New 2D array to put the radial values of each pixel
            r_array = numpy.empty((num_steps, num_steps))
            
            # Calculate the radial distance from each point to center
            for i in range(num_steps):
                x_pos = x_voltages[i] - x_coord
                for j in range(num_steps):
                    y_pos = y_voltages[j]  - y_coord
                    r = numpy.sqrt(x_pos**2 + y_pos**2)
                    r_array[i][j] = r
            
            # define bound on each ring radius, which will be one pixel in size
            low_r = 0
            high_r = pixel_size
            
            # step throguh the radial ranges for each ring, add pixel within ring to list
            while high_r <= (x_high + half_pixel_size):
                ring_counts = []
                for i in range(num_steps):
                    for j in range(num_steps): 
                        radius = r_array[i][j]
                        if radius >= low_r and radius < high_r:
                            ring_counts.append(dif_img_array[i][j])
                # average the counts of all counts in a ring
                counts_r.append(numpy.average(ring_counts))
                # advance the radial bounds
                low_r = high_r
                high_r = high_r + pixel_size
            
            # define the radial values as center values of pizels along x, convert to um
            radii = numpy.array(x_voltages[int(num_steps/2):])*35
            radial_counts_list.append(counts_r)
            wait_time_list.append(wait_time)
            # add to the plot
            ax.plot(radii, counts_r, label = '{} s wait time'.format(wait_time))
      
        except Exception:
            continue

            
    ax.set_xlabel('Radius (um)')
    ax.set_ylabel('Avg counts around ring (kcps)')
    ax.set_title('{} s, {} mW green pulse'.format(green_pulse_time/10**9, '%.2f'%opt_power))
    ax.legend()
    
    # save data from this file
    rawData = {'timestamp': timestamp,
               'nv_sig': nv_sig,
               'nv_sig-units': tool_belt.get_nv_sig_units(),
               'num_steps': num_steps,
               'green_pulse_time': green_pulse_time,
               'green_pulse_time-units': 'ns',
               'green_optical_voltage': opt_volt,
               'green_optical_voltage-units': 'V',
               'green_opt_power': opt_power,
               'green_opt_power-units': 'mW',
               'readout': readout,
               'readout-units': 'ns',
               'wait_time_list': wait_time_list,
               'wait_time_list-units': 's',
               'radii': radii.tolist(),
               'radii-units': 'um',
               'radial_counts_list': radial_counts_list,
               'radial_counts_list-units': 'kcps'}
    
    filePath = tool_belt.get_file_path("image_sample", timestamp, nv_sig['name'], subfolder = sub_folder)
    tool_belt.save_raw_data(rawData, filePath + '_radial_dist_vs_wait_time')
    
    tool_belt.save_figure(fig, filePath + '_radial_dist_vs_wait_time')

# %% 
if __name__ == '__main__':
    parent_folder = "image_sample/branch_Spin_to_charge/2020_07/"
    
#    sub_folder = "hopper_50s_power"
#    sub_folder = "hopper_10s_power"
#    sub_folder = "hopper_1s_power"
#    sub_folder = "spatial variation"
#    folder_name = parent_folder + sub_folder
#    
#    radial_distrbution_power(folder_name, sub_folder)
    
#    sub_folder = "hopper_0.3mw_time"
#    sub_folder = "hopper_2mw_time"
    sub_folder = "hopper_15mw_time"
    folder_name = parent_folder + sub_folder 
    
    radial_distrbution_time(folder_name, sub_folder)
#
    
#    radial_distrbution_wait_time(folder_name, sub_folder)



    # %% Manual data fitting for power

    # Determined by eye from radial plots    
    powers = [0.02, 0.04, 0.2, 0.25, 0.44, 0.55, 0.63, 0.75, 0.89, 0.98, 1.1, 1.3, 3.2, 15.1]   # mW
    power_err = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.07, 0.1]
    
    radius_1 = [0.1, 1.5, 4, 7, 11, 13.5, 14.5, 14.0, 15.2, 15.0, 16.0,16.3, 19.4, 22.1 ] # um
    radius_1_err = [0.1, 1, 3, 5, 2, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ]
    radius_10 = [0.5, 2.5, 8, 9.6, 12.8, 14.7, 15.3, 16.7, 17.3, 17.7, 17.8, 18.7, 21.8, 24.2] # um
    radius_10_err = numpy.ones(len(radius_10))
    radius_10_err = radius_10_err[:]*0.5 # um
    radius_50 = [1.5, 3.5, 10.5, 13.6, 15.4, 16.5, 17.4, 19.3, 18.4, 19.4, 20.2,20.5, 23.5,25  ] # um
    radius_50_err =  radius_10_err
    
    popt_1s, pcov = curve_fit(power_law, powers, radius_1, p0= (20, 1, 1))
    popt_10s, pcov = curve_fit(power_law, powers, radius_10, p0= (20, 1, 1))
    popt_50s, pcov = curve_fit(power_law, powers, radius_50, p0= (20, 1, 1))
    
    power_linspace = numpy.linspace(powers[0], powers[-1], 1000) 

#    fig, ax = plt.subplots(1,1, figsize = (8, 8))
#    ax.errorbar(powers, radius_1, xerr = power_err, yerr = radius_1_err, fmt = 'bo', label = '1 s green pulse')
#    ax.plot(power_linspace, power_law(power_linspace, *popt_1s), 'b-', label = '1 s fit')
#    ax.errorbar(powers, radius_10, xerr = power_err, yerr = radius_10_err, fmt = 'ro',  label = '10 s green pulse')
#    ax.plot(power_linspace, power_law(power_linspace, *popt_10s), 'r-', label = '10 s fit')   
#    ax.errorbar(powers, radius_50, xerr = power_err, yerr = radius_50_err, fmt = 'go',  label = '50 s green pulse')
#    ax.plot(power_linspace, power_law(power_linspace, *popt_50s), 'g-', label = '50 s fit')
#
##    ax.set_xscale('log')
#    ax.set_xlabel('Green optical power (mW)')
#    ax.set_ylabel('Charge ring radius (um)')
#    ax.legend()
#    ax.set_title('Charge ring radius vs green pulse power')
#  
#    text_eq = r'$C - A_0/P^a})$'
#    
#    text_1s = '\n'.join(('1 s fit',
#                     r'$C = {}$ um'.format('%.1f'%popt_1s[0]),
#                      r'$A_0 = {}$ um'.format('%.1f'%popt_1s[1]),
#                      r'$a = {}$'.format('%.1f'%popt_1s[2])))
#    text_10s= '\n'.join(('10 s fit',
#                     r'$C = {}$ um'.format('%.1f'%popt_10s[0]),
#                      r'$A_0 = {}$ um'.format('%.1f'%popt_10s[1]),
#                      r'$a = {}$'.format('%.1f'%popt_10s[2])))
#    text_50s = '\n'.join(('50 s fit',
#                          r'$C = {}$ um'.format('%.1f'%popt_50s[0]),
#                      r'$A_0 = {}$ um'.format('%.1f'%popt_50s[1]),
#                      r'$a = {}$'.format('%.1f'%popt_50s[2])))
#    
#    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#    ax.text(0.3, 0.6, text_eq, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#    ax.text(0.3, 0.52, text_1s, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#    ax.text(0.3, 0.35, text_10s, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#    ax.text(0.3, 0.18, text_50s, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
    
    # %% Manual data fitting for green time
    
    times = [ 
            #0.25, 0.75,
            1, 2.5, 5, 7.5 ,
            10, 25, 50, 75,
            100, 250, 500, 
            750,
            1000
            ]
    times_12 = [ 
            #0.25, 0.75,
            1, 2.5, 5, 7.5 ,
            10, 25, 50, 75,
            100, 250, 500, 
            750,
            1000
            ]
  
    radius_4mW = [ 5, 6.8, 8.5, 10, 11.2,  11.9, 14.2, 15, 15.5, 16.5, 17.1, 17.5, 18]
    radius_4mW_err = [ 1.0, 1, 1, 1, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
    radius_8mW = [ 14.5, 16.2, 16.6, 17.1, 17.3, 18, 19, 19.4, 20, 20.5, 20.9, 21.3, 21.2 ]
    radius_8mW_err = [ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
    radius_12mW = [#11, 12, 
                   16, 17, 18.2, 18.5, 18.5, 19.7,  20.9, 21, 21.5,  22.2, 22.5, 22.9, 22.5 ]
    radius_12mW_err = [#2, 2, 
                       1.0, 0.5, 0.5, 1, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,0.3]
    
    popt_4mw, pcov = curve_fit(power_law, times, radius_4mW, p0= (20, 30, 0.1))
    popt_8mw, pcov = curve_fit(power_law, times, radius_8mW, p0= (20, 30, 0.1))
    popt_12mw, pcov = curve_fit(power_law, times_12, radius_12mW, p0= (20, 20, 0.1))
    time_linspace = numpy.linspace(0.1, times_12[-1], 1000)         
    
#    fig, ax = plt.subplots(1,1, figsize = (8, 8))
#    ax.errorbar(times, radius_4mW, yerr = radius_4mW_err, fmt = 'bo', label = '0.3 mW green pulse')
#    ax.plot(time_linspace, power_law(time_linspace, *popt_4mw), 'b-', label = '0.3 mW fit')
#    ax.errorbar(times, radius_8mW, yerr = radius_8mW_err, fmt = 'ro', label = '0.75 mW green pulse')
#    ax.plot(time_linspace, power_law(time_linspace, *popt_8mw), 'r-', label = '0.75 mW fit')
#    ax.errorbar(times_12, radius_12mW, yerr = radius_12mW_err, fmt = 'go', label = '1.3 mW green pulse')
#    ax.plot(time_linspace, power_law(time_linspace, *popt_12mw), 'g-', label = '1.3 mW fit')
##    ax.set_xscale('log')
#    ax.set_xlabel('Green pulse time (s)')
#    ax.set_ylabel('Charge ring radius (um)')
#    ax.legend()
#    ax.set_title('Charge ring radius vs green pulse length')
#    text_eq = r'$C - A_0/P^a})$'
#    
#    text_4mw = '\n'.join(('0.3 mW fit',
#                      r'$C = {}$ um'.format('%.1f'%popt_4mw[0]),
#                      r'$A_0 = {}$ um'.format('%.1f'%popt_4mw[1]),
#                      r'$a = {}$'.format('%.1f'%popt_4mw[2])))
#    text_8mw = '\n'.join(('0.75 mW fit',
#                      r'$C = {}$ um'.format('%.1f'%popt_8mw[0]),
#                      r'$A_0 = {}$ um'.format('%.1f'%popt_8mw[1]),
#                      r'$a = {}$'.format('%.1f'%popt_8mw[2])))
#    text_12mw = '\n'.join(('1.3 mW fit',
#                      r'$C = {}$ um'.format('%.1f'%popt_12mw[0]),
#                      r'$A_0 = {}$ um'.format('%.1f'%popt_12mw[1]),
#                      r'$a = {}$'.format('%.1f'%popt_12mw[2])))
#    
#    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#    ax.text(0.3, 0.6, text_eq, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#    ax.text(0.3, 0.52, text_4mw, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#    ax.text(0.3, 0.35, text_8mw, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)
#    ax.text(0.3, 0.18, text_12mw, transform=ax.transAxes, fontsize=12,
#            verticalalignment='top', bbox=props)