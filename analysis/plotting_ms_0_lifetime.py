# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 09:51:31 2019

Plot the (0,0) data

@author: Aedan
"""

# %% Imports

import numpy
from scipy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import utils.tool_belt as tool_belt

# %% Functions

# The exponential function used to fit the data

def exp_eq(t, offset, rate, amp):
    return offset+amp * exp(- rate * t)

# %%
    
# Area A5
data = tool_belt.get_raw_data('t1_double_quantum/data_folders/other_data/bachman-A5-ensemble-B1-234MHz/', 
                              '2020_06_01-07_06_15-bachman-A1')

relaxation_time_range = numpy.array(data['relaxation_time_range'])/10**6
num_steps = data['num_steps']
norm_avg_sig_A1 = data['norm_avg_sig']
#sig_counts_A1 = data["sig_counts"]
#ref_counts_A1 = data["ref_counts"]
taus_A1 = numpy.linspace(relaxation_time_range[0], relaxation_time_range[1], num_steps)

# manipulate the data to normalize 
first_point = norm_avg_sig_A1[0]
last_point = norm_avg_sig_A1[-1]

norm_avg_sig_A1 = (numpy.array(norm_avg_sig_A1) - last_point)/ (first_point - last_point)

offset = 0.9
amplitude = 0.1
decay = 0.6*3 # inverse ns

#popt_A1, pcov = curve_fit(exp_eq, taus_A1, norm_avg_sig_A1,
#                          p0=[offset, decay, amplitude])

linspace_tau_A1 = numpy.linspace(relaxation_time_range[0], relaxation_time_range[1], 1000)

    
# Area A1
data = tool_belt.get_raw_data('t1_double_quantum/data_folders/other_data/bachman-A1-ensemble-B1-232MHz/', 
                              '2020_06_16-11_58_29-Bachman-A1-B1')

relaxation_time_range = numpy.array(data['relaxation_time_range'])/10**6
num_steps = data['num_steps']
norm_avg_sig_B1 = data['norm_avg_sig']
taus_B1 = numpy.linspace(relaxation_time_range[0], relaxation_time_range[1], num_steps)

# manipulate the data to normalize 
first_point = norm_avg_sig_B1[0]
last_point = norm_avg_sig_B1[-3]

norm_avg_sig_B1 = (numpy.array(norm_avg_sig_B1) - last_point)/(first_point - last_point)

offset = 0.9
amplitude = 0.1
decay = 0.6*3 # inverse ns

#popt_B1, pcov = curve_fit(exp_eq, taus_B1, norm_avg_sig_B1,
#                          p0=[offset, decay, amplitude])

linspace_tau_B1 = numpy.linspace(relaxation_time_range[0], relaxation_time_range[1], 1000)
    
fig, ax = plt.subplots(1,1, figsize=(10, 8))
ax.semilogy(taus_A1, norm_avg_sig_A1,'ko',
                    label = 'Area A5')
#ax.plot(linspace_tau_A1,
#                    exp_eq(linspace_tau_A1, *popt_A1),
#                    'k-', label = 'B5 fit')
ax.semilogy(taus_B1, norm_avg_sig_B1, 'ro',
                    label = 'Area A1')
#ax.plot(linspace_tau_B1,
#                    exp_eq(linspace_tau_B1, *popt_B1),
#                    'r-', label = 'B1 fit')
ax.set_xlabel('Wait time (ms)')
ax.set_ylabel('Normalized signal Counts')
ax.set_title('Lifetime of ms = 0')
ax.legend()


