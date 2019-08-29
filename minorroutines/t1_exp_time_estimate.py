# -*- coding: utf-8 -*-
"""

Calculate the estimated time for the t1 experiment, how often we optimize, and
the expected error

Created on Thu Aug  1 17:32:46 2019

@author: agardill
"""
# %%
from utils.tool_belt import States
import numpy
# %%

def expected_st_dev_norm(ref_counts, expected_contrast):
    # Expected contrast is the difference between 1 and the minimum signal
    sig_counts = (1 - expected_contrast) * ref_counts
    rel_std_sig = numpy.sqrt(sig_counts) / sig_counts
    rel_std_ref = numpy.sqrt(ref_counts) / ref_counts
    # Propogate the error
    error = (1 - expected_contrast) * numpy.sqrt((rel_std_sig**2) + (rel_std_ref**2))

    return error

def t1_exp_times(exp_array, contrast, exp_count_rate, readout_window):
    total_exp_time_list = []

    for line in exp_array:
        init = line[0][0]
        read = line[0][1]
        relaxation_time_s = line[1][1] * 10**-9
        num_steps = line[2]
        num_reps = line[3]
        num_runs = line[4]
        extra_seq_time = 20 * 10**-6
        optimize_time = 5

        exp = [init.name, read.name]
        sequence_time = (relaxation_time_s + extra_seq_time) * num_reps
        exp_time_s = (sequence_time * num_steps / 2 + optimize_time) * num_runs # seconds
        exp_time_m = exp_time_s / 60
        exp_time_h = exp_time_m / 60

        opti_time = exp_time_m / num_runs

        total_exp_time_list.append(exp_time_h)

        ref_counts = (exp_count_rate * 10 ** 3) * (readout_window * 10**-9) * num_reps * num_runs

        exp_error = expected_st_dev_norm(ref_counts, contrast)

        snr = (contrast) / exp_error
        print('{}: {} hours, optimize every {} minutes, expected snr: {}'.format(exp, '%.1f'%exp_time_h, '%.2f'%opti_time, '%.1f'%snr))

    total_exp_time = sum(total_exp_time_list)


    print('Total experiment time: {} hrs'.format('%.1f'%total_exp_time))

# %%

num_runs = 75
t1_exp_array = numpy.array([
                    [[States.HIGH, States.LOW], [0, 50*10**3], 51, 3*10**4, num_runs],
                    [[States.HIGH, States.LOW], [0, 2.5*10**6], 26, 1*10**4, num_runs],
                    [[States.HIGH, States.HIGH], [0, 50*10**3], 51, 3*10**4, num_runs],
                    [[States.HIGH, States.HIGH], [0, 2.5*10**6], 26, 1*10**4, num_runs],
                    [[States.ZERO, States.HIGH],[0, 1.2*10**6], 26, 2*10**4, num_runs],
                    [[States.ZERO, States.ZERO],[0, 1.2*10**6], 26, 2*10**4, num_runs]
                    ])

contrast = .15
exp_count_rate = 17 # kcps
readout_window = 510 # ns

t1_exp_times(t1_exp_array, contrast, exp_count_rate, readout_window)
