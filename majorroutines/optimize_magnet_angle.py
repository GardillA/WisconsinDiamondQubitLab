# -*- coding: utf-8 -*-
"""Automatically optimize the magnet angle by recording the splittings
at various magnet orientations.

Created on Sun Jun 16 11:38:17 2019

@author: mccambria
"""


# %% Imports


import utils.tool_belt as tool_belt
import labrad
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit
import majorroutines.pulsed_resonance as pulsed_resonance
from random import shuffle
import copy


# %% Figure functions


def create_fit_figure(splittings, angles, fit_func, popt):
    fig, ax = plt.subplots(figsize=(8.5, 8.5))

    ax.set_title('ESR Splitting Versus Magnet Angle')
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Splitting (MHz)')
    ax.scatter(angles, splittings, c='r')

    x_vals = numpy.linspace(0, 180, 1000)
    y_vals = fit_func(x_vals, *popt)
    ax.plot(x_vals, y_vals)

    fig.canvas.draw()
    fig.set_tight_layout(True)
    fig.canvas.flush_events()

    return fig


# %% Other functions


def fit_data(splittings, angles):

    fit_func = AbsCosNoOff
    amp = 200
    phase = 0
    guess_params = [amp, phase]
    # Check if we have any undefined splittings
    if any(numpy.isnan(splittings)):
        fit_func = None
        popt = None
    else:
        popt, pcov = curve_fit(fit_func, angles, splittings, p0=guess_params)

    return fit_func, popt

def clean_up(cxn):

    tool_belt.reset_cfm()

def save_data(name, raw_data, fig):
    """Save the raw data to a txt file as a json object. Save the figures as
    svgs.
    """

    time_stamp = tool_belt.get_time_stamp()

    file_path = tool_belt.get_file_path(__file__, time_stamp, name)

    tool_belt.save_raw_data(raw_data, file_path)

    if fig is not None:
        tool_belt.save_figure(fig, file_path)

def AbsCos(angle, offset, amp, phase):
    return offset + abs(amp * numpy.cos(angle * numpy.pi / 180 + phase * numpy.pi / 180))

def AbsCosNoOff(angle, amp, phase):
    return abs(amp * numpy.cos(angle * numpy.pi / 180 + phase * numpy.pi / 180))


# %% Main


def main(nv_sig, apd_indices, angle_range, num_angle_steps,
         freq_center, freq_range,
         num_freq_steps, num_freq_reps, num_freq_runs,
         uwave_power, uwave_pulse_dur):
    """When you run the file, we'll call into main, which should contain the
    body of the routine.
    """

    with labrad.connect() as cxn:
        main_with_cxn(cxn, nv_sig, apd_indices, angle_range, num_angle_steps,
                      freq_center, freq_range,
                      num_freq_steps, num_freq_reps, num_freq_runs,
                      uwave_power, uwave_pulse_dur)

def main_with_cxn(cxn, nv_sig, apd_indices, angle_range, num_angle_steps,
                  freq_center, freq_range,
                  num_freq_steps, num_freq_reps, num_freq_runs,
                  uwave_power, uwave_pulse_dur):

    # %% Initial set up here

    angles = numpy.linspace(angle_range[0], angle_range[1], num_angle_steps)
    angle_inds = numpy.linspace(0, num_angle_steps-1, num_angle_steps, dtype=int)
    shuffle(angle_inds)
    resonances = numpy.empty((num_angle_steps, 2))  # Expect 2 resonances
    resonances[:] = numpy.nan
    splittings = numpy.empty(num_angle_steps)
    splittings[:] = numpy.nan

    # %% Collect the data

    nv_sig_copy = copy.deepcopy(nv_sig)
    pesr = pulsed_resonance.main_with_cxn

    for ind in angle_inds:

        angle = angles[ind]
        nv_sig_copy['magnet_angle'] = angle

        angle_resonances = pesr(cxn, nv_sig_copy, apd_indices,
                                freq_center, freq_range,
                                num_freq_steps, num_freq_reps, num_freq_runs,
                                uwave_power, uwave_pulse_dur)
        resonances[ind, :] = angle_resonances
        if all(angle_resonances):
            # We got two resonances so take the difference
            splittings[ind] = (angle_resonances[1] - angle_resonances[0]) * 1000
        elif any(angle_resonances):
            # We got one resonance so consider this a splitting of 0
            splittings[ind] = 0
        else:
            # We failed to find any resonances
            splittings[ind] = None

    # %% Analyze the data

    fit_func = AbsCosNoOff
    amp = 200
    phase = 50
    guess_params = [amp, phase]
    # Check if we have any undefined splittings
    if any(numpy.isnan(splittings)):
        opti_angle = None
        fig = None
    else:
        try:
            popt, pcov = curve_fit(fit_func, angles, splittings, p0=guess_params)
            # Find the angle at the peak within [0, 360]
            opti_angle = (-popt[1]) % 360
            fig = create_fit_figure(splittings, angles, fit_func, popt)
        except Exception:
            opti_angle = None

    # %% Wrap up

    if opti_angle is not None:
        cxn.rotation_stage_ell18k.set_angle(opti_angle)
        print('Optimized angle: {}'.format(opti_angle))

    # Set up the raw data dictionary
    raw_data = {'nv_sig': nv_sig,
                'nv_sig-units': tool_belt.get_nv_sig_units(),
                'apd_indices': apd_indices,
                'angle_range': angle_range,
                'angle_range-units': 'deg',
                'num_angle_steps': num_angle_steps,
                'freq_center': freq_center,
                'freq_center-units': 'GHz',
                'freq_range': freq_range,
                'freq_range-units': 'GHz',
                'num_freq_steps': num_freq_steps,
                'num_freq_runs': num_freq_runs,
                'uwave_power': uwave_power,
                'uwave_power-units': 'dBm',
                'resonances': resonances.tolist(),
                'resonances-units': 'GHz',
                'splittings': splittings.tolist(),
                'splittings-units': 'MHz',
                'opti_angle': opti_angle,
                'opti_angle-units': 'deg'}

    # Save the data and the figures from this run
    save_data(nv_sig['name'], raw_data, fig)


# %% Run the file


# The __name__ variable will only be '__main__' if you run this file directly.
# This allows a file's functions, classes, etc to be imported without running
# the script that you set up here.
if __name__ == '__main__':

    # nv25_2019_07_25
    # splittings = [31.4, 11.0, 43.0, 16.3, 33.0, 42.8]
    # angles = [90, 120, 60, 150, 0, 30]
    
    # nv27_2019_07_25
    splittings = [228.9, None, None, 83.2, 44.5]
    angles = [0.0, 60.0, 120, 150, 90]
    
    splittings = [228.9, 83.2, 44.5, 231.4]
    angles = [0.0, 150, 90, 30]

    fit_func, popt = fit_data(splittings, angles)

    opti_angle = None
    fig = None
    if (fit_func is not None) and (popt is not None):
        fig = create_fit_figure(splittings, angles, fit_func, popt)
        # Find the angle at the peak within [0, 180]
        opti_angle = (-popt[1]) % 180
        print('Optimized angle: {}'.format(opti_angle))

