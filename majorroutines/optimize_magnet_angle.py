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
import majorroutines.resonance as resonance
from random import shuffle
import copy


# %% Figure functions


def create_fit_figure(splittings, angles, fit_func, popt):
    opti_angle = popt[2] % 180

    fig, ax = plt.subplots(figsize=(8.5, 8.5))

    ax.set_title('ESR Splitting Versus Magnet Angle')
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Splitting (MHz)')
    ax.scatter(angles, splittings, c='r')

    x_vals = numpy.linspace(min(angles), max(angles) + 30, 1000)
    y_vals = fit_func(x_vals, *popt)
    ax.plot(x_vals, y_vals)
    text = ('Optimized Angle: {}'.format('%.1f'%opti_angle))

    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(0.70, 0.90, text, transform=ax.transAxes, fontsize=12,
            verticalalignment="top", bbox=props)

    fig.canvas.draw()
    fig.set_tight_layout(True)
    fig.canvas.flush_events()

    return fig


# %% Other functions


def fit_data(splittings, angles):

    fit_func = AbsCos
    amp = max(splittings)
    phase = angles[numpy.argmax(splittings)] % 180
    offset = 0
    guess_params = [offset, amp, phase]
    bounds = ([0, 0, 0], [numpy.inf, numpy.inf, 180])
    # Check if we have any undefined splittings
    if any(numpy.isnan(splittings)):
        fit_func = None
        popt = None
    else:
        try:
            popt, pcov = curve_fit(fit_func, angles, splittings,
                                   p0=guess_params, bounds=bounds)
        except Exception as e:
            print(e)
            fit_func = None
            popt = None

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
    angle_rad = angle * numpy.pi / 180
    phase_rad = phase * numpy.pi / 180
    return offset + abs(amp * numpy.cos(angle_rad - phase_rad))

def AbsCosNoOff(angle, amp, phase):
    angle_rad = angle * numpy.pi / 180
    phase_rad = phase * numpy.pi / 180
    return abs(amp * numpy.cos(angle_rad - phase_rad))


# %% Main


def main(nv_sig, angle_range, num_angle_steps,
         freq_center, freq_range,
         num_freq_steps, num_freq_reps, num_freq_runs,
         uwave_power, uwave_pulse_dur=None):
    """When you run the file, we'll call into main, which should contain the
    body of the routine.
    """

    tool_belt.check_exp_lock()
    tool_belt.set_exp_lock()
    
    with labrad.connect() as cxn:
        main_with_cxn(cxn, nv_sig, angle_range, num_angle_steps,
                      freq_center, freq_range,
                      num_freq_steps, num_freq_reps, num_freq_runs,
                      uwave_power, uwave_pulse_dur)

def main_with_cxn(cxn, nv_sig, angle_range, num_angle_steps,
                  freq_center, freq_range,
                  num_freq_steps, num_freq_reps, num_freq_runs,
                  uwave_power, uwave_pulse_dur):

    # %% Initial set up here
    
    tool_belt.reset_cfm(cxn)

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
    cwesr = resonance.main_with_cxn

    for ind in angle_inds:

        angle = angles[ind]
        nv_sig_copy['magnet_angle'] = angle

        angle_resonances = (None, None)  # Default to Nones
        if uwave_pulse_dur is not None:
            _, _, angle_resonances = pesr(cxn, nv_sig_copy,
                                    freq_center, freq_range, num_freq_steps,
                                    num_freq_reps, num_freq_runs,
                                    uwave_power, uwave_pulse_dur, standalone_exp = False)
        else:
            angle_resonances = cwesr(cxn, nv_sig_copy,
                                     freq_center, freq_range, num_freq_steps,
                                     num_freq_runs, uwave_power, standalone_exp = False)
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

    fit_func, popt = fit_data(splittings, angles)
    opti_angle = None
    fig = None
    if (fit_func is not None) and (popt is not None):
        fig = create_fit_figure(splittings, angles, fit_func, popt)
        # Find the angle at the peak within [0, 180]
        opti_angle = popt[2] % 180
        print('Optimized angle: {}'.format(opti_angle))
        try:
            rotation_stage_server = tool_belt.get_server_magnet_rotation(cxn)
            rotation_stage_server.set_angle(opti_angle)
        except:
            print("trying to set magnet angle with no rotation stage. check config?")

    # %% Wrap up

    # Set up the raw data dictionary
    raw_data = {'nv_sig': nv_sig,
                'nv_sig-units': tool_belt.get_nv_sig_units(cxn),
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

    tool_belt.reset_cfm(cxn)
    tool_belt.set_exp_unlock()

# %% Run the file


# The __name__ variable will only be '__main__' if you run this file directly.
# This allows a file's functions, classes, etc to be imported without running
# the script that you set up here.
if __name__ == '__main__':

    file = '2023_07_07-12_22_46-WiQD-nv1'
    data = tool_belt.get_raw_data(file)
    splittings = data['splittings']
    print(splittings)

    angle_range = data['angle_range']
    num_angle_steps = data['num_angle_steps']
    angles = numpy.linspace(angle_range[0], angle_range[1], num_angle_steps)

    fit_func, popt = fit_data(splittings, angles)

    opti_angle = None
    fig = None
    if (fit_func is not None) and (popt is not None):
        fig = create_fit_figure(splittings, angles, fit_func, popt)
        # Find the angle at the peak within [0, 180]
        opti_angle = popt[2] % 180
        print('Optimized angle: {}'.format(opti_angle))
