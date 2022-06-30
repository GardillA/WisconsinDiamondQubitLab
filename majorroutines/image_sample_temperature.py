# -*- coding: utf-8 -*-
"""
Image the temperature of a sample by zfs thermometry in a raster scan.
Only designed for ensemble.

Created on June 19th, 2022

@author: mccambria
"""


import numpy as np
import utils.tool_belt as tool_belt
import time
import labrad
from majorroutines import pulsed_resonance
from utils.tool_belt import States
import copy
import matplotlib.pyplot as plt
import analysis.temp_from_resonances as temp_from_resonances


# region Functions


def process_resonances(ref_resonances, signal_resonances):

    ref_zfss = [[el[1] - el[0] for el in row] for row in ref_resonances]
    signal_zfss = [[el[1] - el[0] for el in row] for row in signal_resonances]

    ref_temps = [
        [temp_from_resonances(zfs) for zfs in row] for row in ref_zfss
    ]
    signal_temps = [
        [temp_from_resonances(zfs) for zfs in row] for row in signal_zfss
    ]

    ref_zfss = np.array(ref_zfss)
    signal_zfss = np.array(signal_zfss)
    ref_temps = np.array(ref_temps)
    signal_temps = np.array(signal_temps)

    diff_temps = signal_temps - ref_temps

    return diff_temps


# endregion

# region Main


def main(
    nv_sig,
    x_range,
    y_range,
    num_steps,
    apd_indices,
):

    with labrad.connect() as cxn:
        img_array, x_voltages, y_voltages = main_with_cxn(
            cxn,
            nv_sig,
            x_range,
            y_range,
            num_steps,
            apd_indices,
        )

    return img_array, x_voltages, y_voltages


def main_with_cxn(
    cxn,
    nv_sig,
    x_range,
    y_range,
    num_steps,
    apd_indices,
):

    # Some initial setup

    tool_belt.reset_cfm(cxn)

    drift = tool_belt.get_drift()
    coords = nv_sig["coords"]
    adjusted_coords = (np.array(coords) + np.array(drift)).tolist()
    x_center, y_center, z_center = adjusted_coords

    signal_resonances = [
        [
            None,
        ]
        * num_steps
    ] * num_steps
    ref_resonances = [
        [
            None,
        ]
        * num_steps
    ] * num_steps

    freq_range = 0.040
    num_steps = 51
    num_reps = 4e3
    num_runs = 4
    pesr_low_lambda = lambda adj_nv_sig: pulsed_resonance.state(
        adj_nv_sig,
        apd_indices,
        States.LOW,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
    )
    pesr_high_lambda = lambda adj_nv_sig: pulsed_resonance.state(
        adj_nv_sig,
        apd_indices,
        States.HIGH,
        freq_range,
        num_steps,
        num_reps,
        num_runs,
    )

    cxn_power_supply = cxn.power_supply_mp710087

    # Get the voltages for the raster
    x_num_steps = num_steps
    y_num_steps = num_steps
    half_x_range = x_range / 2
    half_y_range = y_range / 2
    x_low = x_center - half_x_range
    x_high = x_center + half_x_range
    y_low = y_center - half_y_range
    y_high = y_center + half_y_range
    x_voltages_1d = np.linspace(x_low, x_high, num_steps)
    y_voltages_1d = np.linspace(y_low, y_high, num_steps)

    # Start rasterin'

    parity = +1  # Determines x scan direction

    for y_ind in range(y_num_steps):
        y_voltage = y_voltages_1d[y_ind]

        for x_ind in range(x_num_steps):
            adj_x_ind = x_ind if parity == +1 else -1 - x_ind
            x_voltage = x_voltages_1d[adj_x_ind]

            adjusted_nv_sig = copy.deepcopy(nv_sig)
            adjusted_nv_sig["coords"] = [x_voltage, y_voltage, z_center]

            cxn_power_supply.output_off()

            time.sleep(11)

            res_low, _ = pesr_low_lambda(adjusted_nv_sig)
            res_high, _ = pesr_high_lambda(adjusted_nv_sig)
            ref_resonances[y_ind][adj_x_ind] = (res_low, res_high)

            cxn_power_supply.output_on()
            cxn_power_supply.set_voltage(1.3)

            time.sleep(10)

            res_low, _ = pesr_low_lambda(adjusted_nv_sig)
            res_high, _ = pesr_high_lambda(adjusted_nv_sig)
            signal_resonances[y_ind][adj_x_ind] = (res_low, res_high)

    # Processing

    diff_temps = process_resonances(ref_resonances, signal_resonances)

    img = plt.imshow(diff_temps)

    # Clean up

    tool_belt.reset_cfm(cxn)
    xy_server = tool_belt.get_xy_server(cxn)
    xy_server.write_xy(x_center, y_center)

    # Save the data

    timestamp = tool_belt.get_time_stamp()
    # print(nv_sig['coords'])
    rawData = {
        "timestamp": timestamp,
        "nv_sig": nv_sig,
        "nv_sig-units": tool_belt.get_nv_sig_units(),
        "drift": drift,
        "x_range": x_range,
        "x_range-units": "V",
        "y_range": y_range,
        "y_range-units": "V",
        "num_steps": num_steps,
        "readout-units": "ns",
        "x_voltages": x_voltages_1d.tolist(),
        "x_voltages-units": "V",
        "y_voltages": y_voltages_1d.tolist(),
        "y_voltages-units": "V",
        "ref_resonances": ref_resonances.astype(float).tolist(),
        "signal_resonances": signal_resonances.astype(float).tolist(),
        "diff_temps": diff_temps.astype(float).tolist(),
        "diff_temps-units": "Kelvin",
    }

    filePath = tool_belt.get_file_path(__file__, timestamp, nv_sig["name"])
    tool_belt.save_raw_data(rawData, filePath)
    tool_belt.save_figure(img, filePath)

    return diff_temps, x_voltages_1d, y_voltages_1d


# endregion


# region Run the file


if __name__ == "__main__":

    pass

# endregion
