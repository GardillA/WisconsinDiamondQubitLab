# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:00:15 2019

@author: matth
"""


# %% Imports


import numpy
import matplotlib
import matplotlib.pyplot as plt
import utils.tool_belt as tool_belt
import utils.common as common
import json
from mpl_toolkits.axes_grid1.anchored_artists import (
    AnchoredSizeBar as scale_bar,
)
from colorutils import Color
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

ms = 7
lw = 1.75


# %% Functions


def process_raw_data(data, ref_range):
    """Pull the relaxation signal and ste out of the raw data."""

    num_runs = data["num_runs"]
    num_steps = data["num_steps"]
    sig_counts = numpy.array(data["sig_counts"])
    ref_counts = numpy.array(data["ref_counts"])
    time_range = numpy.array(data["relaxation_time_range"])

    # Calculate time arrays in ms
    min_time, max_time = time_range / 10 ** 6
    times = numpy.linspace(min_time, max_time, num=num_steps)

    # Calculate the average signal counts over the runs, and ste
    avg_sig_counts = numpy.average(sig_counts[::], axis=0)
    ste_sig_counts = numpy.std(sig_counts[::], axis=0, ddof=1) / numpy.sqrt(
        num_runs
    )

    # Assume reference is constant and can be approximated to one value
    avg_ref = numpy.average(ref_counts[::])

    # Divide signal by reference to get normalized counts and st error
    norm_avg_sig = avg_sig_counts / avg_ref
    norm_avg_sig_ste = ste_sig_counts / avg_ref

    # Normalize to the reference range
    diff = ref_range[1] - ref_range[0]
    norm_avg_sig = (norm_avg_sig - ref_range[0]) / diff

    return norm_avg_sig, norm_avg_sig_ste, times


def relaxation_zero_func(t, gamma, omega, infid):

    return (1 / 3) + (2 / 3) * numpy.exp(-3 * omega * t)


def relaxation_high_func(t, gamma, omega, infid):

    # Times are in ms, but rates are in s^-1
    gamma /= 1000
    omega /= 1000

    first_term = (1 / 3) + (1 / 2) * ((1 - infid) ** 2) * numpy.exp(
        -(2 * gamma + omega) * t
    )
    second_term = (
        (-1 / 2) * (infid - (1 / 3)) * numpy.exp(-3 * omega * t) * (1 - infid)
    )
    third_term = (infid - (1 / 3)) * numpy.exp(-3 * omega * t) * infid
    return first_term + second_term + third_term


def get_ref_range(ref_files, path_from_nvdata):

    ref_range = []
    for el in ref_files:
        data = tool_belt.get_raw_data(el, path_from_nvdata)
        sig_counts = numpy.array(data["sig_counts"])
        ref_counts = numpy.array(data["ref_counts"])
        avg_ref = numpy.average(ref_counts[::])
        avg_sig_counts = numpy.average(sig_counts[::], axis=0)
        ref_range.append(avg_sig_counts[0] / avg_ref)
    return ref_range


def exp_eq(t, rate, amp):
    return amp * numpy.exp(-rate * t)


def exp_eq_offset(t, rate, amp, offset):
    return amp * numpy.exp(-rate * t) + offset


# %% Main


def main(data_sets, image_files):

    # fig, axes_pack = plt.subplots(1,2, figsize=(10,5))
    fig = plt.figure(figsize=(6.5, 7))
    grid_columns = 20
    half_grid_columns = grid_columns // 2
    gs = gridspec.GridSpec(2, grid_columns, height_ratios=(1, 1))

    first_row_sep_ind = 9

    # %% Level structure

    # Add a new axes, make it invisible, steal its rect
    ax = fig.add_subplot(gs[0, 0:first_row_sep_ind])
    ax.set_axis_off()
    ax.text(
        0,  # -0.43,
        0.95,
        "(a)",
        transform=ax.transAxes,
        color="black",
        fontsize=18,
    )

    if False:

        ax = plt.Axes(fig, [0.0, 0.2, 0.5, 0.43])
        ax.set_axis_off()
        fig.add_axes(ax)

        level_structure_file = "/run/media/mccambria/Barracuda/acer2-2019-2020/lab/bulk_dq_relaxation/figures_revision2/main1/level_structure.png"
        img = mpimg.imread(level_structure_file)
        _ = ax.imshow(img)

    # l, b, w, h = ax.get_position().bounds
    # ax.set_position([l, b -1.0, w, h])
    # ax.set_axis_off()
    # ax.axis('off')

    # %% Relaxation out of plots

    temps = [el["temp"] for el in data_sets]

    min_temp = min(temps)
    max_temp = max(temps)
    temp_range = max_temp - min_temp
    normalized_temps = [(val - min_temp) / temp_range for val in temps]
    cmap = matplotlib.cm.get_cmap("coolwarm")
    # Trim the alpha value and convert from 0:1 to 0:255
    colors_rgb = [
        [255 * el for el in cmap(val)[0:3]] for val in normalized_temps
    ]
    colors_Color = [Color(tuple(el)) for el in colors_rgb]
    colors_hex = [val.hex for val in colors_Color]
    colors_hsv = [val.hsv for val in colors_Color]
    facecolors_hsv = [(el[0], 0.3 * el[1], 1.2 * el[2]) for el in colors_hsv]
    facecolors_Color = [Color(hsv=el) for el in facecolors_hsv]
    facecolors_hex = [val.hex for val in facecolors_Color]

    ax = fig.add_subplot(gs[0, first_row_sep_ind:])
    # ax = axes_pack[1]

    ax.set_xlabel(r"Wait time $\tau$ (ms)")
    ax.set_ylabel("Normalized fluorescence")
    ax.set_xticks([0, 5, 10, 15])

    times = [0.0, 15.0]

    # Plot decay curves
    for ind in range(len(data_sets)):

        data_set = data_sets[ind]
        path_from_nvdata = data_set["path_from_nvdata"]
        raw_decay = tool_belt.get_raw_data(
            data_set["decay_file"], path_from_nvdata
        )
        ref_files = [data_set["zero_ref_file"], data_set["unity_ref_file"]]
        # ref_range = get_ref_range(ref_files, path_from_nvdata)
        # MCC remove this after single NV data
        ref_range = [0.65, 0.99]

        signal_decay, ste_decay, times_decay = process_raw_data(
            raw_decay, ref_range
        )
        smooth_t = numpy.linspace(times[0], times[-1], 1000)
        fit_decay = relaxation_high_func(
            smooth_t, data_set["gamma"], data_set["Omega"], 0.0
        )

        color = colors_hex[ind]
        facecolor = facecolors_hex[ind]
        label = "Relaxation \nout of {}".format(r"$\ket{0}$")
        patch = mlines.Line2D(
            [],
            [],
            label=label,
            linewidth=lw,
            marker="o",
            color=color,
            markerfacecolor=facecolor,
            markersize=ms,
        )
        ax.plot(smooth_t, fit_decay, color=color, linewidth=lw)
        try:
            times_15 = numpy.where(times_decay > 15.0)[0][0]
        except:
            times_15 = None
        ax.scatter(
            times_decay[:times_15],
            signal_decay[:times_15],
            label="{} K".format(data_set["temp"]),
            zorder=5,
            marker="o",
            color=color,
            facecolor=facecolor,
            s=ms ** 2,
        )

    ax.legend()
    fig.text(
        -0.18, 0.95, "(b)", transform=ax.transAxes, color="black", fontsize=18
    )

    # l, b, w, h = ax.get_position().bounds
    # ax.set_position([l - 0.01, b, w, h])

    # %% Sample plots

    fig_labels = [r"(c)", r"(d)"]
    sample_labels = ["Sample A", "Sample B", "Sample C"]

    for ind in range(len(image_files)):

        ax = fig.add_subplot(
            gs[1, half_grid_columns * ind : half_grid_columns * (ind + 1)]
        )
        # ax.set_axis_off()
        name = image_files[ind]

        # Load the data from the file
        data = tool_belt.get_raw_data(name)

        # Build the image array from the data
        # Not sure why we're doing it this way...
        img_array = []
        try:
            file_img_array = data["img_array"]
        except:
            file_img_array = data["imgArray"]
        for line in file_img_array:
            img_array.append(line)

        # Get the readout in s
        readout = float(data["readout"]) / 10 ** 9

        try:
            xScanRange = data["x_range"]
            yScanRange = data["y_range"]
        except:
            xScanRange = data["xScanRange"]
            yScanRange = data["yScanRange"]

        kcps_array = (numpy.array(img_array) / 1000) / readout

        # Scaling
        scale = 35  # 35  # galvo scaling in microns / volt, MCC correct?
        num_steps = kcps_array.shape[0]
        v_resolution = xScanRange / num_steps  # resolution in volts / pixel
        resolution = v_resolution * scale  # resolution in microns / pixel
        px_per_micron = 1 / resolution

        # Plot 5 um out from center in any direction
        center = [num_steps // 2, num_steps // 2]
        clip_range = 5 * px_per_micron
        x_clip = [center[0] - clip_range, center[0] + clip_range]
        x_clip = [int(el) for el in x_clip]
        y_clip = [center[1] - clip_range, center[1] + clip_range]
        y_clip = [int(el) for el in y_clip]
        neg_test = [val < 0 for val in x_clip + y_clip]
        if True in neg_test:
            raise ValueError("Negative value encountered in image coordinates")
        # print((x_clip[1] - x_clip[0]) * v_resolution)
        clip_array = kcps_array[x_clip[0] : x_clip[1], y_clip[0] : y_clip[1]]
        print(numpy.array(kcps_array).shape)
        print(numpy.array(clip_array).shape)
        img = ax.imshow(clip_array, cmap="inferno", interpolation="none")
        ax.set_axis_off()
        # plt.axis("off")

        # Scale bar
        trans = ax.transData
        bar_text = r"2 $\upmu$m "  # Insufficient left padding...
        bar = scale_bar(
            trans,
            2 * px_per_micron,
            bar_text,
            "upper right",
            size_vertical=int(num_steps / 100),
            pad=0.25,
            borderpad=0.5,
            sep=4.0,
            # frameon=False, color='white',
        )
        ax.add_artist(bar)

        # Labels
        fig_label = fig_labels[ind]
        ax.text(
            0.035,
            0.88,
            fig_label,
            transform=ax.transAxes,
            color="white",
            fontsize=20,
        )
        sample_label = sample_labels[ind]
        ax.text(
            0.035,
            0.07,
            sample_label,
            transform=ax.transAxes,
            color="white",
            fontsize=18,
        )
        # cbar = plt.colorbar(img)

    # %% Wrap up

    shift = 0.103
    gs.tight_layout(fig, pad=0.3, w_pad=-2.50)
    # gs.tight_layout(fig, pad=0.4, h_pad=0.5, w_pad=0.5, rect=[0, 0, 1, 1])
    # fig.tight_layout(pad=0.5)
    # fig.tight_layout()
    # plt.margins(0, 0)
    # fig.subplots_adjust(hspace=0.5, wspace=0.5)


# %% Run


if __name__ == "__main__":

    tool_belt.init_matplotlib()
    # plt.rcParams.update({'font.size': 18})  # Increase font size
    matplotlib.rcParams["axes.linewidth"] = 1.0

    # nvdata_dir = common.get_nvdata_dir()

    # Data set a: room temp (300 K)
    data_a = {
        "temp": 300,
        "path_from_nvdata": "pc_hahn/branch_cryo-setup/t1_interleave_knill/data_collections/hopper-nv1_2021_03_16-300K/",
        "decay_file": "2021_05_11-21_11_55-hopper-nv1_2021_03_16",  # high to high
        "unity_ref_file": "2021_05_11-21_11_47-hopper-nv1_2021_03_16",  # zero to zero
        "zero_ref_file": "2021_05_11-21_11_53-hopper-nv1_2021_03_16",  # zero to high
        "Omega": 59.87,
        "gamma": 131.57,
    }

    # Data set b: 250 K data
    data_b = {
        "temp": 250,
        "path_from_nvdata": "pc_hahn/branch_cryo-setup/t1_interleave_knill/data_collections/hopper-nv1_2021_03_16-250K/",
        "decay_file": "2021_05_12-18_41_17-hopper-nv1_2021_03_16",  # high to high
        "unity_ref_file": "2021_05_12-18_41_10-hopper-nv1_2021_03_16",  # zero to zero
        "zero_ref_file": "2021_05_12-18_41_12-hopper-nv1_2021_03_16",  # zero to high
        "Omega": 28.53,
        "gamma": 71.51,
    }

    decay_data_sets = [data_a, data_b]

    # 60 x 60, 0.2 x 0.2
    sample_image_files = [
        "2022_01_20-22_48_34-wu-nv6_2021_12_25",  # Sample A, Wu
        "2021_05_19-11_44_47-hopper-nv1_2021_03_16",  # Sample B, Hopper
    ]

    main(decay_data_sets, sample_image_files)

    plt.show(block=True)
