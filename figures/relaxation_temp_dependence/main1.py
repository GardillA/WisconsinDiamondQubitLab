# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:00:15 2019

@author: matth
"""


# %% Imports


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import utils.tool_belt as tool_belt
import majorroutines.rabi as rabi
import utils.common as common
import json
from mpl_toolkits.axes_grid1.anchored_artists import (
    AnchoredSizeBar as scale_bar,
)
from scipy.optimize import curve_fit
from colorutils import Color
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
from figures.relaxation_temp_dependence.temp_dependence_fitting import (
    omega_calc,
    gamma_calc,
)

ms = 7
lw = 1.75


# %% Functions


def get_data_decay(folder):
    """Return the data to plot."""

    start_run = None
    stop_run = None
    start_time_ind = None
    end_time_ind = None

    num_runs = data["num_runs"]
    num_steps = data["num_steps"]
    sig_counts = np.array(data["sig_counts"])
    ref_counts = np.array(data["ref_counts"])
    time_range = np.array(data["relaxation_time_range"])

    # _, ax = plt.subplots()
    # counts_flatten = sig_counts[:, -3]
    # # counts_flatten = ref_counts.flatten()
    # bins = np.arange(0, max(counts_flatten) + 1, 1)
    # ax.hist(counts_flatten, bins, density=True)

    # Calculate time arrays in ms
    min_time, max_time = time_range / 10 ** 6
    times = np.linspace(min_time, max_time, num=num_steps)
    # times[0] = 0.5

    # Calculate the average signal counts over the runs, and ste
    avg_sig_counts = np.average(sig_counts[start_run:stop_run, :], axis=0)
    avg_ref_counts = np.average(ref_counts[start_run:stop_run, :], axis=0)
    std_sig_counts = np.std(
        sig_counts[start_run:stop_run, :],
        axis=0,
        ddof=1,
    )
    std_ref_counts = np.std(
        ref_counts[start_run:stop_run, :],
        axis=0,
        ddof=1,
    )
    # std_sig_counts = np.sqrt(avg_sig_counts)
    ste_sig_counts = std_sig_counts / np.sqrt(num_runs)
    ste_ref_counts = std_ref_counts / np.sqrt(num_runs)
    # print(ste_sig_counts)

    single_ref = False
    if single_ref:
        # Assume reference is constant and can be approximated to one value
        avg_ref = np.average(ref_counts[start_run:stop_run, :])
        # Divide signal by reference to get normalized counts and st error
        norm_avg_sig = avg_sig_counts / avg_ref
        norm_avg_sig_ste = ste_sig_counts / avg_ref
    else:
        # Divide signal by reference to get normalized counts and st error
        norm_avg_sig = avg_sig_counts / avg_ref_counts
        norm_avg_sig_ste = norm_avg_sig * np.sqrt(
            (ste_sig_counts / avg_sig_counts) ** 2
            + (ste_ref_counts / avg_ref_counts) ** 2
        )

    return (
        norm_avg_sig[start_time_ind:end_time_ind],
        norm_avg_sig_ste[start_time_ind:end_time_ind],
        times[start_time_ind:end_time_ind],
    )


def zero_to_one_threshold(val):
    if val < 0:
        return 0
    elif val > 1:
        return 1
    else:
        return val


# %% Main


def main(data_sets, dosave=False, draft_version=True):

    nvdata_dir = common.get_nvdata_dir()

    # fig, axes_pack = plt.subplots(1,2, figsize=(10,5))
    fig = plt.figure(figsize=(6.5, 7.5))
    grid_columns = 30
    half_grid_columns = grid_columns // 2
    gs = gridspec.GridSpec(2, grid_columns, height_ratios=(1, 1))

    first_row_sep_ind = 15

    # %% Level structure

    # Add a new axes, make it invisible, steal its rect
    ax = fig.add_subplot(gs[0, 0:first_row_sep_ind])
    ax.set_axis_off()
    ax.text(
        0,
        0.95,
        "(a)",
        transform=ax.transAxes,
        color="black",
        fontsize=18,
    )

    draft_version = True
    # draft_version = False
    if draft_version:
        ax = plt.Axes(fig, [-0.05, 0.5, 0.5, 0.43])
        ax.set_axis_off()
        fig.add_axes(ax)
        level_structure_file = (
            nvdata_dir
            / "paper_materials/relaxation_temp_dependence/figures/level_structure.png"
        )
        img = mpimg.imread(level_structure_file)
        _ = ax.imshow(img)

    # %% Gamma subtraction curve plots

    temps = [round(el["temp"]) for el in data_sets]

    continuous_colormap = False
    if continuous_colormap:
        min_temp = min(temps)
        max_temp = max(temps)
        temp_range = max_temp - min_temp
        normalized_temps = [(val - min_temp) / temp_range for val in temps]
        # adjusted_temps = [normalized_temps[0], ]
        cmap = matplotlib.cm.get_cmap("coolwarm")
        colors_cmap = [cmap(val) for val in normalized_temps]
    else:
        set1 = matplotlib.cm.get_cmap("Set1").colors
        set2 = matplotlib.cm.get_cmap("Dark2").colors
        colors_cmap = [set1[6], set1[0], set2[5], set1[2], set1[1]]

    # Trim the alpha value and convert from 0:1 to 0:255
    colors_rgb = [[255 * val for val in el[0:3]] for el in colors_cmap]
    colors_Color = [Color(tuple(el)) for el in colors_rgb]
    colors_hex = [val.hex for val in colors_Color]
    colors_hsv = [val.hsv for val in colors_Color]
    facecolors_hsv = [(el[0], 0.3 * el[1], 1.2 * el[2]) for el in colors_hsv]
    # Threshold to make sure these are valid colors
    facecolors_hsv = [
        (el[0], zero_to_one_threshold(el[1]), zero_to_one_threshold(el[2]))
        for el in facecolors_hsv
    ]
    facecolors_Color = [Color(hsv=el) for el in facecolors_hsv]
    facecolors_hex = [val.hex for val in facecolors_Color]

    ax = fig.add_subplot(gs[0, first_row_sep_ind:])
    # ax = axes_pack[1]
    l, b, w, h = ax.get_position().bounds
    shift = 0.02
    ax.set_position([l + shift, b, w - shift, h])

    ax.set_xlabel(r"Wait time $\tau$ (ms)")
    ax.set_ylabel(r"P_{+1,+1}(\tau) - P_{+1,-1}(\tau)")

    min_time = 0.0
    max_time = 15.0
    xtick_step = 5
    # max_time = 12.5
    # max_time = 9
    # xtick_step = 4
    times = [min_time, max_time]
    ax.set_xticks(np.arange(min_time, max_time + xtick_step, xtick_step))

    # Plot decay curves
    for ind in range(len(data_sets)):

        data_set = data_sets[ind]
        color = colors_hex[ind]
        facecolor = facecolors_hex[ind]
        temp = round(data_set["temp"])
        gamma = data_set["gamma"]
        Omega = data_set["Omega"]

        # Plot the fit/predicted curves
        if (gamma is None) and (Omega is None):
            # MCC make sure these values are up to date
            gamma = gamma_calc(temp)
            Omega = omega_calc(temp)
            smooth_t = np.linspace(times[0], 1.1 * times[-1], 1000)
            fit_decay = np.exp(-(2 * gamma + Omega) * smooth_t)
            ax.plot(smooth_t, fit_decay, color=color, linewidth=lw)

        if data_set["skip"]:
            continue

        folder = data_set["folder"]
        data_decay, ste_decay, times_decay = get_data_decay(folder)

        # Clip anything beyond the max time
        try:
            times_clip = np.where(times_decay > max_time)[0][0]
        except:
            times_clip = None
        times_decay = times_decay[:times_clip]
        signal_decay = signal_decay[:times_clip]
        ste_decay = ste_decay[:times_clip]

        plot_errors = False
        if plot_errors:
            ax.errorbar(
                times_decay,
                data_decay,
                yerr=np.array(ste_decay),
                label="{} K".format(temp),
                zorder=5,
                marker="o",
                color=color,
                markerfacecolor=facecolor,
                ms=ms,
                linestyle="",
            )
        else:
            ax.scatter(
                times_decay,
                data_decay,
                label="{} K".format(temp),
                zorder=5,
                marker="o",
                color=color,
                facecolor=facecolor,
                s=ms ** 2,
            )

    ax.legend(handlelength=5)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
    fig.text(
        -0.21, 0.95, "(b)", transform=ax.transAxes, color="black", fontsize=18
    )
    x_buffer = 0.02 * max_time
    ax.set_xlim([-x_buffer, max_time + x_buffer])
    # ax.set_ylim([0.3, 1.06])
    # ax.set_ylim([0.009, 1.1])
    ax.set_ylim([0.05, 1.1])
    ax.set_yscale("log")

    # %% Experimental layout

    # Add a new axes, make it invisible, steal its rect
    ax = fig.add_subplot(gs[1, :])
    ax.set_axis_off()
    fig.text(
        0,
        0.95,
        "(c)",
        transform=ax.transAxes,
        color="black",
        fontsize=18,
    )

    if draft_version:
        ax.set_axis_off()
        fig.add_axes(ax)
        level_structure_file = (
            nvdata_dir
            / "paper_materials/relaxation_temp_dependence/figures/experimental_layout_simplified.png"
        )
        img = mpimg.imread(level_structure_file)
        _ = ax.imshow(img)

    # %% Wrap up

    shift = 0.103
    gs.tight_layout(fig, pad=0.3, w_pad=-2.50)
    # gs.tight_layout(fig, pad=0.3, w_pad=0)
    # gs.tight_layout(fig, pad=0.4, h_pad=0.5, w_pad=0.5, rect=[0, 0, 1, 1])
    # fig.tight_layout(pad=0.5)
    # fig.tight_layout()
    # plt.margins(0, 0)
    # fig.subplots_adjust(hspace=0.5, wspace=0.5)

    if dosave:
        file_path = str(
            nvdata_dir
            / "paper_materials/relaxation_temp_dependence/figures/main1.eps"
        )
        fig.savefig(file_path, dpi=500)


# %% Run


if __name__ == "__main__":

    tool_belt.init_matplotlib()
    # plt.rcParams.update({'font.size': 18})  # Increase font size
    matplotlib.rcParams["axes.linewidth"] = 1.0

    decay_data_sets = [
        {
            "temp": 380.168,
            "skip": False,
            "folder": "pc_hahn/branch_master/t1_interleave_knill/data_collections/hopper-search-400K",
            "Omega": None,
            "gamma": None,
        },
        {
            "temp": 337.584,
            "skip": False,
            "folder": "pc_hahn/branch_time-tagger-speedup/t1_interleave_knill/data_collections/hopper-search-350K",
            "Omega": None,
            "gamma": None,
        },
        {
            "temp": 295,
            "skip": False,
            "folder": "pc_hahn/branch_cryo-setup/t1_interleave_knill/hopper-nv1_2021_03_16-300K",
            "Omega": None,
            "gamma": None,
        },
        {
            "temp": 250,
            "skip": False,
            "folder": "pc_hahn/branch_cryo-setup/t1_interleave_knill/hopper-nv1_2021_03_16-250K",
            "Omega": None,
            "gamma": None,
        },
        {
            "temp": 200,
            "skip": False,
            "folder": "pc_hahn/branch_cryo-setup/t1_interleave_knill/hopper-nv1_2021_03_16-200K-gamma_minus_1",
            "Omega": None,
            "gamma": None,
        },
    ]

    main(decay_data_sets, dosave=False, draft_version=True)

    plt.show(block=True)
