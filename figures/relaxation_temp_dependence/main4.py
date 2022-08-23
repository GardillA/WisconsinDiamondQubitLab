import errno
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.patches as patches
import matplotlib.lines as mlines
from scipy.optimize import curve_fit
import pandas as pd
import utils.tool_belt as tool_belt
from utils.kplotlib import color_mpl_to_color_hex, lighten_color_hex
import utils.common as common
from scipy.odr import ODR, Model, RealData
import sys
from pathlib import Path
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
import temp_dependence_fitting
from temp_dependence_fitting import qutrit_color, qubit_color
import csv
import utils.kplotlib as kpl

marker_size = 7
line_width = 1.5
# line_width = 2.5
marker_edge_width = line_width


def round_base_2(val):
    power = round(np.log2(val))
    rounded_val = 2 ** power
    return rounded_val


def bar_gill_replot(file_name, path):

    data_points = []
    with open(path / file_name, newline="") as f:
        raw_data = csv.reader(f)
        prev_point_ind = -1
        new_point = None
        header = True
        for row in raw_data:
            if header:
                header = False
                continue
            point_ind = int(row[3])
            if point_ind != prev_point_ind:
                prev_point_ind = point_ind
                if new_point is not None:
                    data_points.append(new_point)
                new_point = {
                    "temp": int(row[0]),
                    "num_pulses": round_base_2(float(row[1])),
                }
            row_type = row[4].strip()
            val = float(row[2])
            new_point[row_type] = val

    for point in data_points:
        T2 = point["main"]
        if ("ste_above" in point) and ("ste_below" in point):
            avg_ste = (
                (point["ste_above"] - T2) + (T2 - point["ste_below"])
            ) / 2
            point["ste"] = avg_ste
        elif "ste_above" in point:
            point["ste"] = point["ste_above"] - T2
        elif "ste_below" in point:
            point["ste"] = T2 - point["ste_below"]
        else:
            point["ste"] = None

    colors = {
        300: "blue",
        240: "green",
        190: "purple",
        160: "cyan",
        120: "red",
        77: "yellow",
    }
    fig, ax = plt.subplots(figsize=[6.5, 5.0])
    for point in data_points:
        ax.errorbar(
            point["num_pulses"],
            point["main"],
            point["ste"],
            color=colors[point["temp"]],
            marker="o",
        )

    ax.set_yscale("log")
    ax.set_xscale("log")
    fig.tight_layout()


def main(
    file_name,
    path,
    plot_type,
    rates_to_plot,
    temp_range,
    y_range,
    xscale,
    yscale,
    dosave=False,
):

    # bar_gill_label = "[10]"
    # herbschleb_label = "[11]"
    # abobeih_label = "[30]"
    # bar_gill_label = "[2]"
    # herbschleb_label = "[3]"
    # abobeih_label = "[4]"
    bar_gill_label = "Bar-Gill"
    herbschleb_label = "Herbschleb"
    abobeih_label = "Abobeih"
    # fmt: off
    data_points = [
        #
        #
        {"val": 580e-3, "err": 210e-3, "temp": 77, "author": "Bar-Gill", "label": bar_gill_label},
        {"val": 152e-3, "err": 52e-3, "temp": 120, "author": "Bar-Gill", "label": bar_gill_label},
        {"val": 39.8e-3, "err": 7.7e-3, "temp": 160, "author": "Bar-Gill", "label": bar_gill_label},
        {"val": 17.3e-3, "err": 4.3e-3, "temp": 190, "author": "Bar-Gill", "label": bar_gill_label},
        {"val": 5.92e-3, "err": 1.23e-3, "temp": 240, "author": "Bar-Gill", "label": bar_gill_label},
        # {"val": 3.34e-3, "err": 0.41e-3, "temp": 300, "author": "Bar-Gill 2013"},
        #
        # Spin echo
        # {"val": 183.83e-6, "err": 13.0e-6, "temp": 300, "author": "Lin"},
        # {"val": 158.15e-6, "err": 10.9e-6, "temp": 350, "author": "Lin"},
        # {"val": 125.50e-6, "err": 7.61e-6, "temp": 400, "author": "Lin"},
        # {"val": 80.480e-6, "err": 6.02e-6, "temp": 450, "author": "Lin"},
        # {"val": 59.239e-6, "err": 5.07e-6, "temp": 500, "author": "Lin"},
        # {"val": 38.315e-6, "err": 4.12e-6, "temp": 550, "author": "Lin"},
        # {"val": 30.389e-6, "err": 3.80e-6, "temp": 600, "author": "Lin"},
        #
        # Also report gamma and Omega at room temps
        {"val": 3.3e-3, "err": None, "temp": 300, "author": "Herbschleb", "label": herbschleb_label},
        #
        # Record, T1 exceeds expected value from one-phonon calculations
        {"val": 1.58, "err": 0.07, "temp": 3.7, "author": "Abobeih", "label": abobeih_label},
        #
        # 
        # {"val": 2.193e-3, "err": None, "temp": 300, "author": "Pham"},
        #
        # Isotopically purified, just spin echo
        # {"val": 1.82e-3, "err": 0.16e-3, "temp": 300, "author": "Balasubramanian"},
        #
        # Original DD?
        # {"val": 88e-6, "err": None, "temp": 300, "author": "de Lange"},
        #
        # 
        # {"val": 1.6e-3, "err": None, "temp": 300, "author": "Ryan"},
        #
        # 
        # {"val": 2.44e-3, "err": 0.44e-3, "temp": 300, "author": "Naydenov"},
    ]
    # fmt: on

    # Sekiguchi Dynamical Decoupling of a Geometric Qubit
    # Optimizing a dynamical decoupling protocol for solid-state electronic spin ensembles in diamon
    # Robust Quantum-Network Memory Using Decoherence-Protected Subspaces of Nuclear Spins
    # Randomization of Pulse Phases for Unambiguous and Robust Quantum Sensing, Why not try T2 limits?
    # Robust quantum control for the manipulation of solid-state spins, Likewise

    ret_vals = temp_dependence_fitting.main(
        file_name,
        path,
        plot_type,
        rates_to_plot,
        temp_range,
        y_range,
        xscale,
        yscale,
        dosave=False,
    )

    if plot_type == "T2_max_supp":
        fig, ax1, ax2, leg1, T2_max_qubit_hopper_temp = ret_vals
    else:
        fig, ax1, leg1, T2_max_qubit_hopper_temp = ret_vals

    # colors = {
    #     "Bar-Gill 2013": "green",
    #     "Lin": "red",
    #     "Abobeih 2018": "purple",
    #     "Balasubramanian": "orange",
    #     "Pham": "blue",
    # }Gill
    # markers = [
    #     "o",
    #     "^",
    #     "s",
    #     "X",
    #     "D",
    #     "H",
    # ]
    markers = {
        "Bar-Gill": "o",
        "Lin": "H",
        "Abobeih": "^",
        "Pham": "s",
        "Herbschleb": "D",
        "Balasubramanian": "v",
        "de Lange": "P",
        "Ryan": "p",
        "Naydenov": "d",
    }

    if plot_type == "T2_max_supp":
        sub_plot_types = ["T2_max", "T2_frac"]
    else:
        sub_plot_types = [plot_type]

    ms = marker_size ** 2
    used_authors = []
    ind = 0
    for sub_plot_type in sub_plot_types:
        ind += 1
        ax = eval(f"ax{ind}")
        for point in data_points:
            temp = point["temp"]
            T2 = point["val"]
            frac = T2 / T2_max_qubit_hopper_temp(temp)
            if sub_plot_type == "T2_max":
                val = T2
            elif sub_plot_type == "T2_frac":
                val = frac
            # err = point["err"]
            err = None
            author = point["author"]
            # color = colors[author]
            marker = markers[author]
            label = None
            if author not in used_authors:
                used_authors.append(author)
                label = point["label"]
            ax.errorbar(
                temp,
                val,
                err,
                # color="black",
                # markerfacecolor="gray",
                color=qubit_color,
                markerfacecolor=kpl.lighten_color_hex(qubit_color),
                label=label,
                marker=marker,
                ms=marker_size,
                lw=line_width,
                markeredgewidth=marker_edge_width,
                linestyle="None",
            )

    # Legend without errorbars
    handles, labels = ax1.get_legend_handles_labels()
    errorbar_type = matplotlib.container.ErrorbarContainer
    # handles = [h[0] if isinstance(h, errorbar_type) else h for h in handles]
    handles = [h[0] for h in handles if isinstance(h, errorbar_type)]
    if leg1 is not None:
        labels = labels[2:]
    leg2 = ax1.legend(
        handles,
        labels,
        title="Prior results",
        loc="upper right",
        # bbox_to_anchor=(1.0, 0.82),
        handlelength=1,
    )
    # Add back in original legend
    if leg1 is not None:
        # anchor = leg2.get_bbox_to_anchor()
        # leg1.set_bbox_to_anchor(anchor)
        leg1.set_bbox_to_anchor((0.687, 1.0))
        ax1.add_artist(leg1)

    # ax.legend()


if __name__ == "__main__":

    tool_belt.init_matplotlib()
    # matplotlib.rcParams["axes.linewidth"] = 1.0

    file_name = "compiled_data"
    home = common.get_nvdata_dir()
    path = home / "paper_materials/relaxation_temp_dependence"

    ### Main

    # plot_type = "T2_max"
    # y_range = [7e-4, 30]
    # yscale = "log"
    # temp_range = [-5, 480]
    # xscale = "linear"

    plot_type = "T2_frac"
    y_range = [0, 1]
    yscale = "linear"
    temp_range = [-5, 310]
    xscale = "linear"

    rates_to_plot = ["hopper"]

    ### Supp

    # plot_type = "T2_max_supp"
    # temp_range = [[-5, 480], [-5, 310]]
    # xscale = ["linear", "linear"]
    # yscale = ["log", "linear"]
    # y_range = [[7e-4, 30], [-0.02, 0.62]]
    # rates_to_plot = [["hopper", "wu"], ["hopper"]]

    ###

    main(
        file_name,
        path,
        plot_type,
        rates_to_plot,
        temp_range,
        y_range,
        xscale,
        yscale,
        dosave=False,
    )

    # file_name = "bar_gill_2012-2a.csv"
    # home = common.get_nvdata_dir()
    # path = home / "paper_materials/relaxation_temp_dependence/ripped_T2_plots"
    # bar_gill_replot(file_name, path)

    plt.show(block=True)
