# -*- coding: utf-8 -*-
"""
Trial ellipse fitting

Created October 3rd, 2022

@author: mccambria
"""

### Imports
import numpy as np
from numpy.core.shape_base import block
from scipy.optimize import root_scalar, minimize_scalar, minimize
import utils.tool_belt as tool_belt
import time
import matplotlib.pyplot as plt
from utils import common
from utils import kplotlib as kpl
from utils.kplotlib import KplColors
from scipy.optimize import curve_fit
import csv
import sys
import random
from pathos.multiprocessing import ProcessingPool as Pool

cent = 0.5
amp = 0.65 / 2

# region Functions


def ellipse_point(theta, phi):
    return (cent + amp * np.cos(theta + phi), cent + amp * np.cos(theta - phi))


def theta_cost(phi, ellipse_point, test_point):

    theta_cost = lambda theta: sum(
        [(el[0] - el[1]) ** 2 for el in zip(ellipse_point(theta), test_point)]
    )

    # Guess theta by assuming theta and the ellipse amplitude are free and
    # solving for the position of the point on the ellipse
    # x = cent + amp * np.cos(theta + phi)
    # y = cent + amp * np.cos(theta - phi)
    # x-y = amp * (np.cos(theta + phi) - np.cos(theta - phi))
    #     = -2 amp sin(theta) sin(phi)
    # amp = (x-y) / (-2 sin(theta) sin(phi))
    # x+y = 2cent + 2 amp cos(theta) cos(phi)
    #     = 2cent - (x-y) cot(theta) cot(phi)
    x, y = test_point
    # Avoid divide by zero
    guess_denom = np.tan(phi) * (2 * cent - (x + y))
    if guess_denom < 0.01:
        base_guess = np.pi / 2
        # Check the guess and its compliment
        guesses = [
            base_guess % (2 * np.pi),
            (2 * np.pi - base_guess) % (2 * np.pi),
        ]
    else:
        guess_arg = (x - y) / guess_denom
        base_guess = np.arctan(abs(guess_arg)) % np.pi
        # Check the guess and its compliments
        guesses = [
            base_guess % (2 * np.pi),
            (np.pi - base_guess) % (2 * np.pi),
            (base_guess + np.pi) % (2 * np.pi),
            (2 * np.pi - base_guess) % (2 * np.pi),
        ]

    best_cost = None
    for guess in guesses:
        res = minimize(theta_cost, guess)
        opti_theta = res.x
        opti_cost = res.fun
        if (best_cost is None) or (opti_cost < best_cost):
            best_cost = opti_cost
    return best_cost


def ellipse_cost(phi, points, debug=False):

    if debug:
        test = 1

    ellipse_lambda = lambda theta: ellipse_point(theta, phi)

    # The cost is the rms distance between the point and the ellipse
    # Finding the closest distance between an arbitary point and an ellipse,
    # of course, turns out to be a hard problem, so let's just run another
    # minimization for it
    theta_cost_lambda = lambda point: theta_cost(phi, ellipse_lambda, point)
    with Pool() as p:
        theta_costs = p.map(theta_cost_lambda, points)

    cost = sum(theta_costs)
    num_points = len(points)
    cost = np.sqrt(cost / num_points)

    return cost


def gen_ellipses():
    # phis = [0, np.pi / 2, np.pi / 4]
    phis = [0.0]
    num_points = 100
    ellipses = []
    noise_amp = 0.05
    for phi in phis:
        ellipse = [[phi]]
        ellipse_lambda = lambda theta: ellipse_point(theta, phi)
        points = []
        for ind in range(num_points):
            theta = 2 * np.pi * random.random()
            point = ellipse_lambda(theta)
            noisy_point = (
                point[0] + noise_amp * random.random(),
                point[1] + noise_amp * random.random(),
            )
            points.append(noisy_point)
            ellipse.extend(points)
        ellipses.append(ellipse)
    return ellipses


def populate_imported_phis(path, files):

    phis_by_file = []

    for el in files:
        sub_phis = []
        with open(path / el, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                float_row = [float(el) for el in row]
                sub_phis.append(float_row[0])
        phis_by_file.append(sub_phis)

    # Re-sort by ellipse
    num_files = len(files)
    num_ellipses = len(phis_by_file[0])
    phis = []
    for ellipse_ind in range(num_ellipses):
        sub_phis = []
        for sub_ind in range(num_files):
            sub_phis.append(phis_by_file[sub_ind][ellipse_ind])
        phis.append(sub_phis)

    return phis


def import_ellipses(path):

    x_vals = []
    y_vals = []
    phis = []
    file_x = "testingX.csv"
    file_y = "testingY.csv"
    phi_files = ["testingPhi_d.csv", "phi_LS.csv", "phi_NN.csv"]

    trans_x = zip(*csv.reader(open(path / file_x, newline="")))
    for row in trans_x:
        float_row = [float(el) for el in row]
        x_vals.append(float_row)

    trans_y = zip(*csv.reader(open(path / file_y, newline="")))
    for row in trans_y:
        float_row = [float(el) for el in row]
        y_vals.append(float_row)

    phis = populate_imported_phis(path, phi_files)

    ellipses = []
    for ellipse_x, ellipse_y, ellipse_phi in zip(x_vals, y_vals, phis):
        points = list(zip(ellipse_x, ellipse_y))
        ellipse = [ellipse_phi]
        ellipse.extend(points)
        ellipses.append(ellipse)

    return ellipses


# endregion

# region Main


def main(path):

    # ellipses = import_ellipses(path)
    ellipses = gen_ellipses()
    theta_linspace = np.linspace(0, 2 * np.pi, 100)
    phi_errors = []

    do_plot = True

    # ellipses = ellipses[::10]
    # ellipses = [ellipses[35]]  # 34, 38
    for ind in range(len(ellipses)):

        ellipse = ellipses[ind]
        phi_errors_sub = []
        # for ellipse in ellipses[0]:
        ellipse_phis = ellipse[0]
        true_phi = ellipse_phis.pop(0)
        algo_phis = ellipse_phis
        points = ellipse[1:]

        thresh = 0.1
        if thresh < true_phi < (np.pi / 2) - thresh:
            continue
        print(f"Ellipse index: {ind}")

        # Avoid local minima by using multiple guesses
        # guesses = [0, np.pi / 2, np.pi / 4]
        # for guess in guesses:
        #     res = minimize(ellipse_cost, guess, args=(points,))
        #     opti_phi = res.x[0]
        #     # Remove degeneracies
        #     opti_phi = opti_phi % np.pi
        #     if opti_phi > np.pi / 2:
        #         opti_phi = np.pi - opti_phi
        #     cost = ellipse_cost(opti_phi, points, True)
        #     if cost < 0.1:
        #         break
        opti_phi = 0
        algo_phis.insert(0, opti_phi)

        if do_plot:
            fig, ax = plt.subplots()
            for test_phi in [0.0, 0.03]:
                # Plot the data points
                for point in points:
                    color = KplColors.BLUE.value
                    kpl.plot_data(ax, *point, color=color)
                # Plot the fit
                ellipse_lambda = lambda theta: ellipse_point(theta, test_phi)
                x_vals, y_vals = zip(ellipse_lambda(theta_linspace))
                x_vals = x_vals[0]
                y_vals = y_vals[0]
                kpl.plot_line(ax, x_vals, y_vals)
                kpl.tight_layout(fig)
                cost = ellipse_cost(test_phi, points, True)
                print(f"Phi: {round(test_phi, 3)}; cost: {round(cost, 6)}")

        # Get the costs
        # test_phis = [true_phi]
        # test_phis.extend(algo_phis)
        # for phi in test_phis:
        #     cost = ellipse_cost(phi, points, True)
        #     print(f"{round(phi, 3)}: {round(cost, 6)}")
        # print()
        phi = opti_phi
        cost = ellipse_cost(phi, points, True)
        if cost > 0.1:
            print("Algorithm did poorly...")
            print(f"Phi: {round(phi, 3)}; cost: {round(cost, 6)}")

        # Get the phi errors
        for phi in algo_phis:
            phi_errors_sub.append(phi - true_phi)
        phi_errors.append(phi_errors_sub)

    print(f"Summary for ellipse with true phi {true_phi}")
    # for row in phi_errors:
    #     rounded_errs = [round(el, 6) for el in row]
    #     print(
    #         f"Algo: {rounded_errs[0]}; LS: {rounded_errs[1]}; NN:"
    #         f" {rounded_errs[2]}"
    #     )
    # print()
    # print("RMS phase errors for algorithm, least squares, neural net: ")
    print("RMS phase errors for algorithm")
    phi_errors = np.array(phi_errors)
    rms_phi_errors = np.sqrt(np.mean(phi_errors ** 2, axis=0))
    print([round(el, 6) for el in rms_phi_errors])


# endregion

if __name__ == "__main__":

    kpl.init_kplotlib()

    home = common.get_nvdata_dir()
    path = home / "ellipse_data"

    main(path)

    plt.show(block=True)

    sys.exit()

    ellipses = import_ellipses(path)

    true_phis = [el[0][0] for el in ellipses]

    phi_errs = [
        [-0.002696, -0.006651, -0.014066],
        [-0.005236, -0.025661, 0.013953],
        [-0.002495, -0.003856, -0.172563],
        [-0.010397, 0.01137, 0.013603],
        [-0.005612, 0.011135, 0.01627],
        [0.004582, 0.02866, -0.023009],
        [-0.001166, -0.010305, -0.085988],
        [-0.003142, 0.014189, -0.039819],
        [-0.00672, -0.041953, -0.007913],
        [0.019494, 0.028655, 0.018008],
        [0.006595, 0.006517, -0.004063],
        [0.00969, 0.020179, -0.00584],
        [-0.001045, -0.006231, -0.188581],
        [-0.005577, -0.010907, -0.15544],
        [0.002375, 0.010187, 0.03117],
        [-0.020186, -0.014631, 0.014923],
        [-0.016039, -0.014665, -0.050193],
        [-0.007915, -0.007319, -0.022615],
        [0.011466, 0.027088, -0.002615],
        [-0.007573, -0.006274, -0.021856],
        [-0.012851, -0.027916, -0.001479],
        [-0.00625, 0.017605, -0.020141],
        [0.010959, -0.007834, -0.031364],
        [0.007106, 0.014636, -0.060448],
        [-0.005491, -0.002544, -0.039629],
        [0.013471, 0.019428, -0.105871],
        [0.00267, -0.003474, 0.008865],
        [0.002729, 0.019037, 0.003506],
        [-0.003578, -0.006241, 0.081942],
        [-0.006198, -0.025533, 0.002916],
        [0.014313, -0.006932, 0.022175],
        [-0.007534, -0.007198, -0.034657],
        [0.003608, -0.00379, 0.026548],
        [0.002803, 0.006552, 0.004461],
        [-0.029468, -0.057791, -0.025558],
        [-0.002701, 0.034609, -0.004031],
        [0.004912, 0.007318, 0.024522],
        [-0.004386, 0.00847, 0.016392],
        [0.003009, -0.035114, -0.006209],
        [-0.002204, -0.001717, 0.034737],
        [0.001627, -0.022486, 0.022379],
        [-0.003329, -0.023689, 0.011969],
        [-0.006071, 0.011426, -0.021362],
        [0.001658, -0.002936, -0.026672],
        [-0.005253, 0.006608, 0.02229],
        [-0.004552, 0.017152, -0.016773],
        [-0.018404, -0.019594, -0.032711],
        [-0.000391, 0.013142, 0.008744],
        [0.009777, 0.012604, -0.224467],
        [-0.008527, -0.007891, -0.028047],
        [-0.000372, -0.00768, -0.032582],
        [-0.009247, -0.011221, 0.011764],
        [0.003607, -0.005237, 0.01255],
        [0.004043, 0.008858, -0.008026],
        [-0.002131, 0.015301, 0.024754],
        [0.004303, 0.00297, 0.039537],
        [-0.000208, -0.038235, 0.015747],
        [-0.006015, -0.02749, 0.005508],
        [0.014381, 0.008799, 0.063291],
        [-0.009208, 0.013698, -0.020051],
        [0.011965, -0.013355, 0.018704],
        [-0.002119, 0.014874, -0.006485],
        [0.000371, 0.010341, -0.12445],
        [0.013484, 0.019017, -0.020867],
        [-0.003873, -0.019166, -0.088132],
        [-0.006543, 0.001746, -0.142866],
        [0.001001, -0.007538, 0.064274],
        [0.004758, 0.017759, -0.040068],
        [0.008269, 0.008184, -0.115205],
        [-0.011544, -0.005907, 0.069389],
        [0.005727, 0.006894, -0.097154],
        [0.021964, 0.029796, -0.090545],
        [-0.008503, -0.043181, 0.006695],
        [-0.008904, 0.016486, -0.011963],
        [0.005893, -0.014571, 0.011854],
        [0.00064, 0.012811, -0.003998],
        [-0.005962, -0.025845, 0.018181],
        [-0.009019, -0.025957, -0.140684],
        [-0.004009, -0.008738, -0.00384],
        [0.006963, 0.038078, -0.001536],
        [0.000145, -0.003798, -0.036581],
        [-0.008003, -0.002053, -0.009636],
        [0.010323, 0.037769, -0.003184],
        [0.011814, 0.017541, -0.030662],
        [0.001115, -0.029644, 0.019089],
        [-0.010788, -0.01479, -0.049583],
        [0.003803, 0.035284, -0.010491],
        [-0.013374, -0.017504, -0.099624],
        [0.007496, 0.006683, -0.001683],
        [0.003659, 0.001623, -0.015948],
        [0.007829, -0.023922, -0.029186],
        [0.005227, 0.015435, -0.056644],
        [0.005001, 0.023689, 0.00774],
        [-0.00057, -0.003199, 0.027986],
        [0.024445, 0.048901, 0.015959],
        [0.001972, 0.002686, -0.106259],
        [0.010003, -0.009473, -0.024237],
        [0.011258, 0.012474, -0.018905],
        [0.000721, -0.018711, -0.053612],
        [0.007917, -0.007447, 0.035453],
    ]

    # num_fails = 0
    # for ind in range(len(phi_errs)):
    #     row = phi_errs[ind]
    #     abs_errs = [abs(el) for el in row]
    #     best_ind = abs_errs.index(min(abs_errs))
    #     if best_ind != 0:
    #         print(f"{ind}: {best_ind}")
    #         num_fails += 1
    # print(num_fails)

    algo_errs = [el[0] for el in phi_errs]
    ls_errs = [el[1] for el in phi_errs]
    nn_errs = [el[2] for el in phi_errs]
    algo_errs = np.array(algo_errs)
    ls_errs = np.array(ls_errs)
    nn_errs = np.array(nn_errs)

    fig, ax = plt.subplots()
    colors = kpl.data_color_cycler
    kpl.plot_data(ax, true_phis, algo_errs, label="Algo", color=colors[0])
    # kpl.plot_data(ax, true_phis, ls_errs, label="LS", color=colors[1])
    # kpl.plot_data(ax, true_phis, nn_errs, label="NN", color=colors[2])
    ax.set_xlabel("True phase")
    ax.set_ylabel("Error")
    ax.legend()
    kpl.tight_layout(fig)
    plt.show(block=True)

    inds = []
    for ind in range(len(true_phis)):
        val = true_phis[ind]
        thresh = 0.2
        if (np.pi / 4) - thresh < val < (np.pi / 4) + thresh:
            inds.append(ind)
    print(len(inds))
    algo_rms = np.sqrt(np.mean(algo_errs[inds] ** 2))
    ls_rms = np.sqrt(np.mean(ls_errs[inds] ** 2))
    nn_rms = np.sqrt(np.mean(nn_errs[inds] ** 2))
    print(f"Alorithm rms error: {round(algo_rms, 6)}")
    print(f"Least squares rms error: {round(ls_rms, 6)}")
    print(f"Neural net rms error: {round(nn_rms, 6)}")
