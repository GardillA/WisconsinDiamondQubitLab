# -*- coding: utf-8 -*-
"""
Fit circles to superresolution rings in images demonstrating resolved
images of two NVs separated by less than the diffraction limit.
By full auto I mean you don't have to tell the program how many circles
there are - it'll figure it out for you by finding the best circle,
finding the next best circle, and so on until there are no more good
circles left.

Created on February 25, 2022

@author: mccambria
"""

# region Imports

import utils.tool_belt as tool_belt
import utils.common as common
import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, brute
from numpy import pi
from matplotlib.patches import Circle
import cv2 as cv
import sys
import utils.kplotlib as kpl
import cProfile
import functools

# import multiprocessing
# import pathos
from functools import partial
import time
import superresolution_circles_fake_data as fake_data
from pathos.multiprocessing import ProcessingPool as Pool

# endregion

# region Constants

num_circle_samples = 1000

phi_linspace = np.linspace(0, 2 * pi, num_circle_samples, endpoint=False)
cos_phi_linspace = np.cos(phi_linspace)
sin_phi_linspace = np.sin(phi_linspace)

# endregion

# region Functions


def calc_errors(image_file_name, circle_a, circle_b):

    cost_func = cost1

    # Get the image as a 2D ndarray
    image_file_dict = tool_belt.get_raw_data(image_file_name)
    image = np.array(image_file_dict["readout_image_array"])
    image_domain = image.shape
    image_len_x = image_domain[1]
    image_len_y = image_domain[0]

    ret_vals = process_image(image)
    sigmoid_image = ret_vals[-1]

    fig, axes_pack = plt.subplots(1, 3)
    fig.set_tight_layout(True)

    num_points = 1000
    sweep_half_range = 10
    # for circle in [circle_a]:
    for circle in [circle_a, circle_b]:

        print(circle)
        args = [
            sigmoid_image,
            image_len_x,
            image_len_y,
            False,
            False,
            [],
        ]
        opti_cost = cost_func(circle, *args)

        for param_ind in range(3):

            ax = axes_pack[param_ind]
            sweep_center = circle[param_ind]
            sweep_vals = np.linspace(
                sweep_center - sweep_half_range,
                sweep_center + sweep_half_range,
                num_points,
            )

            cost_vals = []
            for sweep_ind in range(num_points):
                test_circle = list(circle)
                test_circle[param_ind] = sweep_vals[sweep_ind]
                cost_vals.append(0.5 - cost_func(test_circle, *args))

            ax.plot(sweep_vals, cost_vals)

            left_width = None
            right_width = None
            if param_ind == 2:
                target_max_ratio = 0.5
            else:
                target_max_ratio = 1 - (1 / pi)
            half_max = target_max_ratio * (0.5 - opti_cost)
            half_ind = num_points // 2
            for delta in range(half_ind):
                test_ind = half_ind - delta
                if (cost_vals[test_ind] < half_max) and (left_width is None):
                    left_width = sweep_vals[test_ind]
                test_ind = half_ind + delta
                if (cost_vals[test_ind] < half_max) and (right_width is None):
                    right_width = sweep_vals[test_ind]
                if (left_width is not None) and (right_width is not None):
                    break
            if right_width is not None and left_width is not None:
                half_width_at_half_max = (right_width - left_width) / 2
                print(half_width_at_half_max)

    print()


def norm_0_to_1(image):
    normed_image = (image - np.min(image)) / (np.max(image) - np.min(image))

    return normed_image


def cost0(params, image, x_lim, y_lim, debug, invert, excluded_centers):
    """
    Faux-integrate the pixel values around the circle. By faux-integrate I mean
    average the values under a 1000 point, linearly spaced sampling of the circle.

    excluded_centers: list of tuples (x,y) describing potential centers that we don't
    want to consider. Useful for ignoring circles we've already found when optimizing.
    """

    circle_center_x, circle_center_y, circle_radius = params

    # Check if the circle is in the exclusion range
    exclusion_range = 4  # Pixels
    for excluded_x, excluded_y in excluded_centers:
        dist = np.sqrt(
            (excluded_x - circle_center_x) ** 2
            + (excluded_y - circle_center_y) ** 2
        )
        if dist < exclusion_range:
            return 1

    circle_samples_x = np.around(
        circle_center_x + circle_radius * cos_phi_linspace
    )
    circle_samples_x = circle_samples_x.astype(int)
    circle_samples_y = np.around(
        circle_center_y + circle_radius * sin_phi_linspace
    )
    circle_samples_y = circle_samples_y.astype(int)
    circle_samples = zip(circle_samples_x, circle_samples_y)

    check_valid = lambda el: (0 <= el[1] < x_lim) and (0 <= el[0] < y_lim)
    num_valid_samples = 0
    integrand = 0
    for coords in circle_samples:
        if check_valid(coords):
            num_valid_samples += 1
            integrand += image[coords]

    cost = integrand / num_valid_samples
    # cost = np.sum(integrand) / len(integrand)
    if invert:  # Best should be minimum
        cost = 1 - cost

    return cost


def cost1(params, image, x_lim, y_lim, debug, invert, excluded_centers):
    """
    Same as cost0, but better numpy use
    """

    x, y, r = params

    # Check if the circle is in the exclusion range
    exclusion_range = 4  # Pixels
    for excluded_x, excluded_y in excluded_centers:
        dist = np.sqrt((excluded_x - x) ** 2 + (excluded_y - y) ** 2)
        if dist < exclusion_range:
            return 1

    circle_samples_x = np.around(x + r * cos_phi_linspace).astype(int)
    circle_samples_y = np.around(y + r * sin_phi_linspace).astype(int)
    circle_samples = zip(circle_samples_x, circle_samples_y)

    check_valid = lambda el: (0 <= el[1] < x_lim) and (0 <= el[0] < y_lim)
    valid_samples = [check_valid(el) for el in circle_samples]

    integrand = np.sum(
        image[circle_samples_x[valid_samples], circle_samples_y[valid_samples]]
    )

    cost = integrand / np.count_nonzero(valid_samples)
    if invert:  # Best should be minimum
        cost = 1 - cost

    return cost


def sigmoid_quotient(laplacian, gradient):

    # Get the zeros from the gradient so we can avoid divide by zeros.
    # At the end we'll just set the sigmoid to the sign of the Laplacian
    # for these values.
    gradient_zeros = gradient < 1e-10
    gradient_not_zeros = np.logical_not(gradient_zeros)
    masked_gradient = (gradient * gradient_not_zeros) + gradient_zeros
    quotient = laplacian / masked_gradient
    sigmoid = 1 / (1 + np.exp(-1 * quotient))
    # sigmoid = 1 / (1 + np.exp(-5 * quotient - 0.0))
    # sigmoid = quotient
    laplacian_positive = np.sign(laplacian) == 1
    sigmoid = (sigmoid * gradient_not_zeros) + (
        laplacian_positive * gradient_zeros
    )
    return sigmoid


def sigmoid_quotient2(laplacian, gradient, orig):

    # Get the zeros from the gradient so we can avoid divide by zeros.
    # At the end we'll just set the sigmoid to the sign of the Laplacian
    # for these values.
    gradient_zeros = gradient < 1e-10
    gradient_not_zeros = np.logical_not(gradient_zeros)
    masked_gradient = (gradient * gradient_not_zeros) + gradient_zeros
    # quotient = laplacian / masked_gradient
    # quotient = 0.5 * laplacian / masked_gradient
    quotient = orig * laplacian / masked_gradient
    sigmoid = 1 / (1 + np.exp(-1 * quotient))
    laplacian_positive = np.sign(laplacian) == 1
    sigmoid = (sigmoid * gradient_not_zeros) + (
        laplacian_positive * gradient_zeros
    )
    return sigmoid


def process_image(image):

    # Blur
    gaussian_size = 7
    blur_image = cv.GaussianBlur(image, (gaussian_size, gaussian_size), 0)
    normed_blur_image = norm_0_to_1(blur_image)
    # zeroed_blur_image = blur_image - np.min(blur_image)
    # norm = np.percentile(
    #     zeroed_blur_image, 75
    # )  # Normalize to the ring brightness, roughly
    # # print(norm)
    # normed_blur_image = zeroed_blur_image / norm
    # # normed_blur_image = blur_image

    processing_root = normed_blur_image

    laplacian_image = cv.Laplacian(
        processing_root, cv.CV_64F, ksize=gaussian_size
    )
    # laplacian_image = norm_0_to_1(laplacian_image)
    # laplacian_image = 2*(laplacian_image - 0.5)
    # laplacian_image += 3
    # offset = 10 * np.average(np.abs(laplacian_image))
    # print(offset)
    # offset = np.sqrt(np.average(laplacian_image ** 2))
    # print(offset)
    # laplacian_image += offset
    # laplacian_image -= np.min(laplacian_image)

    sobel_x = cv.Sobel(processing_root, cv.CV_64F, 1, 0, ksize=gaussian_size)
    sobel_y = cv.Sobel(processing_root, cv.CV_64F, 0, 1, ksize=gaussian_size)
    gradient_image = np.sqrt(sobel_x ** 2 + sobel_y ** 2)
    # gradient_image = norm_0_to_1(gradient_image)
    # gradient_image += 100
    # gradient_image = gradient_image**2

    # sigmoid_image = sigmoid_quotient(laplacian_image, gradient_image)
    sigmoid_image = sigmoid_quotient2(
        laplacian_image, gradient_image, normed_blur_image
    )
    # sigmoid_image = cv.GaussianBlur(
    #     sigmoid_image, (gaussian_size, gaussian_size), 0
    # )

    return normed_blur_image, laplacian_image, gradient_image, sigmoid_image


def calc_distance(fig, x0, x1, y0, y1, sx0, sx1, sy0, sy1):

    dx = x1 - x0
    dy = y1 - y0
    distance = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    err = np.sqrt(
        ((dx / distance) ** 2) * (sx0 ** 2 + sx1 ** 2)
        + ((dy / distance) ** 2) * (sy0 ** 2 + sy1 ** 2)
    )

    print(distance)
    print(err)

    # 0.0004375 V, for Fig 4 each pixel is 0.0005 V. And the conversion is 34.8 um/V
    if fig == 3:
        conversion = 15.225  # nm / pixel
    elif fig == 4:
        conversion = 17.4  # nm / pixel
    distance_nm = conversion * distance
    err_nm = conversion * err
    print(distance_nm)
    print(err_nm)

    print()


# endregion


def main(
    image_file_name=None,
    image=None,
    passed_circles=None,
    brute_bounds=None,
    # "manual", "publication", "recursive", "fixed_r_reconstruction", "full_auto"
    run_type="full_auto",
):

    # region Setup

    main_text_fig_file_names = [
        "2021_09_30-13_18_47-johnson-dnv7_2021_09_23",
        "2021_10_17-19_02_22-johnson-dnv5_2021_09_23",
    ]

    cost_func = cost1

    if image_file_name == "fake":
        image = fake_data.main()
    elif image_file_name is not None:
        # Get the image as a 2D ndarray
        image_file_dict = tool_belt.get_raw_data(image_file_name)
        image = np.array(image_file_dict["readout_image_array"])

    image_domain = image.shape
    image_len_x = image_domain[1]
    image_len_y = image_domain[0]

    # image *= 10
    ret_vals = process_image(image)
    blur_image, laplacian_image, gradient_image, sigmoid_image = ret_vals
    normed_blur_image = norm_0_to_1(blur_image)
    inv_sigmoid_image = 1 - sigmoid_image
    final_image = normed_blur_image * inv_sigmoid_image

    # mean_blur_image = np.percentile(blur_image, 1)
    # normed_blur_image = (blur_image - mean_blur_image) / mean_blur_image
    # sigmoid_blur_image = 1 / (1 + np.exp(-1 * normed_blur_image))
    # sigmoid_blur_image = norm_0_to_1(sigmoid_blur_image)
    # # sigmoid_blur_image = (sigmoid_blur_image - 0.5) * 2
    # final_image = sigmoid_blur_image * inv_sigmoid_image
    # # final_image = (blur_image > mean_blur_image) * inv_sigmoid_image
    final_image = 1 - final_image

    opti_image = sigmoid_image
    # opti_image = final_image
    plot_image = opti_image

    # opti_image = normed_blur_image
    # plot_image = final_image
    # plot_image = gradient_image
    # plot_image = normed_blur_image
    # plot_image = sigmoid_blur_image
    # plot_image = image

    invert_cost = False
    # invert_cost = True

    # Plot the image
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    img = ax.imshow(plot_image, cmap="inferno")
    _ = plt.colorbar(img)
    # return

    # fig2, ax2 = plt.subplots()
    # fig2.set_tight_layout(True)
    # img2 = ax2.imshow(sigmoid_image, cmap="inferno")
    # _ = plt.colorbar(img2)

    # fig3, ax3 = plt.subplots()
    # fig3.set_tight_layout(True)
    # img3 = ax3.imshow(gradient_image, cmap="inferno")
    # _ = plt.colorbar(img3)

    # fig4, ax4 = plt.subplots()
    # fig4.set_tight_layout(True)
    # img4 = ax4.imshow(laplacian_image, cmap="inferno")
    # _ = plt.colorbar(img4)

    # fig5, ax5 = plt.subplots()
    # fig5.set_tight_layout(True)
    # img5 = ax5.imshow(blur_image, cmap="inferno")
    # _ = plt.colorbar(img5)

    # return

    # endregion

    # region Circle finding

    excluded_centers = []
    args = [
        opti_image,
        image_len_x,
        image_len_y,
        False,
        invert_cost,
        excluded_centers,
    ]
    plot_circles = []

    start = time.time()

    if run_type == "manual":

        plot_circles = passed_circles

    elif run_type == "reconstruction":

        # Determine the best radius by a full auto run desinged to find
        # the single best circle in the image. Then use that radius to
        # plot the full image.
        ret_vals = main(image_file_name, run_type="full_auto")
        fit_circles = ret_vals[0]
        radii = [el[2] for el in fit_circles]

        # Fig. 3
        if image_file_name == "2021_09_30-13_18_47-johnson-dnv7_2021_09_23":
            half_range = 18
            radius = np.mean(radii)
        # Fig. 4
        elif image_file_name == "2021_10_17-19_02_22-johnson-dnv5_2021_09_23":
            half_range = 20
            radius = np.mean(radii)
            # radius = 26.845
        else:
            half_range = round(0.15 * image_len_x)
            radius = radii[0]
            # radius = 45.98

        if image_file_name in main_text_fig_file_names:
            radius_half_range = 0.05 * radius
            num_points_r = 20
            rad_linspace = np.linspace(
                radius - radius_half_range,
                radius + radius_half_range,
                num_points_r,
            )
        else:
            num_points_r = 1
            rad_linspace = [radius]

        num_points = 1000
        half_len_x = image_len_x // 2
        half_len_y = image_len_x // 2
        x_linspace = np.linspace(
            half_len_x - half_range, half_len_x + half_range, num_points
        )
        y_linspace = np.linspace(
            half_len_y - half_range, half_len_y + half_range, num_points
        )

        # reconstruction = []
        # for y in y_linspace:
        #     for x in x_linspace:
        #         cost_lambda = lambda radius: cost_func([y, x, radius], *args)
        #         with Pool() as p:
        #             rad_line = p.map(cost_lambda, radius_linspace)
        #         reconstruction.append(x_line)

        # Make one big list of all the coordinates we need to try so that we can
        # efficiently run through the list on multiple threads
        coords_list = []
        for y in y_linspace:
            for x in x_linspace:
                for r in rad_linspace:
                    coords_list.append([y, x, r])
        cost_lambda = lambda coords: cost_func(coords, *args)

        with Pool() as p:
            reconstruction_list = p.map(cost_lambda, coords_list)

        # Loop through the big list and build the reconstruction out of it,
        # picking the best cost from the available r values at each point
        ind = 0
        reconstruction = []
        for y in y_linspace:
            reconstruction.append([])
            for x in x_linspace:
                r_costs = reconstruction_list[ind : ind + num_points_r]
                opti_val = min(r_costs)
                reconstruction[-1].append(opti_val)
                ind += num_points_r

        fig2, ax = plt.subplots()
        fig2.set_tight_layout(True)
        extent = [
            min(x_linspace),
            max(x_linspace),
            max(y_linspace),
            min(y_linspace),
        ]
        img = ax.imshow(
            reconstruction, cmap="inferno_r", extent=extent  # , vmin=0.35
        )
        cbar = plt.colorbar(img)

        if "circle_centers" in image_file_dict:
            circle_centers = image_file_dict["circle_centers"]
            color = kpl.KplColors.GRAY.value
            for circle in circle_centers:
                circle_patch = Circle(
                    (circle[1], circle[0]),
                    0.25,
                    fill=color,
                    color=color,
                )
                ax.add_patch(circle_patch)

        end = time.time()
        print(f"Elapsed time: {end-start}")

        if image_file_name is not None:
            image_file_name_split = image_file_name.split("-")
            timestamp = "-".join(image_file_name_split[0:2])
            nv_name = "-".join(image_file_name_split[2:4])
            file_path = tool_belt.get_file_path(__file__, timestamp, nv_name)
            tool_belt.save_figure(fig2, file_path)

        return

    elif run_type == "full_auto":

        while True:

            if image_file_name in main_text_fig_file_names:
                rad_bounds = (26, 30)
                num_circles = 2
            else:
                rad_bounds = (38, 50)
                num_circles = 1

            # Define the bounds of the optimization
            if brute_bounds is None:
                bounds = [
                    (0.4 * image_len_y, 0.6 * image_len_y),
                    (0.4 * image_len_x, 0.6 * image_len_x),
                    rad_bounds,
                ]
            else:
                bounds = brute_bounds

            # Anything worse than minimum_cost and we stop searching for circles
            minimum_cost = None
            # num_circles = 1

            # Initialize best_cost
            best_cost = 1

            while True:

                # Update args with the last found opti_circle
                args[-1] = excluded_centers

                opti_circle = brute(
                    cost_func,
                    bounds,
                    Ns=30,
                    args=args,
                    finish=None,
                    workers=-1,  # Multiprocessing: -1 means use as many cores as available
                )
                args[-1] = []
                new_best_cost = cost_func(opti_circle, *args)

                threshold = 0.0001 * best_cost
                # threshold = 0.001 * best_cost
                if best_cost - new_best_cost < threshold:
                    break
                best_cost = new_best_cost

                bounds_span = [el[1] - el[0] for el in bounds]
                half_new_span = [0.1 * el for el in bounds_span]
                bounds = [
                    (
                        opti_circle[ind] - half_new_span[ind],
                        opti_circle[ind] + half_new_span[ind],
                    )
                    for ind in range(3)
                ]
                # print(new_best_cost)
                # print(bounds)

            # Quit if we've found the number of specified circles. ELse quit if
            # our best is worse than the threshold we set and we have at least
            # one circle
            if num_circles is None:
                if (len(plot_circles) > 0) and (new_best_cost > minimum_cost):
                    break
            excluded_centers.append(opti_circle[0:2])
            plot_circles.append(opti_circle)
            if (num_circles is not None) and (
                len(plot_circles) == num_circles
            ):
                break

    end = time.time()
    print(f"Elapsed time: {end-start}")

    # endregion

    # region Circle plotting

    for circle in plot_circles:

        # Debug tweak
        # circle[0] -= 0.5

        # Report what we found
        rounded_circle = [round(el, 2) for el in circle]
        rounded_cost = round(cost_func(circle, *args), 5)
        # print("{} & {} & {} & {}".format(*rounded_circle, rounded_cost))
        print("{}, {}, {}, {}".format(*rounded_circle, rounded_cost))

        # Plot the circle
        circle_patch = Circle(
            (circle[1], circle[0]),
            circle[2],
            fill=False,
            color="w",
        )
        ax.add_patch(circle_patch)

        # Plot the center
        circle_patch = Circle(
            (circle[1], circle[0]),
            0.5,
            fill="w",
            color="w",
        )
        ax.add_patch(circle_patch)

    # endregion
    return plot_circles, fig


# region Run the file

if __name__ == "__main__":

    # # Fig 3
    # calc_distance(3, 36.85, 43.87, 41.74, 39.05, 1.4, 0.9, 1.5, 1.3)
    # # Fig 4
    # calc_distance(4, 45.79, 56.32, 50.98, 51.2, 1.0, 1.1, 1.4, 1.4)

    # sys.exit()

    # tool_belt.init_matplotlib()

    # circles = [3]
    # circles = [4]
    circles = [3, 4]
    for circle in circles:

        # Fig. 3
        if circle == 3:
            image_file_name = "2021_09_30-13_18_47-johnson-dnv7_2021_09_23"
            # Best circles by hand
            # circle_a = [41.5, 37, 27.5]
            # circle_b = [40, 44, 27.75]
            # Recursive brute results, 1000 point circle
            circle_a = [41.73, 36.83, 27.72]  # 0.31941
            # errs_a = [1.4, 1.4, 1.2]
            circle_b = [39.1, 43.9, 27.64]  # 0.36108
            # errs_b = [1.5, 0.9, 1.0]
            circle_c = [39.97, 41.93, 29.46]  # test
            # circle_c = [20.1, 41.93, 29.46]  # test

        # Fig. 4
        elif circle == 4:
            image_file_name = "2021_10_17-19_02_22-johnson-dnv5_2021_09_23"
            # Best circles by hand
            # circle_a = [50, 46, 26]
            # circle_b = [51.7, 56.5, 27.3]
            # Recursive brute results, 1000 point circle
            circle_a = [50.64, 45.64, 26.29]  # 0.3176
            # circle_a = [50.98 - 0, 45.79 + 0, 26.14 - 0]  # 0.3176
            # errs_a = [2.1, 1.1, 1.3]
            circle_b = [51.02, 56.15, 27.48]  # 0.35952
            # errs_b = [1.8, 1.1, 1.2]
            circle_c = None

        # circles = [circle_a, circle_b, circle_c]
        circles = [circle_a, circle_b]
        # main(image_file_name, run_type="full_auto")
        # main(image_file_name, run_type="reconstruction")
        # main(image_file_name, passed_circles=circles, run_type="manual")
        # main(image_file_name, circle_a, circle_b, run_type="full_auto")
        calc_errors(image_file_name, circle_a, circle_b)

    image_file_name = "2022_07_17-20_55_18-johnson-nv0_2021_12_22-faked"
    # image_file_name = "2022_07_19-10_59_00-johnson-nv0_2021_12_22-faked"
    # main(image_file_name, run_type="reconstruction")
    # cProfile.run('main(image_file_name, run_type="reconstruction")')

    # image_file_name = "2021_09_30-13_18_47-johnson-dnv7_2021_09_23"
    # main(image_file_name=image_file_name, run_type="fixed_r_reconstruction")

    plt.show(block=True)

# endregion


#  0.0004375 V, for Fig 4 each pixel is 0.0005 V. And the conversion is 34.8 um/V

#  Fig 3: 15.225 nm / pixel
#  Fig 4: 17.4 nm / pixel


51.42, 56.1, 27.56