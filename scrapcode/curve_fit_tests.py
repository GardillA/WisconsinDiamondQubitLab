# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 17:53:37 2019

@author: mccambria
"""

import numpy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def line(t, a, b):
    return a + (b * t)

def line_d_a(t, a, b):
    return 1 + (0 * t)

def line_d_b(t, a, b):
    return t

def exp_dec(t, rate, amp):
    return amp * numpy.exp(-rate * t)

def exp_dec_d_rate(t, rate, amp):
    return (-t) * amp * numpy.exp(-rate * t)

def exp_dec_d_amp(t, rate, amp):
    return numpy.exp(-rate * t)

def exp_dec_offset(t, rate, amp, offset):
    return offset + (amp * numpy.exp(-rate * t))

t = numpy.linspace(0, 10, 6)
y_vals = [9, 6, 4, 3, 2, 1.5]
y_st_devs = numpy.array([2, 1, 1, 0.5, 1, 0.5])
# y_st_devs *= 10
y_vars = y_st_devs**2
y_vars_inv = 1/y_vars

# plt.plot(t, y_vals)
# plt.errorbar(t, y_vals, yerr=sigmas, linestyle='None')

fit_func = line
# fit_func = exp_dec

# popt, pcov = curve_fit(fit_func, t, y_vals)
popt, pcov = curve_fit(fit_func, t, y_vals, sigma=y_st_devs, absolute_sigma=True)
print(pcov)

# Get the chi squared
# residuals = y_vals - fit_func(t, *popt)
# # chi_squared = sum((residuals) ** 2)  # without vars
# chi_squared = sum((residuals / y_st_devs) ** 2)  # with vars
# reduced_chi_squared = chi_squared / (len(t) - 2)  # chi squared per degree of freedom

# Compute q_matrix=(J^T W^2 J)^-1
q_matrix = numpy.empty((2,2))
if fit_func is line:
    q_matrix[0,0] = sum(line_d_a(t, *popt) * y_vars_inv * line_d_a(t, *popt))
    q_matrix[0,1] = sum(line_d_a(t, *popt) * y_vars_inv * line_d_b(t, *popt))
    q_matrix[1,0] = sum(line_d_b(t, *popt) * y_vars_inv * line_d_a(t, *popt))
    q_matrix[1,1] = sum(line_d_b(t, *popt) * y_vars_inv * line_d_b(t, *popt))
elif fit_func is exp_dec:
    q_matrix[0,0] = sum(exp_dec_d_rate(t, *popt) * y_vars_inv * exp_dec_d_rate(t, *popt))
    q_matrix[0,1] = sum(exp_dec_d_rate(t, *popt) * y_vars_inv * exp_dec_d_amp(t, *popt))
    q_matrix[1,0] = sum(exp_dec_d_amp(t, *popt) * y_vars_inv * exp_dec_d_rate(t, *popt))
    q_matrix[1,1] = sum(exp_dec_d_amp(t, *popt) * y_vars_inv * exp_dec_d_amp(t, *popt))

# Manually calculate the covariances to check what scipy is really doing
# Where f is the fit_func and pi is the ith fit param, the covariances are:
pcov_manual = numpy.linalg.inv(q_matrix)
print(pcov_manual)
print(pcov / pcov_manual)

# smooth_t = numpy.linspace(0, 10, 1000)
# plt.plot(smooth_t, fit_func(smooth_t, *popt))
