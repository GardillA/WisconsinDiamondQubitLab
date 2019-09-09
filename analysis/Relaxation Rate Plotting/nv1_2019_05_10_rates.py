# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:36:23 2019

This file plots the relaxation rate data collected for the nv1_2019_05_10.

The data is input manually, and plotted on a loglog plot with error bars along
the y-axis. A 1/f**2 line is also fit to the gamma rates to show the behavior.

@author: Aedan
"""

'''
nv1_2019_05_10


'''
# %%
def fit_eq_1(f, amp):
    return amp*f**(-1)

def fit_eq_2(f, amp):
    return amp*f**(-2)

def fit_eq_alpha(f, amp, alpha, offset):
    return offset + amp*f**(-alpha)

# %%

import matplotlib.pyplot as plt
from scipy import asarray as ar, exp
from scipy.optimize import curve_fit
import numpy

# The data
nv1_splitting_list = [19.8, 27.8, 28, 30, 32.7, 51.8, 97.8, 116, 268, 563.6, 1016.8]
nv1_omega_avg_list = [1.3, 1.25, 1.7, 1.62, 1.48, 2.3, 1.8, 1.18, 1.07, 1.15, 0.525]
nv1_omega_error_list = [0.2, 0.1, 0.4, 0, 0.09, 0.4, 0.2, 0.13, 0.10, 0.12, 0.041]
nv1_gamma_avg_list = [136, 61, 68, 37, 50, 13.0, 3.5, 4.6, 1.92, 0.78, 0.712]
nv1_gamma_error_list = [10, 1, 7, 6, 3, 0.6, 0.2, 0.3, 0.13, 0.1, 0.122]

# Try to fit the gamma to a 1/f^2

fit_1_params, cov_arr = curve_fit(fit_eq_1, nv1_splitting_list, nv1_gamma_avg_list, 
                                p0 = 100, sigma = nv1_gamma_error_list,
                                absolute_sigma = True)

fit_2_params, cov_arr = curve_fit(fit_eq_2, nv1_splitting_list, nv1_gamma_avg_list, 
                                p0 = 100, sigma = nv1_gamma_error_list,
                                absolute_sigma = True)

fit_alpha_params, cov_arr = curve_fit(fit_eq_alpha, nv1_splitting_list, nv1_gamma_avg_list, 
                                p0 = (100, 1, 2), sigma = nv1_gamma_error_list,
                                absolute_sigma = True)

splitting_linspace = numpy.linspace(10, 2000,
                                    1000)

omega_constant_array = numpy.empty([1000]) 
omega_constant_array[:] = numpy.average(nv1_omega_avg_list)


fig, ax = plt.subplots(1, 1, figsize=(10, 8))

#ax.errorbar(splitting_list, omega_avg_list, yerr = omega_error_list)

axis_font = {'size':'14'}

orange = '#f7941d'
purple = '#87479b'

print(fit_alpha_params)

ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
ax.errorbar(nv1_splitting_list, nv1_gamma_avg_list, yerr = nv1_gamma_error_list, 
            label = r'$\gamma$',  fmt='o',markersize = 12, color = purple)
ax.errorbar(nv1_splitting_list, nv1_omega_avg_list, yerr = nv1_omega_error_list, 
            label = r'$\Omega$', fmt='^', markersize = 12, color=orange)

#ax.plot(splitting_linspace, fit_eq_2(splitting_linspace, *fit_2_params), 
#            label = r'$f^{-2}$', color ='teal')
#ax.plot(splitting_linspace, fit_eq_1(splitting_linspace, *fit_1_params), 
#            label = r'$f^{-1}$', color = 'orange')

# %%

ax.plot(splitting_linspace, fit_eq_alpha(splitting_linspace, *fit_alpha_params), 
             linestyle='dashed', linewidth=3, label = r'fit',color =purple)
ax.plot(splitting_linspace, omega_constant_array, color = orange,
            linestyle='dashed', linewidth=3, label = r'$\Omega$')

text = '\n'.join((r'$1/f^{\alpha} + \gamma_\infty$ fit:',
                  r'$\alpha = $' + '%.2f'%(fit_alpha_params[1]),
                  r'$\gamma_\infty = $' + '%.2f'%(fit_alpha_params[2]),
                  r'$A_0 = $' + '%.0f'%(fit_alpha_params[0])
#                  ,r'$a = $' + '%.2f'%(fit_params[2])
                  ))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.85, 0.7, text, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# %%

ax.tick_params(which = 'both', length=6, width=2, colors='k',
                grid_alpha=0.7, labelsize = 18)

ax.tick_params(which = 'major', length=12, width=2)

ax.grid()

ax.set_xlim([10,1200])
ax.set_ylim([0.1,300])

plt.xlabel('Splitting (MHz)', fontsize=18)
plt.ylabel('Relaxation Rate (kHz)', fontsize=18)
#plt.title('NV 1', fontsize=18)
#ax.legend(fontsize=18)

#fig.savefig("C:/Users/Aedan/Creative Cloud Files/Paper Illustrations/Magnetically Forbidden Rate/fig_3c.pdf", bbox_inches='tight')
