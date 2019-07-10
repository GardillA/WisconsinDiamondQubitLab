# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:59:18 2019

Plotting the beam waist to a gaussian beam

@author: Aedan
"""

import numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def beam_waist_eq(z, w_0, z_offset):
    return w_0*numpy.sqrt(1 + ( ((z + z_offset) * 0.589) / (numpy.pi * w_0**2) )**2 )

z = (0, 76.2*10**3, 101.6*10**3, 127*10**3, 152.4*10**3, 177.8*10**3 ) #  um

w = (2900 / 2, 1440 / 2, 900 / 2, 430 / 2, 100 / 2, 390 / 2) # um

init_params = (40, 150000)
opti_params, cov_arr = curve_fit(beam_waist_eq, z, w, init_params)

z_linspace = numpy.linspace(0, 200000, 1000)

print('Beam waist at focus: '+ str(opti_params[0]) + ' um')


fig, ax = plt.subplots(1, 1, figsize=(10, 8))

ax.plot(z_linspace, beam_waist_eq(z_linspace, *opti_params), label = 'fit')
ax.plot(z, w, 'r.', label = 'data')
ax.set_xlabel('z (um)')
ax.set_ylabel('beam waist (um)')
ax.legend()
