# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 22:33:25 2018

@author: Matt
"""

import matplotlib.pyplot as plt
import numpy
import matplotlib.animation as animation
import time

# %% Data

xDim = 4
yDim = 3

resolution = .1

offset = [0.0,  0.0]

vals = numpy.array([1, 2, 10, 5,
                    0, 1, 1, 9,
                    10, 10, 5, 4])

# %% Figure and image setup

imgVals = numpy.zeros((yDim, xDim))

left = offset[0]
right = left + (xDim * resolution)
top = offset[1]
bottom = top + (yDim * resolution)

imageExtent = [left, right, bottom, top]
centImageExtent = [x - (resolution / 2) for x in imageExtent]

fig, ax = plt.subplots()
img = ax.imshow(imgVals, cmap="gray", extent=tuple(centImageExtent))
img.autoscale()
fig.canvas.draw()
fig.canvas.flush_events()
background = fig.canvas.copy_from_bbox(ax.bbox)

# %% Update data
tstart = time.time()
for index in range(yDim):
    start = index * xDim
    if index % 2 == 0:
        extension = vals[start: start+xDim]
        imgVals[index, 0: xDim] = extension
    else:
        extension = vals[start: start+xDim]
        extension = extension[::-1]  # Reverse the list
        imgVals[index, 0: xDim] = extension

    fig.canvas.restore_region(background)

    img.set_data(imgVals)
    img.autoscale()
    ax.draw_artist(img)
    # Also draw the axes so the updated image doesn't overwrite them
    for spine in ax.spines.values():
        ax.draw_artist(spine)

    fig.canvas.blit(ax.bbox)
    fig.canvas.flush_events()

    #time.sleep(1)

tend = time.time()
print(tend-tstart)
