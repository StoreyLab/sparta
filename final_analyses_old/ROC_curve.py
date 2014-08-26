# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 13:06:31 2014

@author: peter
"""

# based on example from:
# http://www.walkingrandomly.com/?p=5215
#import os
#import math
import numpy as np
#from scipy.optimize import curve_fit
#import matplotlib as mpl
#mpl.use('GTK3Agg')
from matplotlib import pyplot as plt

filename = 'ROC_out'
title = 'Error Rate vs. Sensitivity: Simulated BY/RM Reads'

err = []
sens = []
err_naive = []
sens_naive = []

with open (filename, 'r') as inf:
    for l in inf:
        elements = l.rstrip().split('\t')
        err.append(float(elements[1]))
        sens.append(float(elements[2]))
        err_naive.append(float(elements[3]))
        sens_naive.append(float(elements[4]))
        
xdata = np.array(err)
ydata = np.array(sens)

err_naive = err_naive[1:4]
sens_naive = sens_naive[1:4]
xdata_naive = np.array(err_naive)
ydata_naive = np.array(sens_naive)

fig, ax = plt.subplots()

s = ax.plot(xdata, ydata, label='SPARTA')
s2 = ax.scatter(xdata_naive, ydata_naive, label='Separation Based on\n Number of Mismatches', color = 'r', s= 100)

#ax.plot(xdata, f_of_x)
#ax.plot(xdata, xdata, 'grey', label = 'expected')
ax.grid(True)
ax.set_xlabel('Error Rate (incorrect / classified)', fontsize=15)
ax.set_ylabel('Sensitivity (correct / total reads)', fontsize=15)
ax.set_title(title, fontsize=20)
ax.set_xlim(0, 0.0206280)
ax.set_ylim(0, 0.24)

plt.legend(loc='lower right')

labels = ['1 Mismatch\nDifference\nThreshold','2 Mismatch\nDifference\nThreshold','3 Mismatch\nDifference\nThreshold']
ix = 0
for label, x, y in zip(labels, err_naive, sens_naive):
    xytext = (125, 20) if ix > 0 else (60, -70)    
    plt.annotate(
        label, xy = (x, y), xytext = xytext,
        textcoords = 'offset points', ha = 'right', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'skyblue', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    ix += 1
fig.tight_layout()
plt.show()
#plt.savefig('all_mismatch_probs.png',format='png')
