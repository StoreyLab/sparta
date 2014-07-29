# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 15:15:19 2014

@author: peter
"""

# based on example from:
# http://www.walkingrandomly.com/?p=5215
import os
import math
import numpy as np
from scipy.optimize import curve_fit
#import matplotlib as mpl
#mpl.use('GTK3Agg')
from matplotlib import pyplot as plt

filename = '../final_analyses/skelly_mismatch_prob.txt'
title = 'Skelly et al Paired End Reads:\nPhred Scores vs Actual Mismatch Probability'

phred2 = []
mismatch_prob2 = []

with open (filename, 'r') as inf:
    for l in inf:
        elements = l.rstrip().split('\t')
        phred2.append(float(elements[0]))
        mismatch_prob2.append(float(elements[1]))

#phred2 = phred2[1:]
#mismatch_prob2 = mismatch_prob[1:]
assert len(phred2) == len(mismatch_prob2)

xdata2 = np.array([pow(10, (x-33)/-10.0) for x in phred2])
ydata2 = np.array(mismatch_prob2)
#xdata = np.array([x-33 for x in phred])
#ydata = np.array([math.log10(x)*-10.0 for x in mismatch_prob])

def func(x, u1, u2):
    return np.power(10, np.log10(x)*u1 + u2)

'''
popt, pcov = curve_fit(func,xdata,ydata,p0=(1,0.001))

res =  sum((ydata - func(xdata,popt[0], popt[1]))**2)
mean_y = np.mean(ydata)
mean_y_array = np.array([mean_y]*len(phred))
tot = sum((ydata - mean_y_array)**2)
R2 = 1-(res/tot)
print ('R^2 equals: {}'.format(R2))

f_of_x = [func(x, popt[0], popt[1]) for x in xdata]
'''
fig, ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')
s2 = ax.scatter(xdata2, ydata2, color='r', label='Skelly et al', s=40)
#ax.plot(xdata, f_of_x)
ax.plot(xdata, xdata, 'grey', label = 'expected')
ax.grid(True)
ax.set_xlabel('quality score mismatch prob', fontsize=20)
ax.set_ylabel('calculated mismatch prob', fontsize=20)
ax.set_title(title, fontsize=20)
ax.set_xlim(0.0001, 0.7)
ax.set_ylim(0.0001, 0.7)
plt.legend(loc='lower right')

fig.tight_layout()
plt.savefig('all_mismatch_probs.png',format='png')
