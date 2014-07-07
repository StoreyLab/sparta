# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 15:15:19 2014

@author: peter
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('Qt4Agg')
from matplotlib import pyplot as plt

actual_phred = [6.955214743808923, 7.849136533789498, 9.768489425923496, 11.121652762022357, 10.984179814741575, 14.215415131440974, 16.650038244190352, 15.731433334665972, 15.02974323481132, 14.0716337728756, 16.229146335386428, 14.297593028062598, 15.80976492643922, 16.116055325778238, 18.8730144307791, 20.551791958191053, 20.77418284942069, 20.98546787616192, 18.18291345409805, 19.597385558452117, 21.365618145177887, 25.760610496949955, 18.845960792973276, 25.730625876248787, 22.200009726648542, 21.027711306853135, 21.025292833867706, 22.866453246047563, 20.187883487692034, 20.742748066414958, 21.11621109905019, 30.551677811035596, 30.08745842343016]#, 35.64460722276653, 36.07531054192317, 36.75821199662197, 37.04523586956818]
qual_score_phred = range(0, len(actual_phred))

xdata = np.array(qual_score_phred)
ydata = np.array(actual_phred)

def func(x, u1, u2, u3):
    return u1*(x**2) + u2*x + u3
    
popt, pcov = curve_fit(func,xdata,ydata,p0=(0,0.75, 8))

res =  sum((ydata - func(xdata,popt[0], popt[1], popt[2]))**2)
mean_y = np.mean(ydata)
mean_y_array = np.array([mean_y]*len(actual_phred))
tot = sum((ydata - mean_y_array)**2)
R2 = 1-(res/tot)
print ('R^2 equals: {}'.format(R2))

f_of_x = [func(x, popt[0], popt[1], popt[2]) for x in xdata]

fig, ax = plt.subplots()
ax.scatter(xdata, ydata)
ax.plot(xdata, f_of_x)
ax.plot(xdata, xdata, 'r')
ax.grid(True)
ax.set_xlim(-1, 51)
ax.set_ylim(-1, 51)
fig.tight_layout()
fig.show()