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
'''
title = 'E1a12 Quality Score vs. Calculated Probability of Mismatch'
actual_phred = ([8.281213204366239, 6.955214743808923, 7.849136533789498, 9.768489425923496, 11.121652762022357,
                 10.984179814741575, 14.215415131440974, 16.650038244190352, 15.731433334665972, 15.02974323481132,
                 14.0716337728756, 16.229146335386428, 14.297593028062598, 15.80976492643922, 16.116055325778238,
                 18.8730144307791, 20.551791958191053, 20.77418284942069, 20.98546787616192, 18.18291345409805,
                 19.597385558452117, 21.365618145177887, 25.760610496949955, 18.845960792973276, 25.730625876248787,
                 22.200009726648542, 21.027711306853135, 21.025292833867706, 22.866453246047563, 20.187883487692034,
                 20.742748066414958, 21.11621109905019, 30.551677811035596, 30.08745842343016, 35.64460722276653,
                 36.07531054192317, 36.75821199662197, 37.04523586956818])
qual_score_phred = ([2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                     26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41])
'''
'''
title = 'G1a60 Quality Score vs. Calculated Probability of Mismatch'
actual_phred = ([8.40460777100269, 6.493495727635167,7.268787670256572, 9.699990810481605, 10.8119836271688,
                 10.992650558342117, 14.300852984671135, 16.732392004141673, 15.348912832040671, 14.835303769813613,
                 13.647407386796772, 16.350420363493818, 14.125568052424526, 16.593114347465935, 16.069558075087954,
                 19.269701763841173, 20.844657334837102, 20.69031615629621, 20.799670835728246, 18.75012771115153,
                 19.937305913329173, 21.891599899473416, 26.269940366821366, 19.46389184978097, 26.267873164006318,
                 22.705892093195818, 21.61744667758954, 21.28332916941215, 23.70383765148049, 20.7292170246182,
                 20.89831046923191, 21.79857911787673, 30.860509096348267, 30.01946791319427, 35.50253517998259,
                 35.929220007136, 36.87472808781002, 36.839017249308526])
qual_score_phred = ([2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                     26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41])
'''
title = 'simulated read mismatch prob'
#filename = 'dwgsim/mismatch_prob_info.txt'
dir_list = os.listdir('all_sample_analysis/mismatch_prob_reports')
dir_list.remove('Unassigned_lane0_index')
dir_list.remove('Unassigned_lane1_index')
path_list = [os.path.join('all_sample_analysis/mismatch_prob_reports',x,'mismatch_prob_info.txt') for x in dir_list]

phred = []
mismatch_prob = []

for path in path_list:
    
    with open(path,'r') as infile:
        for line in infile:
            e = line.rstrip().split('\t')
            phred.append(float(e[0]))
            mismatch_prob.append(float(e[1]))

'''
with open (filename, 'r') as inf:
    for l in inf:
        elements = l.rstrip().split('\t')
        phred.append(float(elements[0]))
        mismatch_prob.append(float(elements[1]))
        
'''
assert len(phred) == len(mismatch_prob)

xdata = np.array([pow(10, (x-33)/-10.0) for x in phred])
ydata = np.array(mismatch_prob)
#xdata = np.array([x-33 for x in phred])
#ydata = np.array([math.log10(x)*-10.0 for x in mismatch_prob])

def func(x, u1, u2, u3):
    return np.log10(u1*(x**2) + u2*x + u3)
    
popt, pcov = curve_fit(func,xdata,ydata,p0=(0,0.75, 8))

res =  sum((ydata - func(xdata,popt[0], popt[1], popt[2]))**2)
mean_y = np.mean(ydata)
mean_y_array = np.array([mean_y]*len(phred))
tot = sum((ydata - mean_y_array)**2)
R2 = 1-(res/tot)
print ('R^2 equals: {}'.format(R2))

f_of_x = [func(x, popt[0], popt[1], popt[2]) for x in xdata]

fig, ax = plt.subplots()
ax.scatter(xdata, ydata)
ax.plot(xdata, f_of_x)
ax.plot(xdata, xdata, 'r')
ax.grid(True)
ax.set_xlabel('quality score mismatch prob', fontsize=20)
ax.set_ylabel('calculated mismatch prob', fontsize=20)
ax.set_title(title)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0, 0.2)
ax.set_ylim(0, 0.2)

fig.tight_layout()
plt.savefig('all_mismatch_probs.png',format='png')
