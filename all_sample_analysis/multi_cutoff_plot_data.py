# -*- coding: utf-8 -*-
"""
Create ROC data for sparta

Created on Mon Jul 28 11:48:22 2014

@author: peter
"""
from __future__ import print_function
import os

cutoff_init = 0.5000000000001
cutoffs = [cutoff_init + float(x)/100.0 for x in range(0,50)]
cutoffs.extend([1.0 - 10**(-x) for x in range(3, 6)])
cutoffs.append(1.0)

results = [[0]*54]*54

run_names = os.listdir('reports')
run_names.remove('Unassigned_lane0_index')
run_names.remove('Unassigned_lane1_index')

run_ix = 0

os.mkdir ('cutoff_plot_data')

for name in run_names:
    
    os.system('qsub -cwd -o cutoff_plot_data/{}.o -e cutoff_plot_data/{}.e python cutoff_subroutine.py {}'.format(name, name, name))    

