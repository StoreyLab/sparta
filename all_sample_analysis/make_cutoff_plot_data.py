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

results = [[0]*54]*54

for name in run_names:
    
    full_path = os.path.join('reports', name, 'supplementary_output.txt')
    
    with open (full_path, 'r') as inf:

        header = inf.readline()
        
        classified = [0.0]*54
        total = 0.0
        
        for line in inf:
            
            for i in range(0,54):
            
            
                elements = line.rstrip().split()
                prob_BY = float(elements[2]) if elements[2] != 'NA' else 0.5      
                prob_RM = 1.0-prob_BY
                #100 is to deal with the fact that one-sided mappings have NA as num_errs for unmapped
    
                if prob_BY > cutoffs[i] or prob_RM > cutoffs[i]:
                    classified[i] += 1
                
            total += 1
        
        for i in range(0,54):
            
            results[run_ix][i] = classified / total

    run_ix += 1

with open ('cutoff_plot_data', 'w') as opf:
    for i in range(0,54):
        print (''.join(results[i]))
