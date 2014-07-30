# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 14:11:25 2014

@author: peter
"""
from __future__ import print_function
import os
import fileinput


cutoff_init = 0.5000000000001
cutoffs = [cutoff_init + float(x)/100.0 for x in range(0,50)]
cutoffs.extend([1.0 - 10**(-x) for x in range(3, 6)])
cutoffs.append(1.0)

name = fileinput.input()

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
    
    with open(os.path.join('cutoff_plot_data', name)) as outf:
        for i in classified:
            frac = classified[i]/total
            print(frac, file=outf)