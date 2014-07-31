# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 14:11:25 2014

@author: peter
"""
from __future__ import print_function
import sys, os
import argparse
import numpy as np
def parseargs():
    
    parser = argparse.ArgumentParser(description='')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('name', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    
    args = parser.parse_args()
    return args

cutoffs = [0.500000001]
cutoffs.extend([1.0 - 10**(-x) for x in list(np.arange(0.4, 10, 0.2))])
cutoffs.append(1.0)

args = parseargs()
name = args.name

full_path = os.path.join('reports', name, 'supplementary_output.txt')

with open (full_path, 'r') as inf:

    header = inf.readline()
    
    classified = [0.0]*50
    total = 0.0
    
    for line in inf:
        
        for i in range(0,50):
        
        
            elements = line.rstrip().split()
            prob_BY = float(elements[2]) if elements[2] != 'NA' else 0.5      
            prob_RM = 1.0-prob_BY
            #100 is to deal with the fact that one-sided mappings have NA as num_errs for unmapped

            if prob_BY > cutoffs[i] or prob_RM > cutoffs[i]:
                classified[i] += 1
            
        total += 1
    
    with open(os.path.join('cutoff_plot_data', name), 'w') as outf:
        for i in classified:
            frac = i/total
            print(frac, file=outf)
