# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 14:11:25 2014

@author: peter
"""
from __future__ import print_function
import sys, os
import fileinput
import argparse

def parseargs():
    
    parser = argparse.ArgumentParser(description='')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('name', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    
    args = parser.parse_args()
    return args

cutoffs = [0.5000000000001, 0.5100000000001, 0.5200000000001, 0.5300000000001001, 0.5400000000001001,
 0.5500000000001001, 0.5600000000001, 0.5700000000001, 0.5800000000001, 0.5900000000001,
 0.6000000000001, 0.6100000000001, 0.6200000000001, 0.6300000000001, 0.6400000000001,
 0.6500000000001, 0.6600000000001001, 0.6700000000001001, 0.6800000000001001, 0.6900000000001001,
 0.7000000000001001, 0.7100000000001, 0.7200000000001, 0.7300000000001, 0.7400000000001,
 0.7500000000001, 0.7600000000001, 0.7700000000001, 0.7800000000001001, 0.7900000000001,
 0.8000000000001, 0.8100000000001, 0.8200000000001, 0.8300000000001, 0.8400000000001,
 0.8500000000001, 0.8600000000001, 0.8700000000001, 0.8800000000001, 0.8900000000001,
 0.9000000000001, 0.9100000000001001, 0.9200000000001001, 0.9300000000001001, 0.9400000000001001,
 0.9500000000001001, 0.9600000000001001, 0.9700000000001, 0.9800000000001, 0.9900000000001,
 0.9930000000001,0.9960000000001, 0.999, 0.9993000000001, 0.9996000000001, 0.9999, 0.9999300000001,
 0.9999600000001, 0.99999, 0.999993, 0.999996, 0.999999, 0.9999993, 0.9999996, 0.9999999, 1.0]

args = parseargs()
name = args.name

full_path = os.path.join('reports', name, 'supplementary_output.txt')

with open (full_path, 'r') as inf:

    header = inf.readline()
    
    classified = [0.0]*66
    total = 0.0
    
    for line in inf:
        
        for i in range(0,66):
        
        
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
