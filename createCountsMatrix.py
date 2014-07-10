# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 18:17:13 2014

@author: peter
"""
from __future__ import print_function
import os
import collections
from copy import copy

all_count_files = os.listdir('all_sample_analysis/count')
count_output_files = filter((lambda x: x[-1:] == 'o'), all_count_files)

counts_dict = collections.defaultdict(lambda: collections.defaultdict(int))
all_names = []

for count_output_file in count_output_files:
    
    name = copy(count_output_file[-10:-7]) + copy(count_output_file[:-11])
    all_names.append(name)    
    full_name = os.path.join('all_sample_analysis/count', count_output_file)
    
    with open(full_name, 'r') as ipf: 
        for line in ipf:
            fields = line.rstrip().split('\t')
            counts_dict[fields[0]][name] += int(fields[1])

all_genes = counts_dict.keys()
all_genes.sort()
all_names.sort()

with open('all_sample_analysis/count_matrix','w') as matrix_out:
    print(str.join('\t', all_names), file=matrix_out)
    
    for gene in all_genes:
        line_list = [gene]
        for name in all_names:
            line_list.append(str(counts_dict[gene][name]))
        print(str.join('\t', line_list), file=matrix_out)