# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 17:29:49 2014

@author: peter
"""
import os

if not os.path.exists('all_sample_analysis/split_by_genotypes'):
    os.makedirs('all_sample_analysis/split_by_genotypes')
    
if not os.path.exists('all_sample_analysis/reports'):
    os.makedirs('all_sample_analysis/reports')

if not os.path.exists('all_sample_analysis/samfiles'):
    os.makedirs('all_sample_analysis/samfiles')

if not os.path.exists('all_sample_analysis/out'):
    os.makedirs('all_sample_analysis/out')
    
if not os.path.exists('all_sample_analysis/err'):
    os.makedirs('all_sample_analysis/err')
    
input_files = os.listdir('/Genomics/grid/users/emilysn/ase/analysis/RNA-Seq/RUN_111213OB/split_fastq')

for fastq in input_files:
    name = fastq[:-6]
    out = os.path.join('all_sample_analysis/out', name)
    err = os.path.join('all_sample_analysis/err',name)
    os.system('qsub -o {} -e {} -cwd python runBowtieSeqsort.py {}'.format(out, err, fastq))