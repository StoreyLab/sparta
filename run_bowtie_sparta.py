# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 20:43:15 2014

@author: Peter Edge
"""

import argparse
import sys, os

def parseargs():
    
    parser = argparse.ArgumentParser('parse args')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('fastqfile', nargs='?', type = str, help='path to fastq', default=sys.stdin)
    args = parser.parse_args()
    return args

def run_seqsort(fastqfile):
    by_index = 'all_sample_analysis/index/BY_index'
    rm_index = 'all_sample_analysis/index/RM_index'
    
    raw_name = fastqfile[:-6]
    full_path = os.path.join('/Genomics/grid/users/emilysn/ase/analysis/RNA-Seq/RUN_111213OB/split_fastq', fastqfile)
    by_sam = os.path.join('all_sample_analysis/samfiles','{}_BY.sam'.format(raw_name))
    rm_sam = os.path.join('all_sample_analysis/samfiles','{}_RM.sam'.format(raw_name))
    os.system('bowtie2-align-s -x {} -U {} -S {}'.format(by_index, full_path, by_sam))
    os.system('bowtie2-align-s -x {} -U {} -S {}'.format(rm_index, full_path, rm_sam))
    report = os.path.join('all_sample_analysis/reports', raw_name)
    by_out = os.path.join('all_sample_analysis/split_by_genotypes', '{}_BY_sorted.sam'.format(raw_name))
    rm_out = os.path.join('all_sample_analysis/split_by_genotypes', '{}_RM_sorted.sam'.format(raw_name))
    os.system('python sparta.py {} {} -n1 BY -n2 RM -e -s1 {} -s2 {} -o {}'.format(by_sam, rm_sam, by_out, rm_out, report))    

    by_count_out = os.path.join('all_sample_analysis/count', '{}_BY_count.o'.format(raw_name))
    by_count_err = os.path.join('all_sample_analysis/count', '{}_BY_count.e'.format(raw_name))
    rm_count_out = os.path.join('all_sample_analysis/count', '{}_RM_count.o'.format(raw_name))
    rm_count_err = os.path.join('all_sample_analysis/count', '{}_RM_count.e'.format(raw_name))  
    
    os.system('htseq-count -i ID -s no -t gene {} genomes/BY.gff > {} 2> {}'.format(by_out,by_count_out, by_count_err))
    os.system('htseq-count -i ID -s no -t gene {} genomes/RM.gff > {} 2> {}'.format(rm_out, rm_count_out, rm_count_err))
    
if __name__ == '__main__':
    
    args = parseargs()
    fastqfile = args.fastqfile
    run_seqsort(fastqfile)