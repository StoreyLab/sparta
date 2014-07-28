# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 20:58:11 2014

@author: peter
"""
from __future__ import print_function
import argparse
from compatibility import izip
import pysam
import re
import sys


# parse sparta program input.
def parseargs():
    
    parser = argparse.ArgumentParser(description="parse args")
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('samfile1', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)
    parser.add_argument('samfile2', nargs='?', type = str, help='path to samfile 2', default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type = str, help='path to outfile', default=sys.stdin)

    args = parser.parse_args()
    return args
    
def sort_samfiles(samfile1, samfile2, outfile):    

    sam1 = pysam.Samfile(samfile1)
    sam2 = pysam.Samfile(samfile2)
    
    n_regex = re.compile('\d[nN]\d')
    MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)")

    def count_errors(aligned_read):
        if not aligned_read.is_unmapped:               
            err = re.findall(MD_REGEX, aligned_read.opt("MD"))
            len_err = len(err) - len(re.findall(n_regex, aligned_read.opt("MD")))
            assert len_err >= 0
            return err, len_err
        else:
            return None, None    

    with open(outfile,'w') as opf:
        for aligned1, aligned2 in izip(sam1, sam2):
            
            assert aligned1.qname == aligned2.qname
            
            err1, len_err1 = count_errors(aligned1)
            err2, len_err2 = count_errors(aligned2)
            
            if aligned1.is_unmapped and aligned2.is_unmapped:

                print('0',file=opf)
                    
            elif (not aligned1.is_unmapped) and aligned2.is_unmapped:

                print('1',file=opf)            
            
            elif aligned1.is_unmapped and (not aligned2.is_unmapped):
                
                print('2',file=opf)

            
            elif len_err1 == len_err2:

                print('0',file=opf)

            else:
                
                if len_err1 > len_err2:
                    print('2',file=opf)
                else:
                    print('1',file=opf)
        
if __name__ == '__main__':
    
    # Get command line args
    args = parseargs()
    samfile1 = args.samfile1
    samfile2 = args.samfile2
    outfile = args.outfile

    sort_samfiles(samfile1, samfile2, outfile)