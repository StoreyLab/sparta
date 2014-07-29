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
from util import fix_read_mate_order

# parse sparta program input.
def parseargs():
    
    parser = argparse.ArgumentParser(description="parse args")
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('samfile1', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)
    parser.add_argument('samfile2', nargs='?', type = str, help='path to samfile 2', default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type = str, help='path to outfile', default=sys.stdin)
    parser.add_argument('-pe', '--paired_end', nargs='?', type = int, help='set this flag to specify that reads are paired end (default: False)', const=1)

    args = parser.parse_args()
    return args
    
def sort_samfiles(samfile1, samfile2, paired_end, outfile):    

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
        if not paired_end:
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
        else: # paired end
            
            zipped_samfiles = izip(sam1, sam2)
            
            for aligned_pair in zipped_samfiles:
                
                aligned1, aligned2 = aligned_pair
                
                # take aligned pairs in sets of 2 to also get the mate
                # don't forget to assign aligned_pair to next tuple before end of iter
                next_tuple = next(zipped_samfiles)       
                aligned1_mate, aligned2_mate = next_tuple
                
                aligned1, aligned2, aligned1_mate, aligned2_mate = fix_read_mate_order(aligned1, aligned2, aligned1_mate, aligned2_mate)
                
                err_a1, len_err_a1 = count_errors(aligned1)
                err_a1_mate, len_err_a1_mate = count_errors(aligned1_mate)
                err_a2, len_err_a2 = count_errors(aligned2)
                err_a2_mate, len_err_a2_mate = count_errors(aligned2_mate)

                if (not aligned1.is_unmapped and not aligned1_mate.is_unmapped and
                    not aligned2.is_unmapped and not aligned2_mate.is_unmapped):
                    
                    # both the read and its mate mapped to both genomes
                    if err_a1 == err_a2 and err_a1_mate == err_a2_mate:
                        pass
                    else:
                        
                        sort_fate1 = None
                        if err_a1 == err_a2:
                            sort_fate1 = '0'
                        elif err_a1 < err_a2:
                            sort_fate1 = '1'
                        else:
                            sort_fate1 = '2'
                            
                        sort_fate2 = None
                        if err_a1_mate == err_a2_mate:
                            sort_fate2 = '0'
                        elif err_a1_mate < err_a2_mate:
                            sort_fate2 = '1'
                        else:
                            sort_fate2 = '2'
                        
                        print('{}\t{}'.format(sort_fate1, sort_fate2), file=opf)
                # CRUCIAL STEP: skip the mate read
                aligned_pair = next_tuple
                
    sam1.close()
    sam2.close()
    
if __name__ == '__main__':
    
    # Get command line args
    args = parseargs()
    samfile1 = args.samfile1
    samfile2 = args.samfile2
    paired_end = args.paired_end
    outfile = args.outfile

    sort_samfiles(samfile1, samfile2, paired_end, outfile)