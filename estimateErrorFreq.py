#! /usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Created on Mon Jun 23 16:21:32 2014

@author: Peter Edge
"""
desc = '''
This file performs analysis of two SAM files in order to estimate the number
of RNAseq reads with errors to the original sequence. This will play a role in
the SeqSorter software by modifying the original quality scores to reflect the
actual observed scores at each phred score.
'''

import argparse
import copy
import collections
from pprint import pprint
import pysam
import re
import sys
import time

# regex for the MD string that specifies errors from the reference.
# more information about the MD string: page 7 of http://samtools.github.io/hts-specs/SAMv1.pdf
MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)")

def parseargs():
    
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument('samfile1', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)
    parser.add_argument('samfile2', nargs='?', type = str, help='path to samfile 2', default=sys.stdin)
    # optional nicknames for the genomes used in the two samfiles (recommended)    
    parser.add_argument('-n1', '--name1', nargs='?', type = str, help='name for genome 1 (reference for samfile 1)', default='genome1')
    parser.add_argument('-n2', '--name2', nargs='?', type = str, help='name for genome 2 (reference for samfile 2)', default='genome2')
    # this is an optional verbose outfile with useful values for statistical analysis

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# take an aligned read, return the genomic seq EXCEPT for deletions (^)
def create_genome_seq(aligned):
    
    genome_seq = list(copy.copy(aligned.seq))
    
    # see samtools documentation for MD string
    err = re.findall(MD_REGEX, aligned.opt("MD"))
    
    seq_ix = 0
    
    # step through sequence
    for matched_bases, curr_err in err:
        
        seq_ix += int(matched_bases)
        
        if '^' not in curr_err:
            
            genome_seq[seq_ix] = curr_err
            seq_ix += 1
            
    return genome_seq 
    
def create_mismatch_prob_dict(samfile1, samfile2, genome1_name, genome2_name):
    
    sam1 = pysam.Samfile(samfile1)
    sam2 = pysam.Samfile(samfile2)
    
    results = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))

    for aligned1, aligned2 in zip(sam1, sam2):
        
        assert aligned1.qname == aligned2.qname
        
        if not aligned1.is_unmapped and not aligned2.is_unmapped:
            
            assert len(aligned1.seq) == len(aligned2.seq)
            
            # both alignments mapped
            pos_dict1 = dict(aligned1.aligned_pairs)
            pos_dict2 = dict(aligned2.aligned_pairs)
            
            qual = bytearray(aligned1.qual)
            genome_seq1 = create_genome_seq(aligned1)
            genome_seq2 = create_genome_seq(aligned2)
            
            for i in range(0, len(aligned1.seq)):
                            
                if genome_seq1[i] == genome_seq2[i]:
                    
                    chrom1 = sam1.getrname(aligned1.tid)
                    pos1 = pos_dict1[i]      
                    chrom2 = sam2.getrname(aligned2.tid)
                    pos2 = pos_dict2[i]
                    
                    results[(chrom1, chrom2, pos1, pos2, genome_seq1[i])][aligned1.seq[i]][qual[i]] += 1
    
    quality_score_match_counter = collections.defaultdict(int)
    quality_score_mismatch_counter = collections.defaultdict(int)
    
    for coordinate_pairs, nuc_to_qual_dict in results.iteritems():
        # iterate through genome coordinate pairs and their nested dictionaries
        
        # want to know: what is the consensus base at this coord pair?
        # keep a count of how many bases seen at this coord pair
        base_count = collections.defaultdict(int)
        
        for nuc, qual_dict in nuc_to_qual_dict.iteritems():
            # for a given coordinate pair, iterate through its nucleotides and
            # the dictionary that pairs qualities to number of matches
            # iterate through nucleotides and their quality dictionaries
            for qual, num in qual_dict.iteritems():
                # iterate over qualities and number of matches
                
                # AT THIS POINT we have:
                # a coordinate pair, a nucleotide at that coordinate pair,
                # a quality score, and the number of matches at that quality score
                base_count[nuc] += 1
        
        total_bases = sum([v for k,v in base_count.iteritems()])
        cutoff = 0.75
        consensus = 'N'

        if total_bases >= 20:
            # if more than 75% of reads agree on a base, and if the genomic sequence
            # in that position is that same base, then it is the consensus.       
            for base, count in base_count.iteritems():
                if count > cutoff * total_bases and base == coordinate_pairs[4]:
                    consensus = base
               
        if consensus != 'N':
            for nuc, qual_dict in nuc_to_qual_dict.iteritems():
                # for a given coordinate pair, iterate through its nucleotides and
                # the dictionary that pairs qualities to number of matches
                # iterate through nucleotides and their quality dictionaries
                for qual, num in qual_dict.iteritems():
                    # iterate over qualities and number of matches
                    
                    # AT THIS POINT we have:
                    # a coordinate pair, a nucleotide at that coordinate pair,
                    # a quality score, and the number of matches at that quality score
                    if nuc == consensus: 
                        quality_score_match_counter[qual] += num
                    else:
                        quality_score_mismatch_counter[qual] += num
                    
    mismatch_prob_dict = {}
    
    for qual, mismatch_count in quality_score_mismatch_counter.iteritems():

        mismatch_prob = mismatch_count * 1.0 / ((mismatch_count + quality_score_match_counter[qual])*1.0)
        mismatch_prob_dict[qual] = mismatch_prob
        
    pprint(mismatch_prob_dict)
    
    return mismatch_prob_dict

# main logic
# call compare_mappings() on samfile1 and samfile2 from standard input 
def main():
    t1 = time.time()
    # get command line args
    args = parseargs()
    samfile1 = args.samfile1
    samfile2 = args.samfile2
    genome1_name = args.name1
    genome2_name = args.name2
    
    # compare mappings between samfiles
    results = create_mismatch_prob_dict(samfile1, samfile2, genome1_name, genome2_name)
    
    t2 = time.time()
    pprint(results)
    print('TOTAL TIME: {}'.format(t2-t1))
    return 0

if __name__ == '__main__':
    main()