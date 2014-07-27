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
from compatibility import compatibility_dict
from compatibility import izip
from compatibility import rev_comp
from pprint import pprint
import pysam
import re
import sys
import time
from util import fix_read_mate_order
            
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
    
    genome_seq = list(copy.copy(aligned.seq)) if type(aligned.seq) == str else list(copy.copy(aligned.seq.decode('UTF-8')))
    
    # see samtools documentation for MD string
    
    err = re.findall(MD_REGEX, aligned.opt("MD"))
    
    seq_ix = 0
    
    # step through sequence
    for matched_bases, curr_err in err:
        
        seq_ix += int(matched_bases)
        
        if '^' not in curr_err:
            
            genome_seq[seq_ix] = curr_err
            seq_ix += 1

    return ''.join(genome_seq)

def add_to_pileup_dict(sam1, sam2, aligned1, aligned2, pileup_dict):
    
    assert aligned1.qname == aligned2.qname
    
    if not aligned1.is_unmapped and not aligned2.is_unmapped:
        # both alignments mapped
    
        assert len(aligned1.seq) == len(aligned2.seq)
        
        same_strand = (aligned1.is_reverse == aligned2.is_reverse)
            
        pos_dict1 = dict(aligned1.aligned_pairs)
        pos_dict2 = dict(aligned2.aligned_pairs)
            
        qual1 = bytearray(aligned1.qual)
        qual2 = bytearray(aligned2.qual) if same_strand else bytearray(aligned2.qual)[::-1]
        
        genome_seq1 = create_genome_seq(aligned1)
        genome_seq2 = create_genome_seq(aligned2) if same_strand else rev_comp(create_genome_seq(aligned2))
        
        # The UTF-8 stuff is Python 3 shenanigans (Py3 aligned.seqs are bytearrays, != strings)                
        aligned1_seq = aligned1.seq if type(aligned1.seq) == str else aligned1.seq.decode('UTF-8')
        aligned2_seq_temp = aligned2.seq if same_strand else rev_comp(aligned2.seq)
        aligned2_seq = aligned2_seq_temp if type(aligned2_seq_temp) == str else aligned2_seq_temp.decode('UTF-8')
        
        assert aligned1_seq == aligned2_seq
        assert qual1 == qual2
        
        for i in range(0, len(aligned1.seq)):
                        
            chrom1 = sam1.getrname(aligned1.tid)
            pos1 = pos_dict1[i]      
            chrom2 = sam2.getrname(aligned2.tid)
            # pos_dict2 is the only thing that wasn't reversed, so if not same strand
            # we simply reverse how we index into it (len - i - 1) to go right-to-left
            pos2 = pos_dict2[i] if same_strand else pos_dict2[len(aligned2.seq) - i - 1]
            
            pileup_dict[(chrom1, chrom2, pos1, pos2, genome_seq1[i], genome_seq2[i])][aligned1_seq[i]][qual1[i]] += 1

def create_mismatch_prob_dict(samfile1, samfile2, genome1_name, genome2_name, outfile_name, paired_end, pileup_height=20, sample_every=10):
    
    sam1 = pysam.Samfile(samfile1)
    sam2 = pysam.Samfile(samfile2)
    
    pileup_dict = compatibility_dict(lambda: compatibility_dict(lambda: compatibility_dict(int)))
    
    i = 0
        
    if not paired_end:
        
        for aligned1, aligned2 in izip(sam1, sam2):

            # sample every 10th read        
            if i % sample_every != 0:
                i += 1
                continue
            i += 1
            
            add_to_pileup_dict(sam1, sam2, aligned1, aligned2, pileup_dict)

                    
    else:
        
        # the main difference with paired end reads is that bowtie should output
        # a read-mate pair, followed by another read-mate pair, in the same order
        # regardless of whether it is the genome1 mapping or the genome2 mapping
        # but, for a given read-mate pair we don't necessarily know if we are 
        # getting the read we want or its mate.
        # so we have to check, and switch them if necessary
        zipped_samfiles = izip(sam1, sam2)
        for aligned_pair in zipped_samfiles:

            # take aligned pairs in sets of 2 to also get the mate
            # don't forget to assign aligned_pair to next tuple before end of iter
            aligned1, aligned2 = aligned_pair
            next_tuple = next(zipped_samfiles)       
            aligned1_mate, aligned2_mate = next_tuple

            # sample every 10th read + mate combo       
            if i % sample_every != 0:
                aligned_pair = next_tuple
                i += 1
                continue
            i += 1

            aligned1, aligned2, aligned1_mate, aligned2_mate = fix_read_mate_order(aligned1, aligned2, aligned1_mate, aligned2_mate)
            

            for a1, a2 in [(aligned1, aligned2),(aligned1_mate, aligned2_mate)]:
                
                add_to_pileup_dict(sam1, sam2, a1, a2, pileup_dict)
            
            # skip the next pair because we already grabbed it as the mate
            aligned_pair = next_tuple
        
    sam1.close()
    sam2.close()
    
    quality_score_match_counter = compatibility_dict(int)
    quality_score_mismatch_counter = compatibility_dict(int)
        
    with open(outfile_name, 'w') as outfile:
        
        print('#min_pileup_height={}, sample_every={}'.format(pileup_height, sample_every), file=outfile)
        print('chrom1\tchrom2\tpos1\tpos2\tgenome1_seq(i)\tgenome2_seq(i)\tbase_count[A]\tbase_count[C]\tbase_count[G]\tbase_count[T]\tbase_count[N]', file=outfile)
           
        for coordinate_pair, nuc_to_qual_dict in pileup_dict.items():
            # iterate through genome coordinate pairs and their nested dictionaries
            
            # want to know: what is the consensus base at this coord pair?
            # keep a count of how many bases seen at this coord pair
            base_count = compatibility_dict(int)
            
            for nuc, qual_dict in nuc_to_qual_dict.items():
                # for a given coordinate pair, iterate through its nucleotides and
                # the dictionary that pairs qualities to number of matches
                # iterate through nucleotides and their quality dictionaries
                for qual, num in qual_dict.items():
                    # iterate over qualities and number of matches
                    
                    # AT THIS POINT we have:
                    # a coordinate pair, a nucleotide at that coordinate pair,
                    # a quality score, and the number of matches at that quality score
                    base_count[nuc] += 1
            
            total_bases = sum([v for k,v in base_count.items()])
            cutoff = 0.75
            chrom1, chrom2, pos1, pos2, genome1_seq_i, genome2_seq_i = coordinate_pair
            consensus = 'N'
    
            if total_bases >= 20:
                # if more than 75% of reads agree on a base, and if both genomic sequences
                # in that position is that same base, then it is the consensus.       
                for base, count in base_count.items():
                    if count > cutoff * total_bases and base == genome1_seq_i and base == genome2_seq_i:
                        consensus = base
                                
            if total_bases >= pileup_height:
                log_msg = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom1, chrom2,
                pos1, pos2, genome1_seq_i, genome2_seq_i, base_count['A']+base_count['a'], base_count['C']+base_count['c'],
                base_count['G']+base_count['g'],base_count['T']+base_count['t'], base_count['N']+base_count['n'])
                print(log_msg, file=outfile)
            
            for nuc, qual_dict in nuc_to_qual_dict.items():
                # for a given coordinate pair, iterate through its nucleotides and
                # the dictionary that pairs qualities to number of matches
                # iterate through nucleotides and their quality dictionaries

                for qual, num in qual_dict.items():
                    # iterate over qualities and number of matches
                    
                    # AT THIS POINT we have:
                    # a coordinate pair, a nucleotide at that coordinate pair,
                    # a quality score, and the number of matches at that quality score          
                    
                    if consensus != 'N':
                        if nuc == consensus: 
                            quality_score_match_counter[qual] += num
                        else:
                            quality_score_mismatch_counter[qual] += num
                        
    mismatch_prob_dict = {}
    mismatch_prob_total_values = {}
    
    for qual, mismatch_count in quality_score_mismatch_counter.items():

        mismatch_prob = mismatch_count * 1.0 / ((mismatch_count + quality_score_match_counter[qual])*1.0)
        mismatch_prob_dict[qual] = mismatch_prob
        mismatch_prob_total_values[qual] = mismatch_count + quality_score_match_counter[qual]
        
    return mismatch_prob_dict, mismatch_prob_total_values

# main logic
# call compare_mappings() on samfile1 and samfile2 from standard input 
def main():
    '''
    t1 = time.time()
    # get command line args
    args = parseargs()
    samfile1 = args.samfile1
    samfile2 = args.samfile2
    genome1_name = args.name1
    genome2_name = args.name2
    
    # compare mappings between samfiles
    pileup_dict = create_mismatch_prob_dict(samfile1, samfile2, genome1_name, genome2_name)
    
    t2 = time.time()
    pprint(pileup_dict)
    print('TOTAL TIME: {}'.format(t2-t1))
    return 0
    '''
if __name__ == '__main__':
    main()
