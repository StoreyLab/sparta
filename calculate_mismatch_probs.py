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
the SPARTA software by modifying the original quality scores to reflect the
actual observed scores at each phred score.
'''

import argparse
from Bio.Seq import MutableSeq
import os
import pysam
import re
import sys
import time

from util import compatibility_dict
from util import fix_read_mate_order
from util import izip

# default arguments
from sparta import default_output_dir
from sparta import default_pileup_height
from sparta import default_sample_every

# regex for the MD string that specifies errors from the reference.
# more information about the MD string: page 7 of http://samtools.github.io/hts-specs/SAMv1.pdf
MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)")

def parseargs():
    
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument('samfiles', nargs='+', type = str, help='input samfiles', default=[])
    parser.add_argument('-o', '--output_dir', nargs='?', type = str, help='directory to write output to', default=default_output_dir)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# take an aligned read, return the genomic seq EXCEPT for deletions (^)
def create_genome_seq(aligned):

    aligned_seq = aligned.seq if type(aligned.seq) == str else aligned.seq.decode('UTF-8')
        
    genome_seq = MutableSeq(aligned_seq)

    # see samtools documentation for MD string
    
    err = re.findall(MD_REGEX, aligned.opt("MD"))
    
    seq_ix = 0
    
    # step through sequence
    for matched_bases, curr_err in err:
        
        seq_ix += int(matched_bases)
        
        assert '^' not in curr_err
        assert curr_err != genome_seq[seq_ix]
        
        genome_seq[seq_ix] = curr_err
        seq_ix += 1
    
    if aligned.is_reverse:
        genome_seq.reverse_complement()

    return genome_seq

def add_to_pileup_dict(sams, aligned_read_set, pileup_dict):
    
    # sanity check that all the qnames (RNA read IDs) are the same
    for read in aligned_read_set:
        assert read.qname == aligned_read_set[0].qname

    if not True in [read.is_unmapped for read in aligned_read_set]:
        
        # all alignments mapped
        for read in aligned_read_set:
            
            for op, op_len in read.cigar:
                
                if op > 0 and op < 7:
                    # do not sample reads where there are insertions or deletions   
                    return
                    
            assert len(read.seq) == len(aligned_read_set[0].seq)
          
        # if aligned reads are reversed, we reverse them and hold on to that info.
        
        pos_dicts = [dict(read.aligned_pairs) for read in aligned_read_set]
        genome_seqs = [create_genome_seq(read) for read in aligned_read_set]
        qual = bytearray(aligned_read_set[0].qual)
        seq = MutableSeq(aligned_read_set[0].seq if type(aligned_read_set[0].seq) == str else aligned_read_set[0].seq.decode('UTF-8'))  
        if aligned_read_set[0].is_reverse:
            seq.reverse_complement()
            qual = qual[::-1]
        
        for genome_seq in genome_seqs:
            assert len(genome_seq) == len(seq)

        for i in range(0, len(seq)):
                        
            # need (chrom, pos, genome_seq[i]) tuples for each aligned_read
            chroms = [sam.getrname(a.tid) for sam, a in izip(sams, aligned_read_set)]
            positions = [d[i] if not a.is_reverse else d[len(seq) - i - 1] for d, a in zip(pos_dicts, aligned_read_set)]
            genome_seq_i = [g[i] for g in genome_seqs]
            
            genomic_locs = tuple(zip(chroms, positions, genome_seq_i))
            
            pileup_dict[genomic_locs][seq[i]][qual[i]] += 1

def create_mismatch_prob_dict(samfiles, output_dir = default_output_dir, paired_end=False, pileup_height=default_pileup_height, sample_every=default_sample_every):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    sams = tuple(pysam.Samfile(i) for i in samfiles)
    
    pileup_dict = compatibility_dict(lambda: compatibility_dict(lambda: compatibility_dict(int)))
            
    if not paired_end:
        
        for i, aligned_read_set in enumerate(izip(*sams)):

            # sample every 10th read        
            if i % sample_every != 0:
                continue
            
            add_to_pileup_dict(sams, aligned_read_set, pileup_dict)
            
    else:
        
        # the main difference with paired end reads is that bowtie should output
        # a read-mate pair, followed by another read-mate pair, in the same order
        # regardless of whether it is the genome1 mapping or the genome2 mapping
        # but, for a given read-mate pair we don't necessarily know if we are 
        # getting the read we want or its mate.
        # so we have to check, and switch them if necessary
        zipped_samfiles = izip(*sams)
        for i, aligned_read_set in enumerate(zipped_samfiles):

            # take aligned pairs in sets of 2 to also get the mate
            # don't forget to assign aligned_pair to next tuple before end of iter
            aligned_read_mate_set = next(zipped_samfiles)       

            # sample rate      
            if i % sample_every != 0:
                aligned_read_set = aligned_read_mate_set
                continue

            aligned_read_set, aligned_read_mate_set = fix_read_mate_order(aligned_read_set, aligned_read_mate_set)
            

            for aligned_read_set_generic in [aligned_read_set,aligned_read_mate_set]:
                
                add_to_pileup_dict(sams, aligned_read_set_generic, pileup_dict)
            
            # skip the next pair because we already grabbed it as the mate
            aligned_read_set = aligned_read_mate_set
    
    for sam in sams:
        sam.close()
    
    quality_score_match_counter = compatibility_dict(int)
    quality_score_mismatch_counter = compatibility_dict(int)
    # transition matrix at [qual][X][Y] stores, for a given quality score qual,
    # the number of observed transions from X to Y
    transition_matrix = compatibility_dict(lambda: compatibility_dict(int))

    with open(os.path.join(output_dir, 'pileup_counts'), 'w') as outfile:
        
        print('genomic_coordinates\tbase_count[A]\tbase_count[C]\tbase_count[G]\tbase_count[T]\tbase_count[N]', file=outfile)
           
        for genomic_locs, nuc_to_qual_dict in pileup_dict.items():
            # iterate through genome coordinates and their nested dictionaries
            
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
            consensus = 'N'
    
            if total_bases >= pileup_height:
                # if more than 75% of reads agree on a base, and if both genomic sequences
                # in that position is that same base, then it is the consensus.       
                for base, count in base_count.items():
                    if count > cutoff * total_bases and False not in [base == loc[2] for loc in genomic_locs]:
                        consensus = base
                                
                log_msg = '{}\t{}\t{}\t{}\t{}\t{}'.format(genomic_locs, base_count['A']+base_count['a'], base_count['C']+base_count['c'],
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
                            transition_matrix[consensus][nuc] += num
                        
    mismatch_prob_dict = {}
    mismatch_prob_total_values = {}
    transition_prob_dict = {}
    
    for qual, mismatch_count in quality_score_mismatch_counter.items():

        mismatch_prob = mismatch_count * 1.0 / ((mismatch_count + quality_score_match_counter[qual])*1.0)
        mismatch_prob_dict[qual] = mismatch_prob
        mismatch_prob_total_values[qual] = mismatch_count + quality_score_match_counter[qual]
        
        for base1, base_to_count_dict in transition_matrix.items():
            
            for base2, num_trans_base1_base2 in base_to_count_dict.items():
                
                # we have a qual, base1, base2, and a count of transitions
                # from base1 to base2 at that qual
                total_trans_from_base1 = sum(transition_matrix[base1].values())
                transition_prob_dict[(base1, base2)] = num_trans_base1_base2 / (total_trans_from_base1 * 1.0)
    
    # For each phred, print observed probability of mismatch and number of bases observed in creating that probability
    if mismatch_prob_dict and mismatch_prob_total_values:
        
        with open(os.path.join(output_dir, 'mismatch_prob_info.txt'), 'w') as outputfile:
            
            for k in mismatch_prob_dict.keys():
                print ('{}\t{}\t{}'.format(k,mismatch_prob_dict[k],mismatch_prob_total_values[k]), file=outputfile)

        with open(os.path.join(output_dir, 'transition_prob_info.txt'), 'w') as outputfile:
            
            for k, v in transition_prob_dict.items():
                base1, base2 = k
                print ('{}\t{}\t{}'.format(base1, base2, v), file=outputfile)   
                
    return mismatch_prob_dict, mismatch_prob_total_values, transition_prob_dict

# main logic
# call compare_mappings() on samfile1 and samfile2 from standard input 
def main():
    
    t1 = time.time()
    # get command line args
    args = parseargs()
    samfiles = args.samfiles
    output_dir = args.output_dir
    mismatch_prob_dict, mismatch_prob_total_values, transition_prob_dict = create_mismatch_prob_dict(samfiles, output_dir)  
    
    t2 = time.time()
    print('TOTAL TIME: {}'.format(t2-t1))
    
if __name__ == '__main__':
    main()
