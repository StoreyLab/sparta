#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
'''
Created on Mon Jun 16 10:56:22 2014

@author: Peter Edge
adapted from compare_mappings() in Emily Nelson's script snp_pipeline.py
'''
desc = '''
This file contains the multimapped_read_sorter class, which can sort two SAM
alignment files which map the same RNAseq reads to two different parental
genomes. The RNAseq reads should map to either one parental allele type or the
other. This program sorts the reads based on the errors to each genome.
'''
'''
COMMAND LINE USAGE:
python2 SeqSorter.py <samfile1> <samfile2> <optional arguments>

usage: SeqSorter.py [-h] [-n1 [NAME1]] [-n2 [NAME2]] [-o [OUTPUT_DIR]]
                    [-p [PROCESSES]] [-e [ESTIMATE_ERR]]
                    [-c [COMPARE_RESULTS]] [-gp [GENOME1_PRIOR]]
                    [-pc [POSTERIOR_CUTOFF]]
                    [samfile1] [samfile2]

positional arguments:
  samfile1              path to samfile 1
  samfile2              path to samfile 2

optional arguments:
  -h, --help            show this help message and exit
  -n1 [NAME1], --name1 [NAME1]
                        name for genome 1 (reference for samfile 1)
  -n2 [NAME2], --name2 [NAME2]
                        name for genome 2 (reference for samfile 2)
  -o [OUTPUT_DIR], --output_dir [OUTPUT_DIR]
                        directory to write output to
  -p [PROCESSES], --processes [PROCESSES]
                        number of processes to use for sorting step, default =
                        number of CPU cores available
  -e [ESTIMATE_ERR], --estimate_err [ESTIMATE_ERR]
                        set this flag to calculate actual random mismatch
                        probabilities for more accurate mapping. WARNING: very
                        slow
  -c [COMPARE_RESULTS], --compare_results [COMPARE_RESULTS]
                        compare the results of SeqSorter w/ quality scores and
                        SeqSorter w/ estimated error rates. WARNING: very slow
  -gp [GENOME1_PRIOR], --genome1_prior [GENOME1_PRIOR]
                        prior probability that a read belongs to genome1
  -pc [POSTERIOR_CUTOFF], --posterior_cutoff [POSTERIOR_CUTOFF]
                        lower-bound cutoff for probability that a read belongs
                        to a genome for it to be classified as that genome



EXAMPLE:

python2 SeqSorter.py sample_data/S288C_bowtie_test/BY_bowtie_out.sam \
    sample_data/RM_bowtie_test/RM_bowtie_out.sam -n1 BY -n2 RM \
    -o my_output_dir

'''

import argparse
import collections
from compatibility import compatibility_dict
from compatibility import izip
from copy import copy
import estimateErrorFreq
import itertools
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import pysam
import re
import sys
import time

# parse SeqSorter program input.
def parseargs():
    
    parser = argparse.ArgumentParser(description=desc)
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('samfile1', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)
    parser.add_argument('samfile2', nargs='?', type = str, help='path to samfile 2', default=sys.stdin)
    # optional nicknames for the genomes used in the two samfiles (recommended)    
    parser.add_argument('-n1', '--name1', nargs='?', type = str, help='name for genome 1 (reference for samfile 1)', default='genome1')
    parser.add_argument('-n2', '--name2', nargs='?', type = str, help='name for genome 2 (reference for samfile 2)', default='genome2')
    # various other parameters
    parser.add_argument('-o', '--output_dir', nargs='?', type = str, help='directory to write output to', default='output/')
    parser.add_argument('-p', '--processes', nargs='?', type = str, help='number of processes to use for sorting step, default = number of CPU cores available', default=mp.cpu_count())
    parser.add_argument('-e', '--estimate_err', nargs='?', type = int, help='set this flag to calculate actual random mismatch probabilities for more accurate mapping. WARNING: very slow', const=1)
    parser.add_argument('-c', '--compare_results', nargs='?', type = int, help='compare the results of SeqSorter w/ quality scores and SeqSorter w/ estimated error rates. WARNING: very slow', const=1)
    parser.add_argument('-gp', '--genome1_prior', nargs='?', type = float, help='prior probability that a read belongs to genome1', default=0.5)
    parser.add_argument('-pc', '--posterior_cutoff', nargs='?', type = float, help='lower-bound cutoff for probability that a read belongs to a genome for it to be classified as that genome', default=0.9)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# MULTIMAPPED READ SORTER
# This is a sorter class that sorts an RNAseq read to one parental allele or the other
# Each child process processes reads using a separate instance of this class,
# ensuring complete independence of data until recombination at the end
class multimapped_read_sorter():
    
    # The following methods compute the probabilities that matched and mismatched bases contribute.
    # Let M be the probability of a miscalled base, computed as:
    #    10^((number specified by phred char)/-10)
    #    NOTE: the ascii phred char is offset by 33 for sanger format reads
    # A matched base contributes (1 - M): the probability that the base call was correct
    # and the genome in question generated the observed base.
    # A mismatched base contributes M / 3: the probability that the base call was wrong
    # and the called base was the observed 1 of 3 possible other bases.
  
    # lists that hold precomputed matched and missed base probabilities
    # index into the list using the dec value of the ascii char, which is simply
    # the value in the qual string once it is represented as a byte array
     
    log10_matched_base_prob = [0]*127
    log10_mismatched_base_prob = [0]*127
    # hold the result of 
    # Store the results of sorting in a string of chars    
    
    # The init method fills the log10 match/mismatch probability lists
    def __init__(self, samfile1=None, samfile2=None, genome1_name='genome1', genome2_name='genome2', mismatch_prob_dict=None, genome1_prior=0.5, posterior_cutoff=0.9, interleave_ix=0, output_dir='output/'):
        
        # initialize variables        
        self.samfile1 = samfile1
        self.samfile2 = samfile2
        self.genome1_name = genome1_name
        self.genome2_name = genome2_name
        self.genome1_prior = genome1_prior
        self.posterior_cutoff = posterior_cutoff
        self.category_counter = compatibility_dict(int)
        self.sort_fates = []
        self.logs = []
        self.interleave_ix = 0
        self.output_dir = output_dir
        # set up regexes
        # regex for the MD string that specifies errors from the reference.
        # more information about the MD string: page 7 of http://samtools.github.io/hts-specs/SAMv1.pdf
        self.MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)")
        # regex for MD deletions, for tracking unequal deletions
        self.DEL_REGEX = re.compile("\^[A-Z]+")

        # hardcode values for 33 since log10(0) gives an error
        self.log10_matched_base_prob[33] = -float('Inf')
        self.log10_mismatched_base_prob[33] = - 0.47712125471966244 # log10(3) == 0.47712125471966244
        
        # phred scores can be ascii 33-126
        for i in range(34,127):
            self.log10_matched_base_prob[i] = math.log10(1.0 - math.pow(10, (i - 33) * -0.1))
            self.log10_mismatched_base_prob[i] = ((i - 33) * -0.1) - 0.47712125471966244 # log10(3) == 0.47712125471966244

        # if mismatch probabilities have been computed empirically, use those
        # to overwrite the defaults
        if mismatch_prob_dict:
            with open(os.path.join(self.output_dir, 'mismatch_probs_info.txt'), 'w') as outputfile:
                print('mismatch_prob_dict', file=outputfile)
                print(mismatch_prob_dict, file=outputfile)
                print('\n',file=outputfile)
                # overwrite the old, approximate values with empirically determined ones
                for i in range(33,127):
                    if i in mismatch_prob_dict:
                        if mismatch_prob_dict[i] != 1:
                            self.log10_matched_base_prob[i] = math.log10(1.0 - mismatch_prob_dict[i])
                        self.log10_mismatched_base_prob[i] = math.log10(mismatch_prob_dict[i] / 3)

                if interleave_ix == 0:
                    
                    measured_phred = [-10*math.log10(math.pow(10,i)*3) for i in self.log10_mismatched_base_prob[33:83]]
                    qual_score_phred = range(0, 50)
                    
                    print('measured_phred', file=outputfile)
                    print(measured_phred, file=outputfile)
                    print('\n',file=outputfile)                    
                    
                    print('qual_score_phred', file=outputfile)
                    print(qual_score_phred, file=outputfile)
                    print('\n',file=outputfile) 

                    try:                    
                        fig, ax = plt.subplots()                    
                        ax.scatter(qual_score_phred, measured_phred)
                        # add an x=y line
                        ax.plot(qual_score_phred, qual_score_phred, 'r')
                        ax.set_xlabel('Quality Score', fontsize=20)
                        ax.set_ylabel('-10*log10(Prob of Mismatch)', fontsize=20)
                        ax.set_title('Quality Score vs. Calculated Probability of Mismatch')
                        ax.grid(True)
                        ax.set_xlim(-1, 51)
                        ax.set_ylim(-1, 51)
                    
                        fig.tight_layout()
        
                        plt.savefig(os.path.join(self.output_dir, 'qual_scores_vs_actual_mismatch.png'), format='png')
                    except:
                        pass            

    # LOG
    # save strings that will later be written to the verbose output file        
    def log(self, err1, err2, prob1, category, sort_fate):
        self.sort_fates.append(sort_fate)
        msg = '{}\t{}\t{}\t{}'.format(err1, err2, prob1, category)          
        self.logs.append(msg)
    
    # ALIGNED READ PROB
    # Computes the probability of a whole RNAseq read being generated by a genome, 
    # given a pysam aligned_read object of the alignment
    def aligned_read_prob(self, aligned):
        
        # create a list of (num_matched_bases, error) tuples
        # see samtools documentation for MD string
        err = re.findall(self.MD_REGEX, aligned.opt("MD"))
        
        total = 0
        seq_ix = 0
        qual = bytearray(aligned.qual)
        
        # step through sequence
        for matched_bases, curr_err in err:
            
            # step through the matched bases and sum log10 probabilities that the base was called correctly
            for match in range(0, int(matched_bases)):
                total += self.log10_matched_base_prob[qual[seq_ix]]
                seq_ix += 1

            # if there is a deletion, skip forward to after the deletion
            if '^' in curr_err:
                #this is where we should handle deletions
                pass
                
            # step through mismatched bases and sum log10 probabilities of that mismatch occuring
            else:
                total += self.log10_mismatched_base_prob[qual[seq_ix]]                
                seq_ix += 1

        # step through the remaining matched bases
        while seq_ix < len(qual):
            total += self.log10_matched_base_prob[qual[seq_ix]]
            seq_ix += 1
            
        # sum the probabilities and exponentiate to convert back from log10 scale
        return pow(10, total)

    # given two aligned read objects mapping the same read to different genomes,
    # return the name of the most likely genome and probability of genome1
    def untangle_two_mappings(self, aligned1, aligned2):
        
        genome2_prior = 1.0 - self.genome1_prior
        
        # probability of the read given that genome1 generated it
        prob_read_genome1 = self.aligned_read_prob(aligned1)

        # probabiltiy of the read given that genome2 generated it
        prob_read_genome2 = self.aligned_read_prob(aligned2)

        # apply baiyes rule: compute probability that each genome generated
        # the read given our priors for genome1 and genome2
        prob_genome1 = (prob_read_genome1 * self.genome1_prior /
                        (prob_read_genome1 * self.genome1_prior + prob_read_genome2 * genome2_prior))
                        
        prob_genome2 = 1.0 - prob_genome1    

        if (prob_genome1 >= self.posterior_cutoff):
            return 'classified1', prob_genome1, 1
        elif (prob_genome2 >= self.posterior_cutoff):
            return 'classified2', prob_genome1, 2
        else:
            return 'unclassified', prob_genome1, 0

    # UNTANGLE TWO SAMFILES
    # Given two samfile objects mapping the same RNAseq reads to different genomes,
    # sort each alignedread object to one genome or the other.
    # To distribute the work across multiple cores, each process has a unique
    # interleave_ix, and only processes the alignedread pairs where the index modulo
    # the num_processes is the interleave_ix
    def untangle_two_samfiles(self, interleave_ix=0, num_processes=1):
        
        sam1 = pysam.Samfile(self.samfile1)
        sam2 = pysam.Samfile(self.samfile2)
        
        ix = 0
        for aligned1, aligned2 in izip(sam1, sam2):
            
            if (ix % num_processes != interleave_ix):
                ix += 1
                continue
            ix += 1

            assert aligned1.qname == aligned2.qname
            
            if aligned1.is_unmapped and aligned2.is_unmapped:
                # the read does not map to either genome
                # this read is probably junk
                self.category_counter['no_match'] += 1
                self.log('NA','NA','NA', 'unmapped', 0)
                    
            elif (not aligned1.is_unmapped) and aligned2.is_unmapped:
                # the read maps to alignment1 but not alignment2; either junk or
                # an unshared gene bw BY and RM
                self.category_counter['match1'] += 1
                num_err1 = len(re.findall(self.MD_REGEX, aligned1.opt("MD")))
                self.log(num_err1, 'NA', '1', 'mapped {}'.format(self.genome1_name), 1)
            
            elif aligned1.is_unmapped and (not aligned2.is_unmapped):
                # the read maps to alignment2 but not alignment1; either junk or
                # an unshared gene bw BY and RM
                self.category_counter['match2'] += 1
                num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                self.log('NA', num_err2, '0', 'mapped {}'.format(self.genome2_name), 2)
            
            elif aligned1.opt("MD") == aligned2.opt("MD"):
                # the read has the same errors to both genomes, so it is impossible
                # to sort it one way or the other
                self.category_counter['same_errors'] += 1
                num_err1 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                self.log(num_err1, num_err2, '0.5', 'unclassified: same errors', 0)
            else:
            
                classification, prob_genome1, sort_fate = self.untangle_two_mappings(aligned1, aligned2)
                
                num_err1 = len(re.findall(self.MD_REGEX, aligned1.opt("MD")))
                num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))                
                    
                self.log(num_err1, num_err2, prob_genome1, classification, sort_fate)
                self.category_counter[classification] += 1
    
    def print_sorted_samfiles(self):
        
        # open up old, unsorted samfiles
        sam1 = pysam.Samfile(self.samfile1)
        sam2 = pysam.Samfile(self.samfile2)
        
        new_sam1_path = os.path.join(self.output_dir, self.genome1_name+'_sorted.sam')
        new_sam2_path = os.path.join(self.output_dir, self.genome2_name+'_sorted.sam')
        
        # open up samfiles to print sorted aligned_reads to
        new_sam1 = pysam.Samfile(new_sam1_path, 'wh', template=sam1)
        new_sam2 = pysam.Samfile(new_sam2_path, 'wh', template=sam2)
        
        for aligned1, aligned2, sort_fate in izip(sam1, sam2, self.sort_fates):
            
            if sort_fate == 1:
                # the current RNAseq read was classified as genome1
                new_sam1.write(aligned1)
            elif sort_fate == 2:
                # the current RNAseq read was classified as genome2
                new_sam2.write(aligned2)
            elif sort_fate == 0:
                # the current RNAseq read was unclassified due to same errors or under cutoff
                pass
            else:
                # This should never ever happen
                assert(False)

# WORKER PROCEDURE
# The procedure called by different processes using apply_async.
# Processes aligned reads, moving in skips of size interleave_ix. It is recommended
# that you do not call this from an external module
def _worker_procedure(samfile1, samfile2, interleave_ix, num_processes, mismatch_prob_dict, genome1_name, genome2_name, genome1_prior, posterior_cutoff, output_dir):

    sorter = multimapped_read_sorter(samfile1, samfile2, genome1_name, genome2_name, mismatch_prob_dict, genome1_prior, posterior_cutoff, interleave_ix, output_dir)
    sorter.untangle_two_samfiles(interleave_ix, num_processes)
    
    return sorter

# MERGE SORTERS
# Takes a list of multimapped_read_sorter objects and combines their results
# into a single sorter. Used to combine results computed by different processes
# into one.
def merge_sorters(sorter_list):
    new_sorter = sorter_list.pop(0)
    
    # the 'iters' solution allows merging logs in order of original reads, credit to Mark Byers:     
    # http://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
    log_iters = [iter(copy(new_sorter.logs))]
    sort_fate_iters = [iter(copy(new_sorter.sort_fates))]
    
    # merge info
    while(sorter_list):
        old_sorter = sorter_list.pop(0)
        log_iters.append(iter(old_sorter.logs))
        sort_fate_iters.append(iter(old_sorter.sort_fates))
        
        for category, count in old_sorter.category_counter.items():
            new_sorter.category_counter[category] += count
    
    # create the combined logs and sort fate list by cycling through the individual
    # ones and repeatedly taking the first thing off
    new_sorter.logs = list(next(it) for it in itertools.cycle(log_iters))
    new_sorter.sort_fates = list(next(it) for it in itertools.cycle(sort_fate_iters))    
    return new_sorter

# MAIN FUNCTION
# Create a pool of worker processes and have them process aligned_reads from samfile1
# and samfile2 in an interleaved fashion 
def main():
    
    # Start timer
    t1 = time.time()
    
    # Get command line args
    args = parseargs()
    samfile1 = args.samfile1
    samfile2 = args.samfile2
    genome1_name = args.name1
    genome2_name = args.name2
    num_processes = args.processes
    compare_results = args.compare_results
    estimate_error_prob = True if compare_results else args.estimate_err
    genome1_prior = args.genome1_prior
    posterior_cutoff = args.posterior_cutoff
    output_dir = args.output_dir
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
                
    # If estimate_error_prob is True, then calculate actual probabilities of
    # random mismatch for each phred score.
    if estimate_error_prob:
        mismatch_prob_dict, total_values_dict = estimateErrorFreq.create_mismatch_prob_dict(samfile1, samfile2, genome1_name, genome2_name)
        with open(os.path.join(output_dir, 'mismatch_probs_info.txt'), 'w') as outputfile:
            print('Total_values_dict', file=outputfile)
            print(total_values_dict, file=outputfile)
    else:
        mismatch_prob_dict = None

    # Create a pool of worker processes to do the sorting
    worker_pool = mp.Pool(processes=num_processes)
    async_results = []
    
    # Create a set of interleave indices, to allow each process to only work on
    # one alignment every x alignments, where x is the number of processes running    
    for interleave_ix in range(0, num_processes):
        args = (samfile1, samfile2, interleave_ix, num_processes, mismatch_prob_dict, genome1_name, genome2_name, genome1_prior, posterior_cutoff, output_dir)
        async_results.append(worker_pool.apply_async(_worker_procedure, args))
    
    # Join the processes back to the main thread    
    worker_pool.close()
    worker_pool.join()
    
    # Wait and retrieve the results from each child process
    unpacked_results = []
    for result in async_results:
        result.wait()
        unpacked_results.append(result.get())
    
    # Combine the results from each child process into the parent process    
    combined_sorter = merge_sorters(unpacked_results)
    
    # redo the analysis without calculating mismatch probabilities empirically,
    # and compare the results
    if compare_results:
        
        # Create a pool of worker processes to do the sorting
        worker_pool = mp.Pool(processes=num_processes)
        async_results = []    
        
        # Repeat analysis without mismatch_prob_dict  
        for interleave_ix in range(0, num_processes):
            args = (samfile1, samfile2, interleave_ix, num_processes, None, genome1_name, genome2_name, genome1_prior, posterior_cutoff, output_dir)
            async_results.append(worker_pool.apply_async(_worker_procedure, args))
        
        # Join the processes back to the main thread    
        worker_pool.close()
        worker_pool.join()
        
        # Wait and retrieve the results from each child process
        unpacked_results = []
        for result in async_results:
            result.wait()
            unpacked_results.append(result.get())
        
        # Combine the results from each child process into the parent process    
        combined_sorter_qual_scores = merge_sorters(unpacked_results)
        
        # compare the results of calculated error probs vs. quality score derived error probs
        compare_counter = collections.defaultdict(int)        
        for calc_result, qual_result in zip(combined_sorter.sort_fates, combined_sorter_qual_scores.sort_fates):
            if (calc_result != qual_result):
                compare_counter['different'] += 1
            else:
                compare_counter['same'] += 1
        
        print('Comparison of quality score error probs to calculated error probs:')
        print('{} reads : same sorting decision'.format(compare_counter['same']))
        print('{} reads : different sorting decision'.format(compare_counter['different']))

    category_counter = combined_sorter.category_counter
    
    labels =([
            'unmapped by bowtie',
            'mapped by bowtie uniquely to {}'.format(genome1_name),
            'mapped by bowtie uniquely to {}'.format(genome2_name),
            'double-mapped, but errors are same'.format(genome1_name, genome2_name),
            'double-mapped, classified as {} based on errors'.format(genome1_name),
            'double-mapped, classified as {} based on errors'.format(genome2_name),
            'double-mapped, unclassifiable based on errors',
            ])
    fracs = ([
            category_counter['no_match'],
            category_counter['match1'],
            category_counter['match2'],
            category_counter['same_errors'],
            category_counter['classified1'],
            category_counter['classified2'],
            category_counter['unclassified'],
            ])
    
    for num, label in zip(fracs, labels):
        print('{}\t{}'.format(num, label))
    
    # Print the total time
    t2 = time.time()
    print('TOTAL TIME: {}'.format(t2-t1))
    
    # print to the new sorted samfiles
    combined_sorter.print_sorted_samfiles()
    
    verbose_filepath = os.path.join(output_dir, 'supplementary_output.txt')
    # print all the logs to the verbose output file
    verbose_file = open(verbose_filepath, 'w')
    print('err {}\terr {}\tprob {}\tcategory'.format(genome1_name, genome2_name, genome1_name), file=verbose_file)
    for line in combined_sorter.logs:
        print(line,file=verbose_file)
    verbose_file.close()

if __name__ == '__main__':
    main()
