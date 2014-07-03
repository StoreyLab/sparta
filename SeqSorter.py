# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
sys.path.insert(0,'fake_root/usr/local/python/lib/python2.7/site-packages/')

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

usage: SeqSorter.py [-h] [-n1 [NAME1]] [-n2 [NAME2]] [-v [VERBOSE_FILE]]
                    [-p [PROCESSES]] [-e [ESTIMATE_ERR]] [-gp [GENOME1_PRIOR]]
                    [-pc [POSTERIOR_CUTOFF]] [-po [PLOT_OUTPUT]]
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
  -v [VERBOSE_FILE], --verbose_file [VERBOSE_FILE]
                        specify a filename for verbose output (useful for
                        statistical analysis)
  -p [PROCESSES], --processes [PROCESSES]
                        number of processes to use for sorting step, default =
                        number of CPU cores available
  -e [ESTIMATE_ERR], --estimate_err [ESTIMATE_ERR]
                        set this flag to calculate actual random mismatch
                        probabilities for more accurate mapping. WARNING: very
                        slow
  -gp [GENOME1_PRIOR], --genome1_prior [GENOME1_PRIOR]
                        prior probability that a read belongs to genome1
  -pc [POSTERIOR_CUTOFF], --posterior_cutoff [POSTERIOR_CUTOFF]
                        lower-bound cutoff for probability that a read belongs
                        to a genome for it to be classified as that genome
  -po [PLOT_OUTPUT], --plot_output [PLOT_OUTPUT]
                        filepath to output plot of qual score vs probability
                        of mismatch


EXAMPLE:

python2 SeqSorter.py sample_data/S288C_bowtie_test/BY_bowtie_out.sam \
    sample_data/RM_bowtie_test/RM_bowtie_out.sam -n1 BY -n2 RM \
    -v verbose_outfile

'''

import argparse
import collections
from copy import copy as copy
import estimateErrorFreq
import itertools
import math
from matplotlib import pyplot as plt
import multiprocessing as mp
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
    # this is an optional verbose outfile with useful values for statistical analysis
    parser.add_argument('-v', '--verbose_file', nargs='?', type=argparse.FileType('w'), help='specify a filename for verbose output (useful for statistical analysis)', default=None)    
    # various other parameters
    parser.add_argument('-p', '--processes', nargs='?', type = str, help='number of processes to use for sorting step, default = number of CPU cores available', default=mp.cpu_count())
    parser.add_argument('-e', '--estimate_err', nargs='?', type = int, help='set this flag to calculate actual random mismatch probabilities for more accurate mapping. WARNING: very slow', const=1)
    parser.add_argument('-gp', '--genome1_prior', nargs='?', type = float, help='prior probability that a read belongs to genome1', default=0.5)
    parser.add_argument('-pc', '--posterior_cutoff', nargs='?', type = float, help='lower-bound cutoff for probability that a read belongs to a genome for it to be classified as that genome', default=0.9)
    parser.add_argument('-po', '--plot_output', nargs='?', type = str, help='filepath to output plot of qual score vs probability of mismatch', default='qual_scores_vs_actual_mismatch.png')

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
    
    # The init method fills the log10 match/mismatch probability lists
    def __init__(self, mismatch_prob_dict=None, verbose=False, interleave_ix=0, plot_output='qual_scores_vs_actual_mismatch.png'):
        
        # if verbose file was specified then run in verbose mode

        self.verbose = verbose
        self.category_counter = collections.defaultdict(int)
        self.logs = []
        self.interleave_ix = 0
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
            
            # overwrite the old, approximate values with empirically determined ones
            for i in range(33,127):
                if i in mismatch_prob_dict:
                    if mismatch_prob_dict[i] != 1:
                        self.log10_matched_base_prob[i] = math.log10(1.0 - mismatch_prob_dict[i])
                    self.log10_mismatched_base_prob[i] = math.log10(mismatch_prob_dict[i] / 3)

            for i in range(33,127):
                print ('{}\t{}\t{}'.format(i,self.log10_matched_base_prob[i],self.log10_mismatched_base_prob[i]))
            
            if interleave_ix == 0:
                fig, ax = plt.subplots()
                
                expected_phred = [-10*math.log10(math.pow(10,i)*3) for i in self.log10_mismatched_base_prob[33:83]]
                actual_phred = range(1, 51)
                
                ax.scatter(actual_phred, expected_phred)
                # add an x=y line
                ax.plot(actual_phred, actual_phred, 'r')
                ax.set_xlabel('Quality Score', fontsize=20)
                ax.set_ylabel('-10*log10(Prob of Mismatch)', fontsize=20)
                ax.set_title('Quality Score vs. Calculated Probability of Mismatch')
                ax.grid(True)
                ax.set_xlim(-1, 51)
                ax.set_ylim(-1, 51)
                
                fig.tight_layout()
    
                plt.savefig(plot_output, format='png')
            
    # LOG
    # save strings that will later be written to the verbose output file        
    def log(self, err1, err2, prob1, category):
        if self.verbose:
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

    # UNTANGLE TWO SAMFILES
    # Given two samfile objects mapping the same RNAseq reads to different genomes,
    # sort each alignedread object to one genome or the other.
    # To distribute the work across multiple cores, each process has a unique
    # interleave_ix, and only processes the alignedread pairs where the index modulo
    # the num_processes is the interleave_ix
    def untangle_two_samfiles(self, sam1, sam2, interleave_ix=0, num_processes=1, genome1_name='genome1', genome2_name='genome2', genome1_prior=0.5, posterior_cutoff=0.9, verbose=False):
        
        ix = 0
        for aligned1, aligned2 in itertools.izip(sam1, sam2):
            
            if (ix % num_processes != interleave_ix):
                ix += 1
                continue
            ix += 1

            assert aligned1.qname == aligned2.qname
            
            if aligned1.is_unmapped and aligned2.is_unmapped:
                # the read does not map to either genome
                # this read is probably junk
                self.category_counter['no_match'] += 1
                if verbose:
                    self.log('NA','NA','NA', 'unmapped')
                    
            elif (not aligned1.is_unmapped) and aligned2.is_unmapped:
                # the read maps to alignment1 but not alignment2; either junk or
                # an unshared gene bw BY and RM
                self.category_counter['match1'] += 1
                if verbose:
                    num_err1 = len(re.findall(self.MD_REGEX, aligned1.opt("MD")))
                    self.log(num_err1, 'NA', '1', 'mapped {}'.format(genome1_name))
            
            elif aligned1.is_unmapped and (not aligned2.is_unmapped):
                # the read maps to alignment2 but not alignment1; either junk or
                # an unshared gene bw BY and RM
                self.category_counter['match2'] += 1
                if verbose:
                    num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                    self.log('NA', num_err2, '0', 'mapped {}'.format(genome2_name))
            
            elif aligned1.opt("MD") == aligned2.opt("MD"):
                # the read has the same errors to both genomes, so it is impossible
                # to sort it one way or the other
                self.category_counter['same_errors'] += 1
                if verbose:
                    num_err1 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                    num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                    self.log(num_err1, num_err2, '0.5', 'unclassified: same errors')
            else:
                # bowtie matched the read to both genomes
                # use errors and qual scores to classify read one way or the other
                genome2_prior = 1.0 - genome1_prior
                
                # probability of the read given that genome1 generated it
                prob_read_genome1 = self.aligned_read_prob(aligned1)
            
                # probabiltiy of the read given that genome2 generated it
                prob_read_genome2 = self.aligned_read_prob(aligned2)
            
                # apply baiyes rule: compute probability that each genome generated
                # the read given our priors for genome1 and genome2
                prob_genome1 = (prob_read_genome1 * genome1_prior /
                                (prob_read_genome1 * genome1_prior + prob_read_genome2 * genome2_prior))
                                
                prob_genome2 = 1.0 - prob_genome1    
                num_err1 = len(re.findall(self.MD_REGEX, aligned1.opt("MD")))
                num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))                
                    
                if (prob_genome1 >= posterior_cutoff):
                    self.category_counter['classified1'] += 1
                    self.log(num_err1, num_err2, prob_genome1, 'classified {}'.format(genome1_name))
                elif (prob_genome2 >= posterior_cutoff):
                    self.category_counter['classified2'] += 1
                    self.log(num_err1, num_err2, prob_genome1, 'classified {}'.format(genome2_name))
                else:
                    self.category_counter['unclassified'] += 1
                    self.log(num_err1, num_err2, prob_genome1, 'unclassified: under cutoff')

# WORKER PROCEDURE
# The procedure called by different processes using apply_async.
# Processes aligned reads, moving in skips of size interleave_ix. It is recommended
# that you do not call this from an external module
def _worker_procedure(samfile1, samfile2, interleave_ix, num_processes, mismatch_prob_dict, genome1_name, genome2_name, genome1_prior, posterior_cutoff, verbose, plot_output):

    sorter = multimapped_read_sorter(mismatch_prob_dict, verbose, interleave_ix, plot_output)
    sam1 = pysam.Samfile(samfile1)
    sam2 = pysam.Samfile(samfile2)
    sorter.untangle_two_samfiles(sam1, sam2, interleave_ix, num_processes, genome1_name, genome2_name, genome1_prior, posterior_cutoff, verbose)
    return sorter

# MERGE SORTERS
# Takes a list of multimapped_read_sorter objects and combines their results
# into a single sorter. Used to combine results computed by different processes
# into one.
def merge_sorters(sorter_list):
    new_sorter = sorter_list.pop(0)
    
    # the 'iters' solution allows merging logs in order of original reads, credit to Mark Byers:     
    # http://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
    iters = [iter(new_sorter.logs)]
    
    # merge info
    while(sorter_list):
        old_sorter = sorter_list.pop(0)
        iters.append(iter(old_sorter.logs))
        
        for category, count in old_sorter.category_counter.iteritems():
            new_sorter.category_counter[category] += count
    
    new_sorter.logs = list(it.next() for it in itertools.cycle(iters))
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
    verbose_file = args.verbose_file
    verbose = True if verbose_file else False    
    num_processes = args.processes
    estimate_error_prob = args.estimate_err
    genome1_prior = args.genome1_prior
    posterior_cutoff = args.posterior_cutoff
    plot_output = args.plot_output

    # If estimate_error_prob is True, then calculate actual probabilities of
    # random mismatch for each phred score.
    if estimate_error_prob:
        mismatch_prob_dict = estimateErrorFreq.create_mismatch_prob_dict(samfile1, samfile2, genome1_name, genome2_name)
    else:
        mismatch_prob_dict = None    
    
    # Create a pool of worker processes to do the sorting
    worker_pool = mp.Pool(processes=num_processes)
    async_results = []    
    
    # Create a set of interleave indices, to allow each process to only work on
    # one alignment every x alignments, where x is the number of processes running    
    for interleave_ix in range(0, num_processes):
        args = (samfile1, samfile2, interleave_ix, num_processes, mismatch_prob_dict, genome1_name, genome2_name, genome1_prior, posterior_cutoff, verbose, plot_output)
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
    total_result = merge_sorters(unpacked_results)

    category_counter = total_result.category_counter
    
    # Print the total count of each type of alignment pair    
    print('\n{}\tunmapped by bowtie to either {} or {}'.format(category_counter['no_match'], genome1_name, genome2_name))
    print('{}\tmapped by bowtie to either {} or {}'.format(category_counter['match1'], genome1_name, genome2_name))
    print('{}\tmapped by bowtie to both {} or {}, but errors are same'.format(category_counter['same_errors'], genome1_name, genome2_name))
    total_mapped_to_both = category_counter['classified1'] + category_counter['classified2'] + category_counter['unclassified']
    print('{}\tmapped by bowtie to both {} and {}:'.format(total_mapped_to_both, genome1_name, genome2_name))
    print('   {}\tassigned to {} based on errors'.format(category_counter['classified1'], genome1_name))
    print('   {}\tassigned to {} based on errors'.format(category_counter['classified2'], genome2_name))
    print('   {}\tunmapped based on errors'.format(category_counter['unclassified']))
    
    # Print the total time
    t2 = time.time()
    print('TOTAL TIME: {}'.format(t2-t1))
    
    # print all the logs to the verbose output file
    if verbose_file:
        print('err {}\terr {}\tprob {}\tcategory'.format(genome1_name, genome2_name, genome1_name), file=verbose_file)
        for line in total_result.logs:
            print(line,file=verbose_file)
    
    # close output file
    verbose_file.close()

if __name__ == '__main__':
    main()