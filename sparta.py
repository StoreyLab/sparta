#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
'''
Created on Mon Jun 16 10:56:22 2014

@author: Peter Edge
'''
desc = '''
SPARTA takes a set of SAM format files that each map the same RNA reads to a different
ancestral (or parental) genome. This program classifies each read to one of the
ancestral alleles or deems it unclassifiable, based on the assumption that each
read belongs to one of the ancestral allele types provided.
'''
'''
COMMAND LINE USAGE:
python2 sparta.py <samfile1> <samfile2> <optional arguments>

Run the program without arguments to see the description of optional arguments

'''

#imports
import argparse
from compatibility import compatibility_dict
from compatibility import izip
from copy import copy
import calculate_mismatch_probs
import itertools
import math
import multiprocessing as mp
import os
import pysam
import re
import sys
import time
from util import fix_read_mate_order

# parse sparta program input.
def parseargs():

    parser = argparse.ArgumentParser(description=desc)
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('samfiles', nargs='+', type = str, help='input samfiles', default=[])
    parser.add_argument('-pe', '--paired_end', nargs='?', type = int, help='set this flag to specify that reads are paired end (default: False)', const=1)
    # optional nicknames for the genomes that the samfiles map to (recommended)    
    parser.add_argument('-n', '--names', nargs='+', type = str, help='list of nicknames for genomes corresponding for samfile1,samfile2, etc.', default=[])
    # various other parameters
    parser.add_argument('-o', '--output_dir', nargs='?', type = str, help='directory to write output to', default='output')
    parser.add_argument('-ss', '--separated_samfiles', nargs='+', type = str, help='list of filenames to write separated (classified) sam outputs. default: outputdir/genome1_sorted.sam...', default=[])
    parser.add_argument('-p', '--processes', nargs='?', type = int, help='number of processes to use for separation step, default = number of CPU cores available', default=mp.cpu_count())
    parser.add_argument('-c', '--calculate_mismatches', nargs='?', type = int, help='set this flag to calculate actual mismatch probabilities for more accurate mapping. WARNING: very slow', const=1)
    parser.add_argument('-mp', '--mismatch_prob_inputfile', nargs='?', type = int, help='specify an existing file with mismatch probabilities per quality score for more accurate mapping.', default = None)
    parser.add_argument('-ph', '--pileup_height', nargs='?', type = int, help='if calculate_mismatches is True, specify minimum height of read pileup to consider, default = 20', default=20)
    parser.add_argument('-se', '--sample_every', nargs='?', type = int, help='if calculate_mismatches is True, specify N such that calculate_mismatch_probs only samples every N reads, default = 10', default=10)
    parser.add_argument('-gp', '--genome_priors', nargs='+', type = float, help='list of prior probabilities that a read belongs to each genome', default=[])
    parser.add_argument('-pc', '--posterior_cutoff', nargs='?', type = float, help='lower-bound cutoff for probability that a read belongs to a genome for it to be classified as that genome', default=0.99)
    parser.add_argument('-q', '--quiet', nargs='?', type = int, help='set this flag to supress writing final counts to stdout (default: False)', const=1)
    
    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# MULTIMAPPED READ SEPARATOR
# This is a separator class that sorts an RNAseq read to one parental allele or the other
# Each child process processes reads using a separate instance of this class,
# ensuring complete independence of data until recombination at the end
class multimapped_read_separator():
    

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
    # Store the results of separation in a string of chars    
    
    # The init method fills the log10 match/mismatch probability lists
    def __init__(self, samfiles=None, paired_end=False, mismatch_prob_dict=None,
                 transition_prob_dict=None, genome_priors=None, posterior_cutoff=0.99, interleave_ix=0):
        
        # initialize variables        
        self.samfiles = samfiles
        self.paired_end = paired_end
        # assume even prior proababilties
        self.genome_priors = genome_priors if genome_priors else [1.0/len(samfiles)]*len(samfiles)
        self.posterior_cutoff = posterior_cutoff
        self.category_counter = compatibility_dict(int)
        self.sort_fates = []
        self.mismatch_logs = []
        self.posterior_logs = []
        self.interleave_ix = 0
        # set up regexes
        # regex for the MD string that specifies errors from the reference.
        # more information about the MD string: page 7 of http://samtools.github.io/hts-specs/SAMv1.pdf
        self.MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)")
        # regex for MD deletions, for tracking unequal deletions
        self.DEL_REGEX = re.compile("\^[A-Z]+")
        # regex for N, for tracking ambiguous bases
        self.N_REGEX = re.compile('\d[nN]\d')


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

        if transition_prob_dict:
            
            self.log10_transition_prob_dict = {}
            for k, v in transition_prob_dict.item():
                
                self.log10_transition_prob_dict[k] = math.log10(v)
        else:
            self.log10_transition_prob_dict = None
            

    # LOG
    # save strings that will later be written to the various output files
    def log(self, mismatch_list, posterior_list, sort_fate):
        self.sort_fates.append(sort_fate)
        self.mismatch_logs.append(mismatch_list)
        self.posterior_logs.append(posterior_list)

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
                
            elif 'N' in curr_err or 'n' in curr_err:
                # any base counts as a match to an N on the reference
                total += self.log10_matched_base_prob[qual[seq_ix]]
                seq_ix += 1

            # step through mismatched bases and sum log10 probabilities of that mismatch occuring
            else:
                if self.log10_transition_prob_dict:
                    total += self.log10_transition_prob_dict[(qual[seq_ix], curr_err, aligned.seq[seq_ix])]
                else:
                    total += self.log10_mismatched_base_prob[qual[seq_ix]]
                    
                seq_ix += 1

        # step through the remaining matched bases
        while seq_ix < len(qual):
            total += self.log10_matched_base_prob[qual[seq_ix]]
            seq_ix += 1
            
        # sum the probabilities and exponentiate to convert back from log10 scale
        return pow(10, total)

    def count_mismatches(self, aligned_read_set):

        if not aligned_read_set:
            return None
            
        mismatch_counts = []
        for aligned_read in aligned_read_set:
            if not aligned_read.is_unmapped:               
                mismatches = re.findall(self.MD_REGEX, aligned_read.opt("MD"))
                len_mismatches = len(mismatches) - len(re.findall(self.N_REGEX, aligned_read.opt("MD")))
                assert len_mismatches >= 0
                mismatch_counts.append(len_mismatches)
            else:
                mismatch_counts.append('NA')
                
        return mismatch_counts


    def baiyes_classify(self, read_probs):
        if read_probs == [0.0] * len(read_probs):
            return 'unmapped', read_probs, 0

        # all the priors where the read mapped else 0
        sum_mapped_priors = sum([prior if prob else 0.0 for prob, prior in zip(read_probs, self.genome_priors)])
        # correct the original priors for unmapped reads
        corrected_priors = [prior / sum_mapped_priors if prob else 0.0 for prob, prior in zip(read_probs, self.genome_priors)]
        assert sum(corrected_priors) == 1.0
        
        # apply baiyes rule: compute probability that each genome generated
        # the read given our priors for each genome
        baiyes_denom = sum([prob * prior for prob, prior in zip(read_probs, corrected_priors)])         
        read_posterior_probs = [prob * prior / baiyes_denom for prob, prior in zip(read_probs, corrected_priors)]
        
        for i, posterior in enumerate(read_posterior_probs):
            if posterior > self.posterior_cutoff:
                # add 1 because sort fate of 0 means unsorted
                sort_fate = i + 1
                return 'classified{}'.format(sort_fate), read_posterior_probs, sort_fate
        
        return 'unclassified', [0.0]*len(read_probs) , 0            


    # given a list of aligned read objects mapping the same read to different genomes,
    # return the name of the most likely genome and probability of genome1
    def untangle_mappings(self, aligned_reads, aligned_read_mates=None):
    
        # probability of the read given that genome1 generated it
        if not aligned_read_mates: # single reads
            read_probs = [self.aligned_read_prob(read) if not read.is_unmapped else 0.0 for read in aligned_reads]
            classification, posterior_probs, sort_fate = self.baiyes_classify(read_probs)
            mismatch_counts = self.count_mismatches(aligned_reads)
            self.log(mismatch_counts, posterior_probs, sort_fate)
            return classification, posterior_probs, sort_fate

        else: # paired end reads
            read_probs = [(self.aligned_read_prob(read) if not read.is_unmapped else 0.0) for read in aligned_reads]
            mate_probs = [(self.aligned_read_prob(read_mate) if not read_mate.is_unmapped else 0.0) for read_mate in aligned_read_mates]
            # check if, for any genome, read and mate mapped
            total_read_probs = [read_prob * mate_prob for read_prob, mate_prob in zip(read_probs, mate_probs)]
            mismatch_counts = self.count_mismatches(aligned_reads)
            mismatch_counts_mates = self.count_mismatches(aligned_read_mates)

            if total_read_probs != [0.0] * len(aligned_reads):
                
                # for at least one genome, both read and mate mapped
                # so we classify using the combined probability of reads + mates
                classification, posterior_probs, sort_fate = self.baiyes_classify(total_read_probs)
                
                self.log(mismatch_counts, read_probs, sort_fate)
                self.log(mismatch_counts_mates, read_probs, sort_fate)
                return classification, posterior_probs, sort_fate

            else:
                    
                # this can only occur if data is paired end and there is not
                # a single genome for which read and mate both mapped.
                # solution: we classify based on the strongest posterior from
                # either the read or the mate.
                classification1, posterior_probs1, sort_fate1 = self.baiyes_classify(read_probs)
                classification2, posterior_probs2, sort_fate2 = self.baiyes_classify(mate_probs)
                
                # which yielded a stronger posterior, using reads or their mates?
                # return this result
                max1 = max(posterior_probs1)
                max2 = max(posterior_probs2)
            
                if max(max1, max2) == max1:
                    self.log(mismatch_counts, read_probs, sort_fate1)
                    self.log(mismatch_counts_mates, read_probs, sort_fate1)
                    return classification1, posterior_probs1, sort_fate1
                else:
                    self.log(mismatch_counts, read_probs, sort_fate2)
                    self.log(mismatch_counts_mates, read_probs, sort_fate2)
                    return classification2, posterior_probs2, sort_fate2

    # UNTANGLE SAMFILES
    # Given a list of samfiles mapping the same RNAseq reads to different genomes,
    # classify each read to one of the genomes or deem it unclassifiable.
    # To distribute the work across multiple cores, each process has a unique
    # interleave_ix, and only processes the alignedread pairs where the index modulo
    # the num_processes is the interleave_ix
          
    def untangle_samfiles(self, interleave_ix=0, num_processes=1):
        
        sams_tuple = tuple(pysam.Samfile(sam_name) for sam_name in self.samfiles)
        
        ix = 0
                            
        if not self.paired_end: # default: single reads
                        
            for aligned_read_set in izip(*sams_tuple):

                if (ix % num_processes != interleave_ix):
                    ix += 1
                    continue
                ix += 1
                
                classification, posterior_probs, sort_fate = self.untangle_mappings(aligned_read_set)
                
        else: # (if paired end reads)
            
            zipped_sams = izip(*sams_tuple)
            
            for aligned_read_set in zipped_sams:
                                
                # take aligned pairs in sets of 2 to also get the mate
                # don't forget to assign aligned_pair to next tuple before end of iter
                next_tuple = next(zipped_sams)       
                aligned_read_mate_set = next_tuple
                
                if (ix % num_processes != interleave_ix):
                    aligned_read_set = next_tuple
                    ix += 1
                    continue
                ix += 1
                    
                aligned_read_set, aligned_read_mate_set = fix_read_mate_order(aligned_read_set, aligned_read_mate_set)
                
                classification, posterior_probs, sort_fate = self.untangle_mappings(aligned_read_set, aligned_read_mate_set)
                # CRUCIAL STEP: skip the mate read
                aligned_read_set = next_tuple
        
        for sam in sams_tuple:
            sam.close()
                

    def write_samfiles(self, new_sam_paths):
        
        # open up old, unsorted samfiles
        old_sams = [pysam.Samfile(sam) for sam in self.samfiles]

        # open up samfiles to print sorted aligned_reads to
        new_sams = [pysam.Samfile(new_sam_path, 'wh', template=old_sam) for old_sam, new_sam_path in zip(old_sams, new_sam_paths)]
                    
        for sort_fate, aligned_read_set in izip(self.sort_fates, izip(*old_sams)):

            if sort_fate:
                
                new_sams[sort_fate-1].write(aligned_read_set[sort_fate-1])
            
        for sam in old_sams:
            sam.close()
        
        for sam in new_sams:
            sam.close()


# WORKER PROCEDURE
# The procedure called by different processes using apply_async.
# Processes aligned reads, moving in skips of size interleave_ix. It is recommended
# that you do not call this from an external module
def _worker_procedure(samfiles, paired_end, interleave_ix, num_processes, mismatch_prob_dict, transition_prob_dict, genome_priors, posterior_cutoff):

    separator = multimapped_read_separator(samfiles, paired_end, mismatch_prob_dict, transition_prob_dict, genome_priors, posterior_cutoff, interleave_ix)
    separator.untangle_samfiles(interleave_ix, num_processes)
    
    return separator

# MERGE SEPARATORS
# Takes a list of multimapped_read_separator objects and combines their results
# into a single separator. Used to combine results computed by different processes
# into one.
def merge_separators(separator_list):
    
    new_separator = separator_list.pop(0)
    # the 'iters' solution allows merging logs in order of original reads, credit to Mark Byers:     
    # http://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
    mismatch_log_iters = [iter(copy(new_separator.mismatch_logs))]
    posterior_log_iters = [iter(copy(new_separator.posterior_logs))]
    sort_fate_iters = [iter(copy(new_separator.sort_fates))]

    # merge info
    while(separator_list):
        old_separator = separator_list.pop(0)
        mismatch_log_iters.append(iter(old_separator.mismatch_logs))
        posterior_log_iters.append(iter(old_separator.posterior_logs))
        sort_fate_iters.append(iter(old_separator.sort_fates))
        
        for category, count in old_separator.category_counter.items():
            new_separator.category_counter[category] += count
    
    # create the combined logs and sort fate list by cycling through the individual
    # ones and repeatedly taking the first thing off
    new_separator.mismatch_logs = list(next(it) for it in itertools.cycle(mismatch_log_iters))
    new_separator.posterior_logs = list(next(it) for it in itertools.cycle(posterior_log_iters))
    new_separator.sort_fates = list(next(it) for it in itertools.cycle(sort_fate_iters))
    return new_separator

# SPARTA
# This is the main program logic and is the recommended means for calling from other modules
# Takes a list of samfiles that map the same RNAseq reads (mapped via bowtie from same
# RNA-seq fastq file) to separate genomes, and sorts them to parental allele types
# Creates a pool of worker processes and has them process aligned_reads from samfile1
# and samfile2 in an interleaved fashion
def sparta(samfiles, paired_end=False, genome_names=[],
                  num_processes=mp.cpu_count(), calculate_mismatches=False,
                  pileup_height=20, sample_every=10, genome_priors=[], posterior_cutoff=0.99,
                  output_dir='output', separated_samfiles=[], quiet=False, mismatch_prob_inputfile=None):
            
    # IT IS RECOMMENDED NOT TO MOVE THIS DEFAULT ARGUMENT HANDLING CODE
    # default argument handling has to go here for genome_names, genome_priors, and separated_samfiles
    # because they depend on the length of samfiles
    if len(samfiles) < 2:
        print('ERROR: number of samfiles must be at least 2.', file=sys.stderr)
        sys.exit(-1)
    
    if genome_names == []:
        genome_names = ['genome{}'.format(i) for i in range(1, len(samfiles)+1)]
    elif len(genome_names) != len(samfiles):
        print('ERROR: length of genome name list (-n argument) should be equal to number of samfiles', file=sys.stderr)
        sys.exit(-1)
        
    if genome_priors == []:
        genome_priors = [1.0/len(samfiles)]*len(samfiles)
    elif len(genome_priors) != len(samfiles):
        print('ERROR: length of genome prior probability list (-gp argument) should be equal to number of samfiles', file=sys.stderr)
        sys.exit(-1)

    if separated_samfiles == []:
        separated_samfiles = [os.path.join(output_dir, '{}_sorted.sam'.format(name)) for name in genome_names]
    elif len(separated_samfiles) != len(samfiles):
        print('ERROR: length of separated (classified) sam output list (-ss argument) should be equal to number of samfiles', file=sys.stderr)
        sys.exit(-1)
        
    if calculate_mismatches and mismatch_prob_inputfile:
        print('\nERROR: Incompatible command-line options. Mismatch probability inputfile specified (-mp), and calculate mismatch probabilies'+
        'also specfiied.\n', file=sys.stderr)
        sys.exit(-1)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir) 
    
    # If calculate_mismatches is True, then calculate actual probabilities of
    # random mismatch for each phred score.
    if calculate_mismatches:
        mismatch_prob_dict, mismatch_prob_total_values, transition_prob_dict = calculate_mismatch_probs.create_mismatch_prob_dict(samfiles,
                                                        output_dir, paired_end, pileup_height, sample_every)

    elif mismatch_prob_inputfile:
        
        mismatch_prob_dict = {}
        mismatch_prob_total_values = {}
        with open(mismatch_prob_inputfile, 'r') as inf:
            for line in inf:
                qual, prob, total = line.rstrip.split('\t')
                mismatch_prob_dict[qual] = prob
                mismatch_prob_total_values[qual] = total
                
    else:
        mismatch_prob_dict = None
        mismatch_prob_total_values = None
        transition_prob_dict = None
    
    if paired_end == False:
        # If in single read mode, perform rudimentary test for paired end reads.
        # if first two reads have same qname attribute, print a warning that the reads
        # are likely paired end but sparta mode is single read
        paired_end_test_ix = 0
        prev_qname = None
        test_sam = pysam.Samfile(samfiles[0])
        for aligned in test_sam:
            if paired_end_test_ix >= 1:
                if prev_qname == aligned.qname:
                    print('\nWARNING: In single-read mode, but file appears to contain '+
                    'paired-end reads\n         Consider rerunning with -pe option\n', file=sys.stderr)
                break
            prev_qname = aligned.qname
            paired_end_test_ix += 1
    
        test_sam.close()

    # whether or not multiprocessed, combined_separator will hold the results        
    combined_separator = None
    
    # if multiprocessing, then spawn processes, do work, and recombine
    # else, just use the main thread to make the function call
    # this is WAY EASIER for debugging because errors and exceptions aren't just
    # thrown away by a rogue child thread
    if num_processes > 1:
        # Create a pool of worker processes to do the classifications
        worker_pool = mp.Pool(processes=num_processes)
        async_results = []
        
        # Create a set of interleave indices, to allow each process to only work on
        # one alignment every x alignments, where x is the number of processes running    
        for interleave_ix in range(0, num_processes):
            args = (samfiles, paired_end, interleave_ix,
                    num_processes, mismatch_prob_dict, transition_prob_dict, genome_priors, posterior_cutoff)
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
        combined_separator = merge_separators(unpacked_results)
    
    else:
        # NOT MULTIPROCESSING
        # call the worker procedure within main process, without interleaving
        combined_separator = _worker_procedure(samfiles, paired_end, 0, 1, mismatch_prob_dict,
                                               transition_prob_dict, genome_priors, posterior_cutoff)

    ############################################################################
    # PRINT OUTPUT : All printing occurs here
    ############################################################################
    
    # For each phred, print observed probability of mismatch and number of bases observed in creating that probability
    if mismatch_prob_dict and mismatch_prob_total_values:
        with open(os.path.join(output_dir, 'mismatch_prob_info.txt'), 'w') as outputfile:
            
            for k in mismatch_prob_dict.keys():
                print ('{}\t{}\t{}'.format(k,mismatch_prob_dict[k],mismatch_prob_total_values[k]), file=outputfile)

        with open(os.path.join(output_dir, 'transition_prob_info.txt'), 'w') as outputfile:
            
            for k, v in transition_prob_dict.items():
                qual, base1, base2 = k
                print ('{}\t{}\t{}\t{}'.format(qual, base1, base2, v), file=outputfile)
    
    # Write newly sorted Samfiles, one for each genome
    combined_separator.write_samfiles(separated_samfiles)
        
    # Print a file containing a matrix where columns are genomes and rows contain mismatches per alignment
    mismatches_filepath = os.path.join(output_dir, 'mismatch_counts')
    with open(mismatches_filepath, 'w') as mismatch_file:
        print('\t'.join(genome_names), file=mismatch_file)
        for mismatch_list in combined_separator.mismatch_logs:
            print('\t'.join(str(x) for x in mismatch_list), file=mismatch_file)

    # Print a file containing a matrix where columns are genomes and rows 
    # contain posterior probabilities of genomes per alignment
    posteriors_filepath = os.path.join(output_dir, 'posterior_probabilities')
    with open(posteriors_filepath, 'w') as posteriors_file:
        print('\t'.join(genome_names), file=posteriors_file)
        for posterior_list in combined_separator.mismatch_logs:
            print('\t'.join(str(x) for x in posterior_list),file=posteriors_file)

    # Print a file where each line contains the 'sort fate' for each alignment
    # 0 means unmapped while 1,2,3... mean read is classified as genome1,2,3...
    sort_fate_filepath = os.path.join(output_dir, 'sort_fates')
    with open(sort_fate_filepath, 'w') as sort_fate_file:
        for sort_fate in combined_separator.sort_fates:
            print(sort_fate, file=sort_fate_file)


if __name__ == '__main__':
    
    # Start timer
    t1 = time.time()
    
    # Get command line args
    args = parseargs()
    samfiles = args.samfiles
    paired_end = True if args.paired_end else False
    genome_names = args.names
    num_processes = args.processes
    calculate_mismatches = True if args.calculate_mismatches else False
    mismatch_prob_inputfile = args.mismatch_prob_inputfile if args.mismatch_prob_inputfile else None        
    pileup_height = args.pileup_height
    sample_every = args.sample_every
    genome_priors = args.genome_priors
    posterior_cutoff = args.posterior_cutoff
    output_dir = args.output_dir
    separated_samfiles = args.separated_samfiles
    quiet = True if args.quiet else False
    
    sparta(samfiles, paired_end, genome_names, num_processes, calculate_mismatches,
           pileup_height, sample_every, genome_priors, posterior_cutoff, output_dir,
           separated_samfiles, quiet, mismatch_prob_inputfile)
    
    # Print the total time
    t2 = time.time()
    if not quiet:
        print('TOTAL TIME: {}'.format(t2-t1))