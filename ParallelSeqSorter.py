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
This file contains the function compare_mappings(), which takes in two SAM
alignment files which map the same RNAseq reads to two different parental
genomes. The RNAseq reads should map to either one parental allele type or the
other. This program sorts the reads based on the errors to each genome.
'''
'''
COMMAND LINE USAGE:
python2 SeqSorter.py <samfile1> <samfile2> <optional arguments>

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

EXAMPLE:

python2 SeqSorter.py sample_data/S288C_bowtie_test/BY_bowtie_out.sam \
    sample_data/RM_bowtie_test/RM_bowtie_out.sam -n1 BY -n2 RM \
    -v verbose_outfile

'''

import argparse
import collections
import countErrorOccurences
import math
import multiprocessing as mp
import pysam
import re
import sys
import time

# global variables
estimate_error_prob = False

# parse SeqSorter program input.
def parseargs():
    
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument('samfile1', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)
    parser.add_argument('samfile2', nargs='?', type = str, help='path to samfile 2', default=sys.stdin)
    # optional nicknames for the genomes used in the two samfiles (recommended)    
    parser.add_argument('-n1', '--name1', nargs='?', type = str, help='name for genome 1 (reference for samfile 1)', default=sys.stdin)
    parser.add_argument('-n2', '--name2', nargs='?', type = str, help='name for genome 2 (reference for samfile 2)', default=sys.stdin)
    # this is an optional verbose outfile with useful values for statistical analysis
    parser.add_argument('-v', '--verbose_file', nargs='?', type=argparse.FileType('w'), help='specify a filename for verbose output (useful for statistical analysis)', default=sys.stdin)    

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# sorter class that sorts an RNAseq read to one parental allele or the other
class multimapped_read_sorter():
    '''
    The following methods compute the probabilities that matched and mismatched bases contribute.
    Let M be the probability of a miscalled base, computed as:
        10^((number specified by phred char)/-10)
        NOTE: the ascii phred char is offset by 33 for sanger format reads
    A matched base contributes (1 - M): the probability that the base call was correct
    and the genome in question generated the observed base.
    A mismatched base contributes M / 3: the probability that the base call was wrong
    and the called base was the observed 1 of 3 possible other bases.
    '''
  
    # lists that hold precomputed matched and missed base probabilities
    # index into the list using the dec value of the ascii char, which is simply
    # the value in the qual string once it is represented as a byte array
     
    log10_matched_base_prob = [0]*127
    log10_mismatched_base_prob = [0]*127
    
    # the init method fills the log10 match/mismatch probability lists
    def __init__(self, mismatch_prob_dict, verbose):
        
        # if verbose file was specified then run in verbose mode

        self.verbose = verbose
        self.category_counter = collections.defaultdict(int)
        self.logs = []
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
    
            print ('phred\tmatch_prob\tmismatch_prob')

            for i in range(33,127):
                print ('{}\t{}\t{}'.format(i,self.log10_matched_base_prob[i],self.log10_mismatched_base_prob[i]))
        else:
            print('no empirical match/mismatch vals')
            
    def log(self, err1, err2, prob1, category):
        if self.verbose:
            msg = '{}\t{}\t{}\t{}'.format(err1, err2, prob1, category)        
            self.logs.append(msg)

    # compute probability of a whole RNAseq read being generated by a genome, 
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
    def untangle_two_samfiles(self, sam1, sam2, interleave_ix, num_processes, genome1_name, genome2_name, genome1_prior=0.5, posterior_cutoff=0.9, verbose=False):
        
        ix = 0
        for aligned1, aligned2 in zip(sam1, sam2):
            
            if (ix % num_processes != interleave_ix):
                ix += 1
                continue
            ix += 1

            assert aligned1.qname == aligned2.qname
            
            if aligned1.is_unmapped and aligned2.is_unmapped:
                #this read is probably junk
                self.category_counter['no_match'] += 1
                if verbose:
                    self.log('NA','NA','NA', 'unmapped')
                    
            elif (not aligned1.is_unmapped) and aligned2.is_unmapped:
                #maps to alignment1 but not alignment2; either junk or
                #an unshared gene bw BY and RM
                self.category_counter['match1'] += 1
                if verbose:
                    num_err1 = len(re.findall(self.MD_REGEX, aligned1.opt("MD")))
                    self.log(num_err1, 'NA', '1', 'mapped {}'.format(genome1_name))
            
            elif aligned1.is_unmapped and (not aligned2.is_unmapped):
                #maps to alignment2 but not alignment1; either junk or
                #an unshared gene bw BY and RM
                self.category_counter['match2'] += 1
                if verbose:
                    num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                    self.log('NA', num_err2, '0', 'mapped {}'.format(genome2_name))
            
            elif aligned1.opt("MD") == aligned2.opt("MD"):
                # same errors
                self.category_counter['same_errors'] += 1
                if verbose:
                    num_err1 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                    num_err2 = len(re.findall(self.MD_REGEX, aligned2.opt("MD")))
                    self.log(num_err1, num_err2, '0.5', 'unclassified: same errors')
            
            else:
                # bowtie matched the read to both genomes
                # use sorter to map read to a genome
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
            
                '''        
                ##############################################################
                # Analysis of unequal deletions (maybe remove?)
                ##############################################################
                
                # keep track of whether the classified genome had more deleted seqs/bases
                more_del_seqs = False
                more_del_bases = False
                
                # count number of deleted sequences (number of carets)
                num_del_seqs_genome1 = aligned1.opt("MD").count('^')
                num_del_seqs_genome2 = aligned2.opt("MD").count('^')       
                
                if num_del_seqs_genome1 != num_del_seqs_genome2:
                    self.category_counter['diff_number_del_seqs'] += 1
                
                # check if a classified (winner) mapping has more deleted sequences than the loser
                if (((most_likely_genome == "genome1") and (num_del_seqs_genome1 > num_del_seqs_genome2)) or
                    ((most_likely_genome == "genome2") and (num_del_seqs_genome2 > num_del_seqs_genome1))):
                    del_seqs_mattered += 1
                    more_del_seqs = True
                
                # find all deletions
                del1 = re.findall(DEL_REGEX, aligned1.opt("MD"))
                del2 = re.findall(DEL_REGEX, aligned2.opt("MD"))
                
                total_del1 = 0
                total_del2 = 0
                
                # count the total number of deleted bases in aligned read 1 and 2
                for deletion in del1:
                    total_del1 += (len(deletion) - 1)
                
                for deletion in del2:
                    total_del2 += (len(deletion) - 1)
                    
                if total_del1 != total_del2:
                    diff_number_del_bases += 1
                
                # check if a classified (winner) mapping has more deleted sequences than the loser
                if (((most_likely_genome == "genome1") and (total_del1 > total_del2)) or
                    ((most_likely_genome == "genome2") and (total_del2 > total_del1))):
                    del_bases_mattered += 1
                    more_del_bases = True
                
                # check if having more deleted bases is unrelated to having more deleted seqs
                if (more_del_seqs != more_del_bases):
                    bases_seqs_unrel += 1
            
                ##############################################################
                # End analysis of unequal deletions 
                ##############################################################
                '''
def worker_procedure(samfile1, samfile2, interleave_ix, num_processes, mismatch_prob_dict, genome1_name, genome2_name, genome1_prior, posterior_cutoff, verbose):

    sorter = multimapped_read_sorter(mismatch_prob_dict, verbose)
    sam1 = pysam.Samfile(samfile1)
    sam2 = pysam.Samfile(samfile2)
    sorter.untangle_two_samfiles(sam1, sam2, interleave_ix, num_processes, genome1_name, genome2_name, genome1_prior, posterior_cutoff, verbose)
    return sorter

# merge a list of sorters into one.
def merge_sorters(sorter_list):
    new_sorter = sorter_list.pop()
    print(new_sorter.category_counter)

    # merge info
    while(sorter_list):
        old_sorter = sorter_list.pop()
        new_sorter.logs.append(old_sorter.logs)
        
        print(old_sorter.category_counter)

        for category, count in old_sorter.category_counter.iteritems():
            new_sorter.category_counter[category] += count

    return new_sorter

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
    verbose_file = args.verbose_file    

    genome1_prior = 0.5
    posterior_cutoff = 0.9
    
    if estimate_error_prob:
        # currently broken
        mismatch_prob_dict = countErrorOccurences.count_error_occurrences(samfile1, samfile2, genome1_name, genome2_name, 'count_error_occurences_out')
    else:
        mismatch_prob_dict = None    
    
    if verbose_file:
        verbose = True
    else:
        verbose = False
        
    # get the number of available CPUs
    # default to 1
    num_processes = 1
    try:
        num_processes = mp.cpu_count()
    except:
        pass
    
    worker_pool = mp.Pool(processes=num_processes)
    
    async_results = []    
    
    for interleave_ix in range(0, num_processes):
        args = (samfile1, samfile2, interleave_ix, num_processes, mismatch_prob_dict, genome1_name, genome2_name, genome1_prior, posterior_cutoff, verbose)
        async_results.append(worker_pool.apply_async(worker_procedure, args))
        
    worker_pool.close()
    worker_pool.join()
    
    unpacked_results = []
    for result in async_results:
        result.wait()
        unpacked_results.append(result.get())
        
    total_result = merge_sorters(unpacked_results)
    category_counter = total_result.category_counter
    # counters for different relationships between aligned reads
    
    # print counts of each scenario    
    print('\n{}\tunmapped by bowtie to either {} or {}'.format(category_counter['no_match'], genome1_name, genome2_name))
    print('{}\tmapped by bowtie to either {} or {}'.format(category_counter['match1'], genome1_name, genome2_name))
    print('{}\tmapped by bowtie to both {} or {}, but errors are same'.format(category_counter['same_errors'], genome1_name, genome2_name))
    total_mapped_to_both = category_counter['classified1'] + category_counter['classified2'] + category_counter['unclassified']
    print('{}\tmapped by bowtie to both {} and {}:'.format(total_mapped_to_both, genome1_name, genome2_name))
    print('   {}\tassigned to {} based on errors'.format(category_counter['classified1'], genome1_name))
    print('   {}\tassigned to {} based on errors'.format(category_counter['classified2'], genome2_name))
    print('   {}\tunmapped based on errors'.format(category_counter['unclassified']))
    
    t2 = time.time()
    print('TOTAL TIME: {}'.format(t2-t1))
    
    '''    
    print
    print('{}\tcases of difference in number of deleted seqs'.format(diff_number_del_seqs))
    print('{}\tcases of difference in number of deleted bases'.format(diff_number_del_bases))
    print('   {}\twinners had more deleted seqs'.format(del_seqs_mattered))
    print('   {}\twinners had more deleted bases'.format(del_bases_mattered))
    print('   {}\tcases of the above 2 being different'.format(bases_seqs_unrel))
    '''
    # close output files
    verbose_file.close()
if __name__ == '__main__':
    main()