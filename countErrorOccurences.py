# -*- coding: utf-8 -*-

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
import collections
import pysam
import sys

def parseargs():
    
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument('samfile1', nargs='?', type = str, help='path to samfile 1', default=sys.stdin)
    parser.add_argument('samfile2', nargs='?', type = str, help='path to samfile 2', default=sys.stdin)
    # optional nicknames for the genomes used in the two samfiles (recommended)    
    parser.add_argument('-n1', '--name1', nargs='?', type = str, help='name for genome 1 (reference for samfile 1)', default=sys.stdin)
    parser.add_argument('-n2', '--name2', nargs='?', type = str, help='name for genome 2 (reference for samfile 2)', default=sys.stdin)
    # this is an optional verbose outfile with useful values for statistical analysis
    parser.add_argument('-v', '--output_file', nargs='?', type=argparse.FileType('w'), help='specify a filename for verbose output (useful for statistical analysis)', default=sys.stdin)    

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

'''
class consensus_dict(collections.defaultdict):
    
    # if key is missing then perform a pileup, determine the consensus base,
    # and return the consensus base or NULL if unclear or unmatched to BY/RM
    def __missing__(key):
        chrom1, chrom2, pos1, pos2, samfile1 = key
        for pileupcolumn in samfile.pileup(chr1, 100, 120):
            
            for pileupread in pileupcolumn.pileups:
        
        return
'''
        
def count_error_occurrences(samfile1, samfile2, genome1_name, genome2_name, output_file):
    
    sam1 = pysam.Samfile(samfile1)
    sam2 = pysam.Samfile(samfile2)
    
    results = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    
    for aligned1, aligned2 in zip(sam1, sam2):

        assert aligned1.qname == aligned2.qname
        
        if aligned1.is_unmapped and aligned2.is_unmapped:
            # both unmapped 
            pass
        
        elif not aligned1.is_unmapped and aligned2.is_unmapped:
            # only alignment1 mapped            
            pass
        
        elif aligned1.is_unmapped and not aligned2.is_unmapped:
            # only alignment2 mapped
            pass
        else:
            
            assert len(aligned1.seq) == len(aligned2.seq)
            # both alignments mapped
            pos_dict1 = dict(aligned1.aligned_pairs)
            pos_dict2 = dict(aligned2.aligned_pairs)
                        
            for i in range(0, len(aligned1.seq)):
                            
                if aligned1.seq[i] == aligned2.seq[i]:
                    
                    chrom1 = sam1.getrname(aligned1.tid)
                    pos1 = pos_dict1[i]      
                    chrom2 = sam2.getrname(aligned2.tid)
                    pos2 = pos_dict2[i]
                    
                    results[(chrom1, chrom2, pos1, pos2)][aligned1.seq[i]][aligned1.qual[i]] += 1
            

    return results

# main logic
# call compare_mappings() on samfile1 and samfile2 from standard input 
def main():
    
    # get command line args
    args = parseargs()
    samfile1 = args.samfile1
    samfile2 = args.samfile2
    genome1_name = args.name1
    genome2_name = args.name2
    output_file = args.output_file
    
    # compare mappings between samfiles
    results = count_error_occurrences(samfile1, samfile2, genome1_name, genome2_name, output_file)
    
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
        
        for base, count in base_count.iteritems():
            if count > cutoff * total_bases:
                consensus = base
               
               
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
                    
    for qual, mismatch_count in quality_score_mismatch_counter.iteritems():
        mismatch_prob = mismatch_count * 1.0 / ((mismatch_count + quality_score_match_counter[qual])*1.0)
        print('{}\t{}'.format(qual,mismatch_prob))
        
    # close output files
    output_file.close()
    
if __name__ == '__main__':
    main()