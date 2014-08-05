#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Unit tests for the compare_mappings.py module

Created on Tue Jun 17 16:44:51 2014

@author: Peter Edge
"""

from __future__ import print_function
import unittest
import filecmp
import math
import os
import shutil
import sparta
import sys
import pysam

# take a sequence, quality string, and MD string and return an aligned_read object
def read_gen(seq, qual, MD):
    read = pysam.AlignedRead()
    read.seq = seq
    read.qual = qual
    read.tags = [('MD',MD)]
    return read

# tests for the matched base probability calculator
class test_matched_prob(unittest.TestCase):    
    
    # values were hand-computed in wolframalpha: formula log10(1-(10^((x-33)/-10)))
    # where x was obtained from 'dec' column in the ascii table at http://www.asciitable.com/
    
    def setUp(self):
        self.separator = sparta.multimapped_read_separator(genome_priors=[0.5, 0.5])
    
    def test_F(self):
        prob = self.separator.log10_matched_base_prob[70]
        self.assertAlmostEqual(prob, -0.0000866617872714981203221)

    def test_crunch(self):
        prob = self.separator.log10_matched_base_prob[35]
        self.assertAlmostEqual(prob, -0.4329234333362482973965396)
        
    def test_tilde(self):
        prob = self.separator.log10_matched_base_prob[126]
        self.assertAlmostEqual(prob, -2.1766285001922517006651895e-10)
        

# tests for the mismatched base probability calculator
class test_mismatched_prob(unittest.TestCase):

    # values were hand-computed in wolframalpha: formula log10((10^((x-33)/-10))/3)
    # where x was obtained from the 'dec' column ascii table at http://www.asciitable.com/
    
    def setUp(self):
        self.separator = sparta.multimapped_read_separator(genome_priors=[0.5, 0.5])
    
    def test_F(self):
        prob = self.separator.log10_mismatched_base_prob[70]
        self.assertAlmostEqual(prob, -4.1771212547196624372950279)

    def test_crunch(self):
        prob = self.separator.log10_mismatched_base_prob[35]
        self.assertAlmostEqual(prob, -0.6771212547196624372950279)
        
    def test_tilde(self):
        prob = self.separator.log10_mismatched_base_prob[126]
        self.assertAlmostEqual(prob, -9.7771212547196624372950279) 


# tests for the method that acumulates probabilities for a whole aligned read
class test_aligned_read_prob(unittest.TestCase):
    
    def setUp(self):
        self.separator = sparta.multimapped_read_separator(genome_priors=[0.5, 0.5])

    # test values were hand-computed using wolfram alpha

    # this read is expected to yield a high probability
    def test_likely_read(self):
        read = read_gen('ATGCAAAGGC','FAAAFFF5FF','10')
        prob = self.separator.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 0.98694488490516351120939182)
        
    # this read has matches with probability 0 so result should be 0
    def test_unlikely_read_bc_garbage_quality(self):
        read = read_gen('ATGCAAAGGC','!2222!!111','10')
        prob = self.separator.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 0.0)

    # this read is high quality but has an error; should have low probability
    def test_unlikely_read_bc_high_qual_mismatch(self):
        read = read_gen('ATGCAAAGGC','JJJJJJJJJJ','2A7')
        prob = self.separator.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 0.00002645868511694054719287)

    # this read has two high-quality mismatches and therefore should have very
    # low probability.

    def test_unlikely_read_bc_2_high_qual_mismatch(self):
        read = read_gen('ATGCAAAGGC','JJJJJJJJJJ','0G4G4')
        prob = self.separator.aligned_read_prob(read)
        self.assertAlmostEqual(math.log10(prob), math.log10(7.0061834016176909988601998e-10))

    # slightly longer read, many more mismatches
    def test_longer_unlikely_read(self):
        read = (read_gen('TTTTTTTTTTATGCAAAGGCATGCAAAGGCATGCAAAGGC',
                         'JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ',
                         '0A0A0A0A0A0A0A0A0A0A30'))
        prob = self.separator.aligned_read_prob(read)
        self.assertAlmostEqual(math.log10(prob), math.log10(1.6894778199911023692428023e-46))    
        
    # this read is 60 bp in length
    # half of the read are high-quality mismatches (30 total)
    # this would never happen in real life
    def test_longest_read_many_mismatch(self):
        read = (read_gen('TTTTTTTTTTTTTTTTTTTTATGCAAAGGCTTTTTTTTTTATGCAAAGGCATGCAAAGGC',
                         'JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ',
                         '0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A10A0A0A0A0A0A0A0A0A0A20'))

        prob = self.separator.aligned_read_prob(read)
        # computed on wolframalpha like so:
        # (10^(((74-33)*-0.1)-log10(3)))^30 * (1-10^(((74-33)*-0.1)))^30
        self.assertAlmostEqual(math.log10(prob), math.log10(4.84537506679955566671356787e-138))  
    
    # same test as before but do it 1e5 times to test overhead of prob lookup
    '''
    def test_lookup_overhead(self):
        for i in range(0,10000):
            read = (read_gen('TTTTTTTTTTTTTTTTTTTTATGCAAAGGCTTTTTTTTTTATGCAAAGGCATGCAAAGGC',
                             'JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ',
                             '0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A0A10A0A0A0A0A0A0A0A0A0A20'))
    
            prob = self.separator.aligned_read_prob(read)
            # computed on wolframalpha like so:
            # (10^(((74-33)*-0.1)-log10(3)))^30 * (1-10^(((74-33)*-0.1)))^30
            self.assertAlmostEqual(math.log10(prob), math.log10(4.84537506679955566671356787e-138))
    '''
# tests for the untangle_mappings method, which takes reads that are mapped
# to two different genomes and tries to assign them to the correct one.
class test_untangle(unittest.TestCase):
    
    def setUp(self):
        self.separator = sparta.multimapped_read_separator(genome_priors=[0.5, 0.5])
    
    # variable to hold prob_genome1 in return result.
    NIL = 0.0
    NIL2 = 0
    
    # if neither sequence has errors then it is impossible to map
    # result should be unmapped
    def test_untangle_both_errorfree(self):
        read1 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        read2 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'unclassified')
    
    # if read1 has one error then the result should be classified2
    def test_untangle_one_error(self):
        read1 = read_gen('AAATAAAAAA','FFFFFFFFFF','3T6')
        read2 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'classified2')
        
    # if read1 is all errors then the result should be classified2
    def test_untangle_all_errors(self):
        read1 = read_gen('TTTTTTTTTT','FFFFFFFFFF','0A0A0A0A0A0A0A0A0A0A0')
        read2 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'classified2')
        
    # if read1 has 2 errors but read2 has 4, result should be classified1
    def test_untangle_one_has_more_errors(self):
        read1 = read_gen('AGCAAAAAAA','FFFFFFFFFF','1A0A7')
        read2 = read_gen('AAGGAACTAA','FFFFFFFFFF','2A0A2A0A2')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'classified1')
        
    # if both reads have the same number of errors but genome1's
    # errors have low quality,the result should be classified1
    def test_untangle_based_on_really_low_qual_errors(self):
        read1 = read_gen('AAGCAAAAAA','FF!!FFFFFF','2A0A6')
        read2 = read_gen('AAAAAACTAA','FF!!FFFFFF','6A0A2')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'classified1')

    # if both reads have the same number of errors but genome1's errors have low quality,
    # the result should be genome1
    def test_untangle_based_on_low_qual_errors(self):
        read1 = read_gen('AAGCAAAAAA','FF55FFFFFF','2A0A6')
        read2 = read_gen('AAAAAACTAA','FF55FFFFFF','6A0A2')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'classified1')  
    
    # test posterier cutoff
    # if both reads have the same number of errors and genome1's errors only have
    # a slight dip in quality then read is unsorted based on default posterior
    # probability cutoff of 0.9
    def test_posterior_cutoff(self):
        read1 = read_gen('AAGCAAAAAA','FFBBFFFFFF','2A0A6')
        read2 = read_gen('AAAAAACTAA','FFBBFFFFFF','6A0A2')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'unclassified')
        
    # if both reads have the same number of errors and genome1's errors only have
    # a slight dip in quality then the by severely loosening the posterior
    # cutoff then the result should be genome1
    def test_loosened_posterior_cutoff(self):
        self.separator = sparta.multimapped_read_separator(posterior_cutoff=0.7, genome_priors=[0.5, 0.5])
        read1 = read_gen('AAGCAAAAAA','FFAAFFFFFF','2A0A6')
        read2 = read_gen('AAAAAACTAA','FFAAFFFFFF','6A0A2')
        result, NIL, NIL2 = self.separator.untangle_mappings([read1, read2])
        self.assertTrue(result == 'classified1')

# test that the same result is obtained when sparta is run in single core
# (no child processes) mode and multiprocessing mode
class test_multiprocessing(unittest.TestCase):
    
    def setUp(self):
        if not os.path.exists('unit_test/output'):
            os.makedirs('unit_test/output')
    
    # test single reads on a small (~10K lines) single-read file
    def test_single_reads(self):

        sams = [os.path.join('unit_test','data','single_reads_genome1.sam'), os.path.join('unit_test','data','single_reads_genome2.sam')]
        
        out_dir = os.path.join('unit_test','output')        
        files_to_check = ['sorted1.sam', 'sorted2.sam'] #'supplementary_output.txt'
        
        # no multiprocessing
        sparta.sparta(sams, num_processes=1, output_dir=os.path.join(out_dir, 'no_mp'),
                      separated_samfiles=[os.path.join(out_dir,'no_mp','sorted1.sam'), os.path.join(out_dir,'no_mp','sorted2.sam')], quiet=True)
        
        # processes = cpu_count
        sparta.sparta(sams, output_dir=os.path.join(out_dir,'mp'),
                      separated_samfiles=[os.path.join(out_dir,'mp','sorted1.sam'), os.path.join(out_dir,'mp','sorted2.sam')], quiet=True)

        # processes = 10
        sparta.sparta(sams, num_processes=10, output_dir=os.path.join(out_dir,'mp10'),
                      separated_samfiles=[os.path.join(out_dir,'mp10','sorted1.sam'), os.path.join(out_dir,'mp10','sorted2.sam')], quiet=True)
                      
        match, mismatch, errors = filecmp.cmpfiles(os.path.join(out_dir,'no_mp'), os.path.join(out_dir,'mp'),
                                                   common=files_to_check)
                                                   
        match2, mismatch2, errors2 = filecmp.cmpfiles(os.path.join(out_dir,'no_mp'), os.path.join(out_dir,'mp10'),
                                                      common=files_to_check)

        # check that output was created, and no files mismatched, and no errors occured
        assert match != []
        assert mismatch == []
        assert errors == []
        assert match2 != []
        assert mismatch2 == []
        assert errors2 == []

        # clean up
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)


    # test on a small (~10K lines) paired-end read file
    def test_paired_end_reads(self):
        
        sams = [os.path.join('unit_test','data','paired_end_reads_genome1.sam'), os.path.join('unit_test','data','paired_end_reads_genome2.sam')]
        
        out_dir = os.path.join('unit_test','output')        
        files_to_check = ['sorted1.sam', 'sorted2.sam'] # , 'supplementary_output.txt'
        
        # no multiprocessing
        sparta.sparta(sams, paired_end=True, num_processes=1, output_dir=os.path.join(out_dir, 'no_mp'),
                      separated_samfiles=[os.path.join(out_dir,'no_mp','sorted1.sam'),os.path.join(out_dir,'no_mp','sorted2.sam')], quiet=True)
        
        # processes = cpu_count
        sparta.sparta(sams, paired_end=True, output_dir=os.path.join(out_dir,'mp'),
                      separated_samfiles=[os.path.join(out_dir,'mp','sorted1.sam'), os.path.join(out_dir,'mp','sorted2.sam')], quiet=True)

        # processes = 10
        sparta.sparta(sams, paired_end=True, num_processes=10, output_dir=os.path.join(out_dir,'mp10'),
                      separated_samfiles=[os.path.join(out_dir,'mp10','sorted1.sam'), os.path.join(out_dir,'mp10','sorted2.sam')], quiet=True)
                      
        match, mismatch, errors = filecmp.cmpfiles(os.path.join(out_dir,'no_mp'), os.path.join(out_dir,'mp'),
                                                   common=files_to_check)
                                                   
        match2, mismatch2, errors2 = filecmp.cmpfiles(os.path.join(out_dir,'no_mp'), os.path.join(out_dir,'mp10'),
                                                      common=files_to_check)

        # check that output was created, and no files mismatched, and no errors occured
        assert match != []
        assert mismatch == []
        assert errors == []
        assert match2 != []
        assert mismatch2 == []
        assert errors2 == []
        
        # clean up
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)

if __name__ == '__main__':

    # this try/catch block allows test suite to be run in spyder IDE interpreter
    # without hanging forever afterwards
    #via: http://stackoverflow.com/questions/9202772/tests-succeed-still-get-traceback
    try:
        unittest.main()
    except SystemExit as inst:
        pass
    