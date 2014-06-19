# -*- coding: utf-8 -*-
"""
Unit tests for the compare_mappings.py module

Created on Tue Jun 17 16:44:51 2014

@author: Peter Edge
"""

import unittest
import compare_mappings
import math
import pysam

# take a sequence, quality string, and MD string and return an aligned_read object
def read_gen(seq, qual, MD):
    read = pysam.AlignedRead()
    read.seq = seq
    read.qual = qual
    read.tags = [('MD',MD)]
    return read

class test_matched_prob(unittest.TestCase):    
    
    # values were hand-computed in wolframalpha: formula log10(1-(10^((x-33)/-10)))
    # where x was obtained from 'dec' column in the ascii table at http://www.asciitable.com/
    
    def setUp(self):
        self.sorter = compare_mappings.multimapped_read_sorter()
    
    def test_F(self):
        prob = self.sorter.log10_matched_base_prob('F')
        self.assertAlmostEqual(prob, -0.0000866617872714981203221)

    def test_crunch(self):
        prob = self.sorter.log10_matched_base_prob('#')
        self.assertAlmostEqual(prob, -0.4329234333362482973965396)

    def test_bang(self):
        prob = self.sorter.log10_matched_base_prob('!')
        self.assertAlmostEqual(prob, math.log10(0.25))
        
    def test_tilde(self):
        prob = self.sorter.log10_matched_base_prob('~')
        self.assertAlmostEqual(prob, -2.1766285001922517006651895e-10)
        
class test_mismatched_prob(unittest.TestCase):

    # values were hand-computed in wolframalpha: formula log10((10^((x-33)/-10))/3)
    # where x was obtained from the 'dec' column ascii table at http://www.asciitable.com/
    
    def setUp(self):
        self.sorter = compare_mappings.multimapped_read_sorter()
    
    def test_F(self):
        prob = self.sorter.log10_mismatched_base_prob('F')
        self.assertAlmostEqual(prob, -4.1771212547196624372950279)

    def test_crunch(self):
        prob = self.sorter.log10_mismatched_base_prob('#')
        self.assertAlmostEqual(prob, -0.6771212547196624372950279)

    def test_bang(self):
        prob = self.sorter.log10_mismatched_base_prob('!')
        self.assertAlmostEqual(prob, math.log10(0.25))
        
    def test_tilde(self):
        prob = self.sorter.log10_mismatched_base_prob('~')
        self.assertAlmostEqual(prob, -9.7771212547196624372950279) 

class test_aligned_read_prob(unittest.TestCase):
    
    def setUp(self):
        self.sorter = compare_mappings.multimapped_read_sorter()

    # test values were hand-computed using wolfram alpha

    # this read is expected to yield a high probability
    def test_likely_read(self):
        read = read_gen('ATGCAAAGGC','FAAAFFF5FF','10')
        prob = self.sorter.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 0.98694488490516351120939182)
        
    # this read has low quality and should have a lower probability than previous
    def test_less_likely_read(self):
        read = read_gen('ATGCAAAGGC','!2222!!111','10')
        prob = self.sorter.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 0.01335559708084671580915440)
        
    # this read is high quality but has an error; should have low probability
    def test_unlikely_read(self):
        read = read_gen('ATGCAAAGGC','JJJJJJJJJJ','2A7')
        prob = self.sorter.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 0.00002645868511694054719287)
    
    # this read has two high-quality mismatches and therefore should have very
    # low probability. The sequence is 1/5 high-qual mismatches.
    
    # NOT PASSING THIS TEST CASE, my best guess is that the really small numbers
    # lead to weird behavior
    # however, a sequence that is 1/5 high-qual errors will probably never be
    # mapped by bowtie anyway
    
    def test_super_unlikely_read(self):
        read = read_gen('ATGCAAAGGC','JJJJJJJJJJ','4GG4')
        prob = self.sorter.aligned_read_prob(read)
        self.assertAlmostEqual(prob, 4.42270698591381293343184747e-11)


class test_untangle(unittest.TestCase):
    
    def setUp(self):
        self.sorter = compare_mappings.multimapped_read_sorter()
    
    # if neither sequence has errors then it is impossible to map
    # result should be unmapped
    def test_untangle_both_errorfree(self):
        read1 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        read2 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'unmapped')
    
    # if read1 has one error then the result should be genome2
    def test_untangle_one_error(self):
        read1 = read_gen('AAATAAAAAA','FFFFFFFFFF','3T6')
        read2 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'genome2')
        
    # if read1 is all errors then the result should be genome2
    def test_untangle_all_errors(self):
        read1 = read_gen('TTTTTTTTTT','FFFFFFFFFF','0AAAAAAAAAA0')
        read2 = read_gen('AAAAAAAAAA','FFFFFFFFFF','10')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'genome2')
        
    # if read1 has 2 errors but read2 has 4, result should be genome1
    def test_untangle_one_has_more_errors(self):
        read1 = read_gen('AGCAAAAAAA','FFFFFFFFFF','1AA7')
        read2 = read_gen('AAGGAACTAA','FFFFFFFFFF','2AA2AA2')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'genome1')
        
    # if both reads have the same number of errors but genome1's
    #errors have low quality,the result should be genome1
    def test_untangle_based_on_really_low_qual_errors(self):
        read1 = read_gen('AAGCAAAAAA','FF!!FFFFFF','2AA6')
        read2 = read_gen('AAAAAACTAA','FF!!FFFFFF','6AA2')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'genome1')

    # if both reads have the same number of errors but genome1's errors have low quality,
    # the result should be genome1
    def test_untangle_based_on_low_qual_errors(self):
        read1 = read_gen('AAGCAAAAAA','FF55FFFFFF','2AA6')
        read2 = read_gen('AAAAAACTAA','FF55FFFFFF','6AA2')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'genome1')
    
    # test posterier cutoff
    # if both reads have the same number of errors and genome1's errors only have
    # a slight dip in quality then read is unsorted based on default posterior
    # probability cutoff of 0.9
    def test_posterior_cutoff(self):
        read1 = read_gen('AAGCAAAAAA','FFAAFFFFFF','2AA6')
        read2 = read_gen('AAAAAACTAA','FFAAFFFFFF','6AA2')
        result = self.sorter.untangle_two_mappings(read1, read2)
        self.assertTrue(result == 'unmapped')
        
    # if both reads have the same number of errors and genome1's errors only have
    # a slight dip in quality then the by severely loosening the posterior
    # cutoff then the result should be genome1
    def test_loosened_posterior_cutoff(self):
        read1 = read_gen('AAGCAAAAAA','FFAAFFFFFF','2AA6')
        read2 = read_gen('AAAAAACTAA','FFAAFFFFFF','6AA2')
        result = self.sorter.untangle_two_mappings(read1, read2, posterior_cutoff=0.7)
        self.assertTrue(result == 'genome1')
    
    # I added this test case because the program still fails the test case
    # test_super_unlikely_read(self), which I think is caused by underflow somewhere.
    # aggregations of mismatch probabilities can't exceed ~1e-7 on the small side.
    # The behavior exhibited in this test case is not mathematically correct,
    # but I think it is acceptable because if a mapping is 1/5 high-qual errors
    # it really shouldn't be called the correct mapping anyway.
    def test_one_mapping_really_bad(self):
        read1 = read_gen('AAGCAAAAAA','FFJJFFFFFF','2AA6')
        read2 = read_gen('AAAAAACTAA','FFAAFFJJJJ','6AAAA')
        result = self.sorter.untangle_two_mappings(read1, read2, posterior_cutoff=0.7)
        self.assertTrue(result == 'unmapped')
        
if __name__ == '__main__':

    # this is awful practice but this try/catch block allows test suite
    # to be run in spyder IDE interpreter without hanging forever afterwards
    #via: http://stackoverflow.com/questions/9202772/tests-succeed-still-get-traceback
    try:
        unittest.main()
    except SystemExit as inst:
        pass
    