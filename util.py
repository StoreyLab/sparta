# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 17:46:09 2014

@author: peter
"""

# this clever function via Adam Rosenfield: http://stackoverflow.com/questions/383565/how-to-iterate-over-a-list-repeating-each-element-in-python
def dup_cycle(iterable):
    while True:
        for item in iterable:
            yield item
            yield item

def fix_read_mate_order(aligned_reads, aligned_read_mates):
    
    # qname field should match for all alignedread objects
    # because really it is 2 copies of the same read+mate pair    
    assert aligned_reads[0].qname == aligned_read_mates[0].qname   

    # store the corrected lists to output
    fixed_aligned_reads = []
    fixed_aligned_read_mates = []
    
    # check that across all genomes' alignments, qname for this read is the same
    for aligned_read, aligned_read_mate in zip(aligned_reads, aligned_read_mates):
        assert aligned_read.qname == aligned_reads[0].qname
        assert aligned_read_mate.qname == aligned_read_mates[0].qname
        
        if aligned_read.is_read1:
            fixed_aligned_reads.append(aligned_read)
            fixed_aligned_read_mates.append(aligned_read_mate)
        else:
            fixed_aligned_reads.append(aligned_read_mate)
            fixed_aligned_read_mates.append(aligned_read)
    
    for aligned_read, aligned_read_mate in zip(fixed_aligned_reads, fixed_aligned_read_mates):
        assert aligned_read.is_read1
        assert aligned_read_mate.is_read2
            
    return fixed_aligned_reads, fixed_aligned_read_mates