# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 17:46:09 2014

@author: peter
"""

def fix_read_mate_order(aligned_reads, aligned_read_mates):
    
    # qname field should match for all alignedread objects
    # because really it is 2 copies of the same read+mate pair
    assert aligned_reads[0].qname == aligned_read_mates[0].qname   
    
    # check that across all genomes' alignments, qname for this read is the same
    for aligned_read, aligned_read_mate in zip(aligned_reads, aligned_read_mates):
        assert aligned_read.qname == aligned_reads[0].qname
        assert aligned_read_mate.qname == aligned_read_mates[0].qname
        
        if aligned_read.is_read2:
            temp = aligned_read
            aligned_read = aligned_read_mate
            aligned_read_mate = temp
        
        assert aligned_read.is_read1
        assert aligned_read.is_read2
            
    return aligned_reads, aligned_read_mates