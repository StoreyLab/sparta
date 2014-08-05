# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 17:46:09 2014

@author: peter
"""
from compatibility import rev_comp

def fix_read_mate_order(aligned_reads, aligned_read_mates, cautious=False):
    
    # qname field should match for all 4 alignedread objects
    # because really it is 2 copies of the same read+mate pair
    assert aligned_reads[0].qname == aligned_read_mates[0].qname   
    
    # check that across all genomes' alignments, qname for this read is the same
    for aligned_read, aligned_read_mate in zip(aligned_reads, aligned_read_mates):
        assert aligned_read.qname == aligned_reads[0].qname
        assert aligned_read_mate.qname == aligned_read_mates[0].qname
    
    # if the reads are the same length as their mates, we must judge based on their sequence
    if len(aligned_reads[0].seq) == len(aligned_read_mates[0].seq):
        cautious=True    
    
    fixed_aligned_reads = []
    fixed_aligned_read_mates = []    
    
    # ASSUMING THAT THE SAM FILES WERE GENERATED USING THE SAME FASTQ FILE,
    # CAUTIOUS MODE IS NOT NECESSARY.
    # cautious mode checks that the aligned1 sequence and aligned2 sequence match.
    # this may involve reverse complementation. This is obviously SLOWER.
    # Should only be necessary in the case of paired end reads where mates have identical length
    
    if cautious:
        # The RNA reads in aligned1 might be flip-flopped
        # (As in, aligned1 actually refers to aligned2_mate)
        # in this case, switch aligned2 and aligned2_mate
    
        # make all the reads match the ordering of the first read/mate pair
        # use sequence to check if reads/mates are swapped

        for read, mate in zip(aligned_reads, aligned_read_mates):
            
            if read.is_reverse == aligned_reads[0].is_reverse:
                
                if read.seq == aligned_reads[0].seq:

                    fixed_aligned_reads.append(read)
                    fixed_aligned_read_mates.append(mate)
                else:
                    fixed_aligned_reads.append(mate)
                    fixed_aligned_read_mates.append(read)
            
            else:
                # one read is reversed, need to revcomp one to see if equal                    
                
                read_revcomp = rev_comp(read.seq)
                
                if read_revcomp == aligned_reads[0].seq:
                    
                    fixed_aligned_reads.append(read)
                    fixed_aligned_read_mates.append(mate)
                else:
                    fixed_aligned_reads.append(mate)
                    fixed_aligned_read_mates.append(read)
                    
        for read, mate in zip(fixed_aligned_reads, fixed_aligned_read_mates):
            
            if read.is_reverse == fixed_aligned_reads[0].is_reverse:
                assert read.seq == fixed_aligned_reads[0].seq
            else:
                assert read.seq == rev_comp(fixed_aligned_reads[0].seq)
                
            if mate.is_reverse == fixed_aligned_read_mates[0].is_reverse:
                assert mate.seq == fixed_aligned_read_mates[0].seq
            else:
                assert mate.seq == rev_comp(fixed_aligned_read_mates[0].seq)        
                
    else:
        # shouldn't be necessary to check seqs are identical
        # make all the reads match the ordering of the first read/mate pair
        # use length as the metric to quickly check if reads/mates are swapped
        reads_len = len(aligned_reads[0].seq)
        mates_len = len(aligned_read_mates[0].seq)

        for read, mate in zip(aligned_reads, aligned_read_mates):
            if len(read.seq) == reads_len:
                fixed_aligned_reads.append(read)
                fixed_aligned_read_mates.append(mate)
            else:
                fixed_aligned_read_mates.append(read)
                fixed_aligned_reads.append(mate)
            
        for read, mate in zip(fixed_aligned_reads, fixed_aligned_read_mates):
            assert len(read.seq) == reads_len
            assert len(mate.seq) == mates_len
        
    return fixed_aligned_reads, fixed_aligned_read_mates