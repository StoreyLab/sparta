# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 17:46:09 2014

@author: peter
"""
from compatibility import rev_comp

def fix_read_mate_order(aligned1, aligned2, aligned1_mate, aligned2_mate, cautious=False):
    # qname field should match for all 4 alignedread objects
    # because really it is 2 copies of the same read+mate pair
    assert aligned1.qname == aligned2.qname
    assert aligned1.qname == aligned1_mate.qname
    assert aligned1_mate.qname == aligned2_mate.qname                
        
    if (len(aligned1.seq) == len(aligned1_mate.seq) and
        len(aligned2.seq) == len(aligned2_mate.seq)):
        cautious=True    
    
    # ASSUMING THAT THE SAM FILES WERE GENERATED USING THE SAME FASTQ FILE,
    # CAUTIOUS MODE IS NOT NECESSARY. 
    # and if they weren't, then the whole thing is messed up anyway because that
    # is pretty much the one big assumption made by SPARTA.
    # cautious mode checks that the aligned1 sequence and aligned2 sequence match.
    # this may involve reverse complementation
    # overall this is very slow and should never be necessary
    # except in the case of paired end reads where mates have identical length
    
    if cautious:
        # The RNA reads in aligned1 might be flip-flopped
        # (As in, aligned1 actually refers to aligned2_mate)
        # in this case, switch aligned2 and aligned2_mate
                        
        if aligned1.is_reverse == aligned2.is_reverse:
            
            if aligned1.seq != aligned2.seq:
                # we have the mate instead
                # switch aligned2 with its mate
                temp = aligned2
                aligned2 = aligned2_mate
                aligned2_mate = temp
            
            if aligned1.is_reverse == aligned2.is_reverse:
                assert aligned1.seq == aligned2.seq
            else:
                assert rev_comp(aligned1.seq) == aligned2.seq
            
            if aligned1_mate.is_reverse == aligned2_mate.is_reverse:
                assert aligned1_mate.seq == aligned2_mate.seq
            else:
                assert rev_comp(aligned1_mate.seq) == aligned2_mate.seq
        
        else:
            # one read is reversed, need to revcomp one to see if equal                    
            
            aligned1_revcomp = rev_comp(aligned1.seq)
            
            if aligned1_revcomp != aligned2.seq:
                # we have the mate instead
                # switch aligned2 with its mate
                temp = aligned2
                aligned2 = aligned2_mate
                aligned2_mate = temp
            
            if aligned1.is_reverse == aligned2.is_reverse:
                assert aligned1.seq == aligned2.seq
            else:
                assert aligned1_revcomp == aligned2.seq
            
        
            if aligned1_mate.is_reverse == aligned2_mate.is_reverse:
                assert aligned1_mate.seq == aligned2_mate.seq
            else:
                assert rev_comp(aligned1_mate.seq) == aligned2_mate.seq
    
    else:
        
        # just check that the lengths are correct
        if len(aligned1.seq) != len(aligned2.seq):
            # we have the mate instead
            # switch aligned2 with its mate
            temp = aligned2
            aligned2 = aligned2_mate
            aligned2_mate = temp
            
        assert len(aligned1.seq) == len(aligned2.seq)
        assert len(aligned1_mate.seq) == len(aligned2_mate.seq)
        
    return aligned1, aligned2, aligned1_mate, aligned2_mate