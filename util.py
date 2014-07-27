# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 17:46:09 2014

@author: peter
"""
from compatibility import rev_comp

def fix_read_mate_order(aligned1, aligned2, aligned1_mate, aligned2_mate):
    # qname field should match for all 4 alignedread objects
    # because really it is 2 copies of the same read+mate pair
    assert aligned1.qname == aligned2.qname
    assert aligned1.qname == aligned1_mate.qname
    assert aligned1_mate.qname == aligned2_mate.qname                
    
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
            
    return aligned1, aligned2, aligned1_mate, aligned2_mate