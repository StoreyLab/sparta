# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 17:46:09 2014

@author: peter
"""

from __future__ import print_function
from Bio.Seq import Seq
import collections

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

################################################################################
## Python 3 compatibility specific code
################################################################################

# These functions simply exist to allow SPARTA to be cross-compatible with Python 2.7 and 3.4

# in Python 3 itertools.izip does not exist because the basic zip returns an iterator
# try to import itertools.izip and if it isn't there (Python 3), use zip instead
# credit to: http://codereview.stackexchange.com/questions/26271/cleaner-way-to-import-izip-for-different-versions-of-python
try:
    from itertools import izip
except ImportError:  #python3.x
    izip = zip

# subclass of defaultdict to maintain compatibility
# in python 2 iteritems returns an iterator
# in python 3 iteritems doesn't exist, items returns a view which is like an iterator
class compatibility_dict(collections.defaultdict):
    
    def items(self):
        try:
            return collections.defaultdict.iteritems(self)
        except AttributeError:   #python3.x
            return collections.defaultdict.items(self)
