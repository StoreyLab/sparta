# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 12:28:01 2014

@author: peter
"""
from __future__ import print_function
from Bio.Seq import Seq
import collections
# Python 3 compatibility
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

# the following code is horrendous.
# but, pysam seqs are byte objects in python 3.4
# in order to reverse complement, have to convert
# bytes -> string -> Seq -> reverse comp -> string -> bytes
def rev_comp(inseq):
    
    if type(inseq) == str: # python 2.7
        return str(Seq(inseq).reverse_complement())  
    else: # pyton 2.7
        return bytes(str(Seq(inseq.decode(encoding='UTF-8')).reverse_complement()), encoding="UTF-8")
