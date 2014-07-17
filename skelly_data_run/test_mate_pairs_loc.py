# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 13:04:34 2014

@author: peter
"""
from __future__ import print_function
import pysam
import itertools


byfile = pysam.Samfile('by.sam')
rmfile = pysam.Samfile('rm.sam')

for by, rm in itertools.izip(byfile, rmfile):
    
    assert(by.qname == next(byfile).qname)
    assert(rm.qname == next(rmfile).qname)    

    by = next(byfile)
    rm = next(rmfile)
    
print('done')