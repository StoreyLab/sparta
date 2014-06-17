# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 16:44:51 2014

@author: Peter Edge
"""

import unittest

class test_compare_mappings(unittest.TestCase):
    pass

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(test_compare_mappings)
    unittest.TextTestRunner(verbosity=2).run(suite)