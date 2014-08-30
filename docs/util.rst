util.py
==================================

.. autofunction:: util.dup_cycle

Credit goes to Adam Rosenfield for this function: http://stackoverflow.com/questions/383565/how-to-iterate-over-a-list-repeating-each-element-in-python

It serves to allow iteration over sort_fates twice per item, so that saved sort fates for paired end reads (where only one parental genome 'fate' is saved for each pair) can be zipped properly to the iterable samfile objects.

.. autofunction:: util.fix_read_mate_order

This function

.. autofunction:: util.compatibility_dict

compatibility_dict is meant to allow collections.defaultdict to be used on both Python 2.7 and 3.4. It simply overrides the items() function of the dict so that in python 2.7 it uses collections.defaultdict.iteritems(self), and in python 3.4 (where the items() function returns an iterator-like object) it simply uses collections.defaultdict.items(self).
