util.py
==================================

.. autofunction:: util.dup_cycle

Credit goes to Adam Rosenfield for this function: http://stackoverflow.com/questions/383565/how-to-iterate-over-a-list-repeating-each-element-in-python

It serves to allow iteration over sort_fates twice per item, so that saved sort fates for paired end reads (where only one parental genome 'fate' is saved for each pair) can be zipped properly to the iterable samfile objects.

.. autofunction:: util.fix_read_mate_order

This function is designed to circumvent the fact that in the case of samfiles of paired end reads, there is no guarantee that reads and their mates will appear in the same order that they appeared in the fastq file. For example, imagine a set of paired-end reads is aligned to two different genomes with bowtie2, resulting in two different samfiles. Now, we examine the Nth read/mate pair in both files (reads and their mates appear as consecutive alignments in the samfile). It is entirely possible that the first alignment of the pair is that of the read, or that both are alignments of its mate, or that one is of the read and one is to its mate.

Therefore, since we intend to examine equivalent alignments in multiple files, we must order read/mate pairs as we get them in (read, mate) order. This function takes a list of "first" pysam aligned_read objects (the aligned_reads parameter) in a pairing as well as the list of "second" ones (the aligned_read_mates parameter), and returns a list of correctly ordered reads, and a list of correctly ordered mates.

.. autofunction:: util.compatibility_dict

compatibility_dict is meant to allow collections.defaultdict to be used on both Python 2.7 and 3.4. It simply overrides the items() function of the dict so that in python 2.7 it uses collections.defaultdict.iteritems(self), and in python 3.4 (where the items() function returns an iterator-like object) it simply uses collections.defaultdict.items(self).
