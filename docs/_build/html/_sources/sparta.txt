sparta.py
==================================

sparta.py is the main source file for SPARTA. Its main functionality is to take as input two SAM format files that map the same set of RNA-seq reads (i.e. from the same .fastq file) to different genomes, and sort each read as belonging to one genome or the other. It is assumed that each read belongs to one genome or the other, for instance pooled RNA-seq data from a population of crosses of two haploid yeast strains.

.. autofunction:: sparta.sparta

sort_samfiles takes the pathnames of two SAM format files, each that map

.. autoclass:: sparta.multimapped_read_separator

.. automethod:: sparta.multimapped_read_separator.aligned_read_prob

.. autofunction:: sparta.merge_separators

