sparta.py
==================================

sparta.py is the main source file for SPARTA. Its main functionality is to take as input two SAM format files that map the same set of RNA-seq reads (i.e. from the same .fastq file) to different genomes, and sort each read as belonging to one genome or the other. It is assumed that each read belongs to one genome or the other, for instance pooled RNA-seq data from a population of crosses of two haploid yeast strains.

.. autofunction:: sparta.parseargs

.. autoclass:: sparta.multimapped_read_separator

.. automethod:: sparta.multimapped_read_separator.log

.. automethod:: sparta.multimapped_read_separator.aligned_read_prob

.. automethod:: sparta.multimapped_read_separator.count_mismatches

.. automethod:: sparta.multimapped_read_separator.bayes_classify

.. automethod:: sparta.multimapped_read_separator.untangle_mappings

.. automethod:: sparta.multimapped_read_separator.untangle_samfiles

.. automethod:: sparta.multimapped_read_separator.write_samfiles

.. automethod:: sparta.multimapped_read_separator._worker_procedure

.. autofunction:: sparta.merge_separators

.. autofunction:: sparta.sparta

