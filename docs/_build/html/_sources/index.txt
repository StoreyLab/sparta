.. sparta documentation master file, created by
   sphinx-quickstart2 on Thu Jul 24 16:10:22 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


This is SPARTA's documentation!
==================================

Overview:
**********************************

SPARTA (**S**\ eparate **P**\ arental **A**\ lleles for **R**\ eads in **T**\ angled **A**\ lignments) is a bioinformatics tool for python 2.7 and 3.4. It is designed for RNA-seq experiments involving pooled ("tangled") reads sampled from a population containing two ancestral allele types (e.g. from the cross of two haploid yeast strains). Following alignment of the reads to the ancestral genomes with Bowtie2, SPARTA serves to separate alignments based on the most likely parental allele type.

Features
**********************************

SPARTA includes the following features:

- Support for the widely accepted SAM format for NGS alignments
- Support for both single read and paired end sequencing experiments
- Fast performance on multicore systems by harnessing all available cores with python multiprocessing
- Optional estimation of mismatch probabilities per quality score from read-pileup to allow more accurate read classification
- Produce separate SAM output files for aligned reads that are uniquely classifiable to each genome
- Bicompatibile with Python 2.7 and 3.4

Dependencies
**********************************
SPARTA requires the following to be installed:

1. Python 2.7.x or Python 3.4.x
2. Biopython 1.64-2
3. Pysam 0.7.7-1

Earlier versions of python or libraries may work, but are not guaranteed.
For instance, pysam 0.5 is NOT compatible.

It is recommended to use GNU/Linux or a Unix-like OS. Sparta was tested on Arch Linux and CentOS. Windows compatibility for SPARTA is **untested** and not guaranteed.
The sample analysis scripts run_bowtie_sparta.py and run_sparta_pipeline.py will NOT to work on Windows.

For best performance, use python 2.7. Datatype conversions and other modifications that allow python 3.4 to work may result in slightly slower performance.

Contents:
**********************************


.. toctree::
   :maxdepth: 3

   sparta
   estimate_error_freq
   compatibility

.. automodule:: sparta

.. automodule:: estimate_error_freq

.. automodule:: compatibility

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Github Page with source code:
==================
https://github.com/StoreyLab/sparta
