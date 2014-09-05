.. sparta documentation master file, created by
   sphinx-quickstart2 on Thu Jul 24 16:10:22 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


This is SPARTA's documentation!
==================================

SPARTA (**S**\ eparate **P**\ arental **A**\ lleles for **R**\ eads in **T**\ angled **A**\ lignments) is a bioinformatics tool for python 2.7 and 3.4. It is designed for RNA-seq experiments involving pooled ("tangled") reads sampled from a population of known ancestral allele types (e.g. from the cross of two haploid yeast strains). Following alignment of the reads to the ancestral genomes with Bowtie2, SPARTA serves to separate alignments based on the most likely parental allele type.

Features
**********************************

SPARTA includes the following features:

- Support for an arbitrary number of parental genomes
- Support for the widely accepted SAM format for NGS alignments
- Support for both single read and paired end sequencing experiments
- Fast performance on multicore systems by harnessing all available cores with python multiprocessing
- Optional estimation of mismatch probabilities per quality score from read-pileup to allow more accurate read classification (this is very computationally and memory intensive)
- Produce separate SAM output files for aligned reads that are uniquely classifiable to each parental genome
- Bicompatibile with Python 2.7 and 3.4

Dependencies
**********************************
SPARTA requires the following to be installed:

1. Python 2.7.x or Python 3.4.x
2. Biopython 1.64-2
3. Pysam 0.7.7-1

Earlier versions of python or libraries may work, but are not guaranteed.
For instance, pysam 0.5 is NOT compatible.

Also included is a shell script (sample_analysis_skelly_data.sh) that serves as a sample analysis workflow of real publicly available RNA-seq data (Skelly et al 2011). It downloads the skelly RNA-seq data and runIt depends on the following:

1. a unix-based system
2. `Bowtie2 <https://github.com/BenLangmead/bowtie2>`_ 2.2.3
3. SRA toolkit (fastq-dump tool version 2.3.5)

Older versions of SRA toolkit and bowtie2 MAY work. If you are working on a high-memory system (> ~20 GB) you can add the -c option to the line that runs SPARTA in order to empirically calculate mismatch probabilities.

It is recommended to use GNU/Linux or a Unix-like OS. Sparta was tested on Arch Linux and CentOS. Windows compatibility for SPARTA is **untested** and not guaranteed. For best performance, use python 2.7. Datatype conversions and other modifications that allow python 3.4 to work may result in slightly slower performance.

Contents:
**********************************

.. toctree::
   :maxdepth: 3

   sparta
   calculate_mismatch_probs
   util

.. automodule:: sparta

.. automodule:: calculate_mismatch_probs

.. automodule:: util

Indices and tables
==================================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Github Page with source code:
==================================
https://github.com/StoreyLab/sparta

References:
==================================

Skelly, D. A., Johansson, M., Madeoy, J., Wakefield, J. & Akey, J. M. A powerful and flexible statistical 
framework for testing hypotheses of allele-specific gene expression from RNA-seq data. Genome Res. 
21, 1728–1737 (2011).