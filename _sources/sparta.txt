sparta.py
==================================

sparta.py is the main source file for SPARTA. Its main functionality is to take as input a set of SAM format files that map the same set of RNA-seq or DNA-seq reads (i.e. from the same .fastq file) to different genomes via `Bowtie2 <https://github.com/BenLangmead/bowtie2>`_, and sort each read as belonging to one genome or the other. It is assumed that each read belongs to one genome or the other, for instance pooled RNA-seq data from a population of crosses of two haploid yeast strains. sparta will classify each read as belonging uniquely to one of the genomes, or deem it unclassifiable.

.. autofunction:: sparta.sparta

sparta.sparta is the central logic of the entire program (the main function simply parses command line args and passes them directly as the arguments to this function).

sparta.sparta takes the following arguments:

*samfiles*: a list of paths to samfiles that map the same reads to different genomes

*paired_end*: True if the aligned reads are paired end, otherwise False

*genome_names*: A list of names where the Nth name is associated with the Nth sam alignment's genomic reference

*num_processes*: how many child processes (cores) to use while classifying reads

*calculate_mismatches*: whether or not to calculate the probability of mismatch per quality score as well as transition probabilities empirically using :func:`calculate_mismatch_probs.create_mismatch_prob_dict`

*pileup_height*: The minimum number of bases in a "pileup" to calculate mismatches from (if calculate_mismatches).

*sample_every*: If calculate_mismatches, sample every (this many) reads

*genome_priors*: a list of prior probabilities for each genome that any given read belongs to that genome (default: equal probabilities)

*posterior_cutoff*: the minimum posterior probability that a read belongs to a genome to classify it as such

*unmapped_read_prob*: the probability that a read will be unmapped to the genome that it originated from

*insertion_prob*: the probability that a read will have an inserted base compared to the genome that it originated from

*deletion_prob*: the probability that a read will have a deleted base compared to the genome that it originated from

*softclipped_prob*: the probability that a read will have a soft clipped base compared to the genome that it originated from

*hardclipped_prob*: the probability that a read will have a hard clipped base compared to the genome that it originated from

*output_dir*: directory where output files will be written. Will be created if nonexistent.

*separated_samfiles*: a list of file paths to write newly separated samfiles for each genome

*quiet*: if True, do not print total time to stdout

*mismatch_prob_dict_inputfile*: a file that specifies probability of mismatch at a given quality score (e.g. output/mismatch_prob_info.txt). File should be tab-delimited format with the first column being a quality score in ord(ascii) form (integers between 33 - 127), second column should be probabilities of a base mismatch at that quality, and third column should be the number of total number of bases observed in calculating that probability.

*transition_prob_dict_inputfile*: a file that specifies probability of transition from each base to each other base given that a mismatch has occured (e.g. output/transition_prob_info.txt). The file should be tab-delimited and column 1 should be the reference base, column 2 should be the base being transitioned to, and column 3 should be the probability of transition from the base in column 1 to the base in column 2.

A brief description of what the sparta.sparta function does:

If the names of inputfiles for transition and mismatch probability information is provided, then these are parsed. Alternatively, they can be calculated new using calculate_mismatch_probs.create_mismatch_prob_dict or skipped altogether (in which case mismatch probability is calculated directly from quality score and transition probability is fixed as 1/3 for each transition)

If multiple processes are specified then a pool of worker processes is created, set to process the samfiles in interleaved fashion, and then recombined. Otherwise, _worker_procedure is called directly from the parent process. Finally, all results are printed to the output directory.

**The remaining classes and methods are designed to be relatively 'private' to sparta (despite the distinction of public vs. private not existing in Python) and should be used to understand how the program works or for forking/modifying/extending the code, and are NOT recommended to be called by a user desiring intended SPARTA functionality.**

.. autoclass:: sparta.multimapped_read_separator

 The multimapped_read_separator class encapsulates the procedure of classifying a series of reads to their respective parental genomes.

 It contains lists that store the results of the sort: multimapped_read_separator.sort_fates, multimapped_read_separator.mismatch_logs, and multimapped_read_separator.posterior_logs.

	.. automethod:: sparta.multimapped_read_separator.log

	This method logs the results of one sorting decision into the log lists. mismatch_list is indexed by genome number (genome of first samfile is 0, genome of second samfile is 1, etc.), that records the number of mismatches to each genome. posterior_list is also indexed by genome number, and records the posterior probability that a read belongs to each genome (as calculated by sparta.multimapped_read_separator.untangle_samfiles). sort_fate simply contains the sorting decision of each read: 0 means unsorted whereas 1-N mean that the read was classified as that genome (same genome index scheme as earlier except 1-indexed).

	.. automethod:: sparta.multimapped_read_separator.aligned_read_prob

	aligned_read_prob computes the prior probability that a read was generated given that it was produced by the genome it is aligned to. For this reason, its only argument is a pysam alignedread object.

	A set probability per base for differences from the reference that are not mismatches (0.0001 by default), such as insertions/deletions/etc. Then, the product of the probability of each aligned base is computed. A match has probability (1 - Pr(mismatch given quality)) while a mismatch has probability (Pr(mismatch) * Pr(specific transition given that there is a mismatch)).

	.. automethod:: sparta.multimapped_read_separator.count_mismatches

	This method simply counts the total number of base mismatches present (from the alignedread's MD string) for a list of alignedreads. This is stored and is interesting for later analyses, although it does not directly play a role in classification.

	.. automethod:: sparta.multimapped_read_separator.bayes_classify

	This method applies bayes rule with the law of total probability to a set of prior probabilities (one for each aligned genome). It then determines if any of the probabilities meet the posterior probability cutoff (--posterior_cutoff argument). It returns a string describing the classification decision, the list of posterior probabilities, and the sort fate as an integer.

	.. automethod:: sparta.multimapped_read_separator.untangle_mappings

	untangle_mappings takes a list alignedreads and an optional list of mates and returns a string describing the classification decision, the list of posterior probabilities of each possible genome, and the integer sort fate (0 means unsorted and 1-N are the 1-indexed input genomes). For unmapped reads, the method applies a very small nonzero probability that the read was generated by that genome. For paired end reads, the prior probabilities of the read and the mate are multiplied together before applying bayes rule (essentially treating them as one long read). 

	.. automethod:: sparta.multimapped_read_separator.untangle_samfiles

	untangle_samfiles is one further step up the heirarchy from untangle_mappings; it opens all of the samfiles and iterates over them in parallel and calls untangle_mappings to obtain a classification result. 

	alignedreads are processed at an interval appropriate for the number of processes being used, so that the results (which do not depend on one another by definition) can be recombined into the results for the entire files. For example, if there are 4 processes being distributed across 4 cores, then core 1 might process lines 0,4,8... in the file while core 2 processes lines 1,5,9..., core 3 processes 2,6,10... and on. In this example, the num_processes argument would be 4 and each process would have a different interleave_ix (0-3).

	For paired-end reads, there is the additional task of ensuring that reads and mates are ordered properly before calling untangle_samfiles, and this is handled with the function util.fix_read_mate_order.

	.. automethod:: sparta.multimapped_read_separator.write_samfiles

	This method makes a final pass through the samfiles and writes new samfiles for each genome (where each samfile contains only the reads classified to each genome).

	It is required that the length of sort_fates equals the number of fastq reads (in other words, that results have already been recombined if multiprocessing was used).

	The single argument to the method is a list of pathnames where the new samfiles will be written.


.. autofunction:: sparta._worker_procedure

This function encapsulates the procedure that a worker process executes when multiprocessing is used.

A new multimapped_read_separator is created and used to separate reads in a samfile according to a specific interleave_ix. 

The separator itself is returned.

.. autofunction:: sparta.merge_separators

Merge a list of multimapped_read_separators (specifically, the mismatch/posterior_prob/sort_fate logs) and return the combined separator.

Items in the returned separator are exactly like they would be if the code were run linearly without multiprocessing.
