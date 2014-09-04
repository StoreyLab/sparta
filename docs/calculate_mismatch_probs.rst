calculate_mismatch_probs.py
==================================

**The following classes and methods are designed to be relatively ‘private’ to sparta (despite the distinction of public vs. private not existing in Python) and should be used to understand how the program works or for forking/modifying/extending the code, and are NOT recommended to be called by a user desiring intended SPARTA functionality.**

The functionality of calculate_mismatch_probs is to take a set of samfiles that align the same next-gen sequencing reads to different genomes (preferably aligned with Bowtie2) and calculate the probability that bases of a given quality score will mismatch against the consensus. calculate_mismatch_probs inspects reads that map to all provided genomes. It finds base positions where the read agrees with every genome. Assuming that a likely consensus is observed (greater than 75% of "piled up" reads agree), any aberrant bases in the pileup are labeled as mismatches. Then, the proportion of matches vs. mismatches is tallied for each possible base pair quality score, and from this, the probability of random mismatch at each quality score is obtained. 

.. autofunction:: calculate_mismatch_probs.create_genome_seq

The purpose of create_genome_seq is to take in an alignedread object (*aligned*) and return a representation of the genomic sequence that it aligns to. It is assumed that the alignedread object does not have insertions/deletions/clippings/etc.

If the sequence is reversed, the reverse complement of the genomic sequence will be returned.

In summary, this function simply edits the original read sequence with the mismatches specified by the sam alignment's MD field.

.. autofunction:: calculate_mismatch_probs.add_to_pileup_dict

add_to_pileup_dict logs information about alignments of reads to each genome in a nested dictionary structure. The goal is to record, for a set of positions (one single-base genomic position from each genome), how many times each base (A, C, G, T) is observed aligning to ALL of those positions in an alignment for a given base quality. Saving this information will allow determination of where likely consensuses are (where do all genomes agree, and bases from aligned reads agree?)

If the read did not align to each genome, or if any of the alignments include sequence deviations from the reference other than single-base mismatches, the read is passed over. Disallowing insertions/deletions/etc allows for the code to be much simpler, and is allowable because only a sufficiently large subset of reads needs to be sampled anyway. Then, the original genomic sequence is re-created for each alignment (for ease of recording genomic base identities). The length of the read is iterated, and for each position, the coordinates in each genome and the genomic base identity are recorded in the dictionary, as well as the base quality and base identity of read. The dictionary stores data in the following manner:

level 1 of dictionary: tuple of tuples of genomic coordinates and genomic base identities for each genome
level 2 of dictionary: identity of base on the read
level 3 of dictionary: quality score of read base

the value stored in the dictionary of counts for number of times a given combination of the above values is observed.

.. autofunction:: calculate_mismatch_probs.create_mismatch_prob_dict

This function takes a set of samfiles that align the same reads to different genomes, and returns a dictionary mapping base quality to likelihood of mismatch (as well as other supplementary information).

It first accumulates counts of bases seen in "pileups" (in quotes to distinquish from pysam's pileup feature) using add_to_pileup_dicts. Then, for each "pileup" (referring to any given set of genomic coordinates with counts of aligned bases to that set of coordinates), it requires that the number of total observations reaches a cutoff (default=20).

It determines if, for a given "pileup", a consensus is observed (defined by default as having greater than 75% of bases in the pileup agree on the same base as each genome), in which case, the number of matches and mismatches to that consensus is tallied for each quality score. The transitions between different bases are also tallied to form a separate transition probability table.

Finally, these tallied values are turned into probabilities and returned. Files titled "mismatch_prob_info.txt" and "transition_prob_info.txt" are also written to the output directory with the relevant information.