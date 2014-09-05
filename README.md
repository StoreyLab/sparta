SPARTA
=========

Separate Parental Alleles for Reads from Tangled Alignments (SPARTA)

Divide RNA-Seq reads from pooled runs based on their ancestral genotype. Designed to be used as a post-processing step for Bowtie2, after aligning the pooled reads to both ancestral genomes.

To run unit tests, simply run test_sparta.py

To use sparta.py from the command line:

usage: sparta.py [-h] [-pe [PAIRED_END]] [-n NAMES [NAMES ...]]
                 [-o [OUTPUT_DIR]]
                 [-ss SEPARATED_SAMFILES [SEPARATED_SAMFILES ...]]
                 [-pr [PROCESSES]] [-c [CALCULATE_MISMATCHES]]
                 [-m [MISMATCH_PROB_INPUTFILE]]
                 [-t [TRANSITION_MATRIX_INPUTFILE]] [-ph [PILEUP_HEIGHT]]
                 [-se [SAMPLE_EVERY]] [-g GENOME_PRIORS [GENOME_PRIORS ...]]
                 [-pc [POSTERIOR_CUTOFF]] [-u [UNMAPPED_READ_PROB]]
                 [-i [INSERTION_PROB]] [-d [DELETION_PROB]]
                 [-s [SOFTCLIPPED_PROB]] [-hp [HARDCLIPPED_PROB]] [-q [QUIET]]
                 samfiles [samfiles ...]

SPARTA takes a set of SAM format files that each map the same RNA reads to a
different ancestral (or parental) genome. This program classifies each read to
one of the ancestral alleles or deems it unclassifiable, based on the
assumption that each read belongs to one of the ancestral allele types
provided.

positional arguments:
  samfiles              input samfiles

optional arguments:
  -h, --help            show this help message and exit
  -pe [PAIRED_END], --paired_end [PAIRED_END]
                        set this flag to specify that reads are paired end
                        (default: False)
  -n NAMES [NAMES ...], --names NAMES [NAMES ...]
                        list of nicknames for genomes corresponding for
                        samfile1,samfile2, etc.
  -o [OUTPUT_DIR], --output_dir [OUTPUT_DIR]
                        directory to write output to
  -ss SEPARATED_SAMFILES [SEPARATED_SAMFILES ...], --separated_samfiles SEPARATED_SAMFILES [SEPARATED_SAMFILES ...]
                        list of filenames to write separated (classified) sam
                        outputs. default: outputdir/genome1_separated.sam...
  -pr [PROCESSES], --processes [PROCESSES]
                        number of processes to use for separation step,
                        default = number of CPU cores available
  -c [CALCULATE_MISMATCHES], --calculate_mismatches [CALCULATE_MISMATCHES]
                        set this flag to calculate actual mismatch
                        probabilities for more accurate mapping. WARNING: very
                        slow
  -m [MISMATCH_PROB_INPUTFILE], --mismatch_prob_inputfile [MISMATCH_PROB_INPUTFILE]
                        specify an existing sparta mismatch file (e.g.
                        output/mismatch_prob_info.txt) with mismatch
                        probabilities per quality score for more accurate
                        mapping.
  -t [TRANSITION_MATRIX_INPUTFILE], --transition_matrix_inputfile [TRANSITION_MATRIX_INPUTFILE]
                        specify file with transition matrix in tab-delimited
                        melted format ("A T 0.3" means A to T transition has
                        probability 0.3)
  -ph [PILEUP_HEIGHT], --pileup_height [PILEUP_HEIGHT]
                        if calculate_mismatches is True, specify minimum
                        height of read pileup to consider, default = 20
  -se [SAMPLE_EVERY], --sample_every [SAMPLE_EVERY]
                        if calculate_mismatches is True, specify N such that
                        calculate_mismatch_probs only samples every N reads,
                        default = 10
  -g GENOME_PRIORS [GENOME_PRIORS ...], --genome_priors GENOME_PRIORS [GENOME_PRIORS ...]
                        list of prior probabilities that a read belongs to
                        each genome
  -pc [POSTERIOR_CUTOFF], --posterior_cutoff [POSTERIOR_CUTOFF]
                        lower-bound cutoff for probability that a read belongs
                        to a genome for it to be classified as that genome.
                        default: 0.99
  -u [UNMAPPED_READ_PROB], --unmapped_read_prob [UNMAPPED_READ_PROB]
                        set the (SMALL but NON-ZERO) probability of a read
                        being unmapped (in the SAM) to its genome of origin.
                        default = 0.0001
  -i [INSERTION_PROB], --insertion_prob [INSERTION_PROB]
                        set the (SMALL but NON-ZERO) probability of a read
                        having an inserted base relative to its genome of
                        origin. default = 0.0001
  -d [DELETION_PROB], --deletion_prob [DELETION_PROB]
                        set the (SMALL but NON-ZERO) probability of a read
                        having a deleted base relative to its genome of
                        origin. default = 0.0001
  -s [SOFTCLIPPED_PROB], --softclipped_prob [SOFTCLIPPED_PROB]
                        set the (SMALL but NON-ZERO) probability of a read
                        having a softclipped base relative to its genome of
                        origin. default = 0.0001
  -hp [HARDCLIPPED_PROB], --hardclipped_prob [HARDCLIPPED_PROB]
                        set the (SMALL but NON-ZERO) probability of a read
                        having a hardclipped base relative to its genome of
                        origin. default = 0.0001
  -q [QUIET], --quiet [QUIET]
                        set this flag to supress writing final counts to
                        stdout (default: False)

See the documentation here:
http://storeylab.github.io/sparta
