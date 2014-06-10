#!/usr/bin/env python
"""

    rnaseq_pipeline.py
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]
                        [--jobs]
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]

"""
import sys
import os
import re
import time
import urllib2
import random
import subprocess
import itertools
import collections
import gzip

from Bio import Seq, SeqIO, SeqRecord


module_name = "snp_pipeline"
JOB_NAME = "RUN_111213OB"
#JOB_NAME = "Akey"

if not os.path.exists(JOB_NAME):
    os.mkdir(JOB_NAME)
os.chdir(JOB_NAME)

# input files
GENOME_DATA = "genome_data"
GENOME_FASTA = os.path.join(GENOME_DATA,
                            "S288C_reference_sequence_R64-1-1_20110203.fsa")
GENOME_GFF = os.path.join(GENOME_DATA,
                        "saccharomyces_cerevisiae_R64-1-1_20110208.gff")


# subsampling proportion parameters: on a log10 scale
SUBSAMPLING_START = -3
SUBSAMPLING_END = 0
SUBSAMPLING_STEP = .05


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %progs [options]")



    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", default=0,
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")




    #
    #   pipeline
    #
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        default=1,
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", default=False,
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")

    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true", default=False,
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = [ ]

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have b een defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)

        if not len(missing_options):
            return

        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
    check_mandatory_options (options, mandatory_options, helpstr)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions

def levenshtein_distance(s1, s2):
    """
    Python version of Levenshtein distance for compatability. The third-party
    library Levenshtein is much faster and recommended. Taken from recipe at:
    http://code.activestate.com/recipes/576874-levenshtein-distance/
    """
    l1 = len(s1)
    l2 = len(s2)

    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in xrange(l2 + 1):
        matrix[zz] = range(zz, zz + l1 + 1)
    for zz in range(0, l2):
        for sz in range(0, l1):
            if s1[sz] == s2[zz]:
                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
                                    matrix[zz][sz + 1] + 1, matrix[zz][sz])
            else:
                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
                                    matrix[zz][sz + 1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]


try:
    import Levenshtein
    levenshtein_distance = Levenshtein.distance
except ImportError:
    pass


class BarcodeIndex:
    """Represents a set of indices, with the ability to find the closest one"""
    def __init__(self, index_file, max_distance=3):
        self.cache = {}
        with open(index_file) as inf:
            for l in inf:
                if l.startswith("#") or l == "":
                    continue
                g, barcode = l[:-1].split("\t")
                self.cache[barcode] = g
        self.barcodes = self.cache.keys()
        self.groups = self.cache.values()
        self.max_distance = max_distance
        self.hits = 0
        self.total = 0

    def find_barcode(self, barcode):
        """
        match barcode and return the group it's supposed to be in. If
        there is none within the Levenshtein distance given, or if there
        is a tie, return None
        """
        self.total += 1
        # if there's an exact match, return that
        exact = self.cache.get(barcode)
        if exact is not None or self.max_distance == 0:
            self.hits += 1
            return exact

        # find the Levenshtein distance to each
        distances = [levenshtein_distance(barcode, b) for b in self.barcodes]
        best = min(distances)
        # check if there's a tie or the distance is too great:
        if best > self.max_distance or distances.count(best) > 1:
            self.cache[barcode] = None
            return None
        # otherwise, return the best one, after caching it for future use
        ret = self.groups[distances.index(best)]
        self.cache[barcode] = ret
        return ret



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#       Put pipeline code here

GENOME_BUILD = "S288C_reference_genome_R64-1-1_20110203"
GENOME_FASTA_ORIGINAL = os.path.join(GENOME_BUILD,
                            "S288C_reference_sequence_R64-1-1_20110203.fsa")
GENOME_GFF_ORIGINAL = os.path.join(GENOME_BUILD,
                        "saccharomyces_cerevisiae_R64-1-1_20110208.gff")


@files([], [GENOME_FASTA_ORIGINAL, GENOME_GFF_ORIGINAL])
def wget_genome(infiles, outfile):
    """
    Download the R64 reference genome from SGD using wget and unzip it
    """
    url = ("http://downloads.yeastgenome.org/sequence/S288C_reference/" +
           "genome_releases/%s.tgz" % (GENOME_BUILD))
    os.system("wget %s" % (url))
    os.system("tar -xvf %s.tgz" % (GENOME_BUILD))


@follows(mkdir("bowtie_index"), wget_genome)
@files(GENOME_FASTA_ORIGINAL, os.path.join("bowtie_index", "S288C.fsa"))
def convert_fasta(infile, outfile):
    """
    convert default fasta file to the appropriate format, named by chromosome
    """
    with open(infile) as inf:
        with open(outfile, "w") as outf:
            for l in inf:
                if l.startswith(">"):
                    if "chromosome" in l:
                        l = ">chr%s\n" % re.search("chromosome=(.*?)\]",
                                                l).groups()[0]
                outf.write(l)

@follows(mkdir("bowtie_index"), wget_genome)
@files("RM11_1A/assembly/genome.fa", os.path.join("bowtie_index", "RM11.fsa"))
def convert_fastaRM(infile, outfile):
    """
    convert default fasta file to the appropriate format, named by chromosome
	drop scplasm and chrm
    """
    ignore=False
    with open(infile) as inf:
        with open(outfile, "w") as outf:
            for l in inf:
                if l.startswith(">"):
		    if "scplasm" in l or "m" in l:
			ignore=True
                    else:
			ignore=False
			num = int(l[len(l)-3:])
			l = ">chr"+str(num)+"\n"
                        #l = ">chr%s\n" % re.search("chromosome=(.*?)\]",
                                                #l).groups()[0]
		if(not ignore):               
		    outf.write(l)

@follows(wget_genome)
@files(GENOME_GFF_ORIGINAL, "S288C.gff")
def remove_fasta_gff(infile, outfile):
    """
    Remove the fasta sequence from the end of the .gff file, or htseq-count
    can't use it. Truncate the file
    """
    with open(infile) as inf:
        with open(outfile, "w") as outf:
            for l in inf:
                if l.startswith("##FASTA"):
                    break
                outf.write(l)


@follows(convert_fasta)
@files(os.path.join("bowtie_index", "S288C.fsa"),
         os.path.join("bowtie_index", "S288C.1.bt2"),
	 os.path.join("bowtie_index", "RM11.fsa"),
         os.path.join("bowtie_index", "RM11.1.bt2"))
def build_bowtie(infile_by, outfile_by, infile_rm, outfile_rm):
    out_base_by = outfile_by.split(".1.bt2")[0]
    out_base_rm = outfile_rm.split(".1.bt2")[0]
    command = "bowtie2-build %s %s" % (infile_by, out_base_by)
    command2 = "bowtie2-build %s %s" % (infile_rm, out_base_rm)
    os.system(command)
    os.system(command2)


@follows(build_bowtie, mkdir("truncated"))
@transform(os.path.join("split_fastq", "[EG]*"),
           regex(".*\/([^\/]*)$"),
           os.path.join("truncated", r"\1"))
def truncate(infile, outfile):
    """
    truncate each line to only 67 characters, which can be used to study the
    effect of read length on conclusions
    """
    if PARALLEL:
        cmd = "qsub -cwd -o %s cut -c1-67 %s" % (outfile, infile)
    else:
        cmd = "cut -c1-67 %s > %s" % (infile, outfile)
    print cmd
    os.system(cmd)


#INPUT_DATA = "/Genomics/grid/users/dgrtwo/Storey/ase/data"
INPUT_DATA = "/Genomics/grid/users/emilysn/ase/analysis/RNA-Seq/Akey"


@files((os.path.join(INPUT_DATA, "ASE_lane_design.txt"),
        os.path.join(INPUT_DATA, "indices.fasta")),
        ["lane0_index.txt", "lane1_index.txt"])
def lane_indices(infiles, outfiles):
    design_infile, indices_infile = infiles
    outfs = [open(o, "w") for o in outfiles]

    indices = SeqIO.to_dict(SeqIO.parse(indices_infile, "fasta"))

    with open(design_infile) as inf:
        header = next(inf)
        for l in inf:
            spl = l.strip().split("\t")
            s = "\t".join([spl[4], str(indices[spl[-3]].seq)])
            outfs[int(spl[-1])].write(s + "\n")
    for o in outfs:
        o.close()


LANE0_FASTQ_GZ = os.path.join(INPUT_DATA, "ase-zero-for-148-cycles-h7bh1adxx_1_read_1_passed_filter.fastq.gz")
LANE0_BARCODES = os.path.join(INPUT_DATA, "ase-zero-for-148-cycles-h7bh1adxx_1_read_2_index_read_passed_filter.fastq.gz")

LANE1_FASTQ_GZ = os.path.join(INPUT_DATA, "ase-one-for-148-cycles-h7bh1adxx_2_read_1_passed_filter.fastq.gz")
LANE1_BARCODES = os.path.join(INPUT_DATA, "ase-one-for-148-cycles-h7bh1adxx_2_read_2_index_read_passed_filter.fastq.gz")


@follows(lane_indices, mkdir("split_fastq"))
@files([[[LANE0_FASTQ_GZ, LANE0_BARCODES, "lane0_index.txt"],
         os.path.join("split_fastq", "G1alpha0.fastq")],
        [[LANE1_FASTQ_GZ, LANE1_BARCODES, "lane1_index.txt"],
         os.path.join("split_fastq", "E1a0.fastq")]])
def split_reads(infiles, outfile):
    """
    Given a gzipped fastq file of reads, a gzipped fastq file with
    corresponding barcodes, and a file with the barcode index. Create one fastq
    file for each sample based on the closest matching barcode
    """
    read_file, barcode_file, index_file = infiles
    if "lane1" not in index_file:
        return

    index = BarcodeIndex(index_file, max_distance=2)

    counts = collections.defaultdict(int)

    lane_name = os.path.splitext(index_file)[0]

    # one output for each group
    unassigned = "Unassigned_%s" % lane_name
    infolder = os.path.split(outfile)[0]
    outfs = dict([(g, open(os.path.join(infolder, g + ".fastq"), "w"))
                    for g in index.groups + [unassigned]])

    from itertools import izip
    starting_time = time.time()
    with gzip.GzipFile(read_file) as read_inf:
        with gzip.GzipFile(barcode_file) as barcode_inf:
            for i, (r, b) in enumerate(izip(izip(*[read_inf]*4),
                                            izip(*[barcode_inf]*4))):
                group = index.find_barcode(str(b[1][:-1]))
                group = group if group != None else unassigned
                outfs[group].write("".join(r))
                counts[group] += 1

                if i % 10000 == 0 and i > 0:
                    print i, (i) / (time.time() - starting_time),
                    print float(index.hits) / index.total

    for k, v in sorted(counts.items(), key=lambda t: t[0] == unassigned):
        print k, v

    for o in outfs.values():
        o.close()


PARALLEL = True


@follows(build_bowtie, split_reads, mkdir("mapped"), mkdir("reports"))
@transform(os.path.join("split_fastq", "[EG]*"),
           regex(".*\/([^\/]*).fastq"),
           os.path.join("mapped", r"\1.bam"), r"\1")
def map_reads(infile, outfile, base):
    index = os.path.join("bowtie_index", "S288C")

    command = "bowtie2 %s %s | samtools view -bS - > %s" % (
                    index, infile, outfile)

    if PARALLEL:
        job_file = os.path.join("reports", "bowtie2_%s.txt" % base)
        command = "echo \"%s\" | qsub -cwd -e %s -b n" % (command, job_file)
    print command
    os.system(command)

@follows(build_bowtie, split_reads, mkdir("mapped"), mkdir("reports"), mkdir("mappedRM"))
@transform(os.path.join("split_fastq", "[EG]*"),
           regex(".*\/([^\/]*).fastq"),
           os.path.join("mappedRM", r"\1.bam"), r"\1")
def map_readsRM(infile, outfile, base):
    index = os.path.join("bowtie_index", "RM11")

    command = "bowtie2 %s %s | samtools view -bS - > %s" % (
                    index, infile, outfile)

    if PARALLEL:
        job_file = os.path.join("reports", "bowtie2_%s.txt" % base)
        command = "echo \"%s\" | qsub -cwd -e %s -b n" % (command, job_file)
    print command
    os.system(command)


@follows(map_reads, map_readsRM, mkdir("comparisons"))
@collate(["mapped/*", "mappedRM/*"], regex("mapped.*\/(.*).bam"),
	 r"comparisons/\1.txt")
def compare_mappings(infiles, outfile):
	byfile, rmfile = infiles
	import pysam
	import itertools
	import collections
	name = byfile.split(".")[0]
	print name
	snpname=name+".snps"
	bysam = pysam.Samfile(byfile)
	rmsam = pysam.Samfile(rmfile)
	nummap = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
	commonmap = collections.defaultdict((collections.Counter))

	def getrelpos(spl):
		total = 0
		rel_pos = []
		for n, b in spl:
			total+=int(n)
			rel_pos.append((total, b))
			if not b.startswith("^"):
				total+=1
		return rel_pos

	match1 = 0
	nomatch = 0
	matchboth = 0	
	ident_errors = 0
	weirdness = 0
	both_snps = 0
	difflen = 0
	diffbase = 0
	total = 0
	suspects = 0

	MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)")

	for total, (byread, rmread) in enumerate(itertools.izip(bysam, rmsam)):
		assert byread.qname == rmread.qname
		total+=1
		#print byread.is_unmapped, rmread.is_unmapped
		if byread.is_unmapped and rmread.is_unmapped:
			#this read is probably junk
			nomatch+=1
			continue
		elif byread.is_unmapped != rmread.is_unmapped:
			#maps to one but not the other; either junk or
			#an unshared gene bw BY and RM
			match1+=1
			continue
		else:
			#maps to both; these are actually interesting
			matchboth+=1
			#pull out chromosome numbers and positions
			chrby = bysam.getrname(byread.rname)
			chrrm = rmsam.getrname(rmread.rname)
			posby = byread.pos
			posrm = rmread.pos
			#put chr/pos in dictionary
			#nummap[(chrby, chrrm)] = nummap[(chrby, chrrm)](posby, posrm)
			#get errors
			errby = byread.opt("MD")
			errrm = rmread.opt("MD")

			#get quality score
			qual = byread.qual # will be the same for both

			splby = re.findall(MD_REGEX, errby)
			splrm = re.findall(MD_REGEX, errrm)

			#if the errors are the same I dont care--not a SNP!
			if errby==errrm:
				ident_errors+=1
				commonmap[(chrby, posby)][(chrrm, posrm)] += 1
				continue

			if len(byread.positions)!=len(rmread.positions):
				if (len(splby)==0 and len(splrm)==0):
					ident_errors+=1
					commonmap[(chrby, posby)][(chrrm, posrm)] += 1
				else:
					difflen+=1
				continue

			rel_BY = getrelpos(splby)
			rel_RM = getrelpos(splrm)

			#print "ORIGINAL ERROR CODE: "+str(errby), str(errrm)
			#print "RELATIVE POS: "+str(rel_BY), str(rel_RM)

			inboth = set(n for n, b in rel_BY).intersection(n for n, b in rel_RM)

			full_BY = [(byread.positions[n], rmread.positions[n], ord(qual[n]), b) for n, b in rel_BY if not n in inboth]
			full_RM = [(byread.positions[n], rmread.positions[n], ord(qual[n]), b) for n, b in rel_RM if not n in inboth]

			#print "UNIQUE FULL POSITIONS: "+str(full_BY), str(full_RM)

			errorsby = len(full_BY)
			errorsrm = len(full_RM)

			if not (errorsby==0 or errorsrm==0): #read has both BY and RM SNPs
				both_snps+=1
				continue

			if (full_BY==[] and full_RM==[]): #this is the case where both reads differ at the same position, by a different base
				diffbase+=1
				continue

			#dont double count SNPs! one SNP per read should be counted
			#just take the first one

			if(errorsby==0): #SNPs map to RM
				full_RM =  [full_RM[0]]
			if(errorsrm==0):
				full_BY = [full_BY[0]]

			for posby, posrm, q, b in full_BY:
				suspects+=1
				nummap[(chrby, posby)][(chrrm, posrm)][(q, b, "RM")] += 1
			for posby, posrm, q, b in full_RM:
				suspects+=1
				nummap[(chrby, posby)][(chrrm, posrm)][(q, b, "BY")] += 1

			#stick in dict
			#if this chr combination is not in nummap yet:			
			#if (chrby, posby) in nummap and (chrrm, posrm) in nummap[(chrby, posby)]:
			#	nummap[(chrby, posby)][(chrrm, posrm)].append((full_BY, full_RM)) 
			#chr in but not positions
			#elif (chrby, posby) in nummap and not (chrrm, posrm) in nummap[(chrby, posby)]:
			#	nummap[(chrby, posby)][(chrrm, posrm)] = [(full_BY, full_RM)]
			#if one or both is missing:
			#else:
			#	nummap[(chrby, posby)] = {(chrrm, posrm) : [(full_BY, full_RM)]}



	print "Read match only 1 genome: "+str(float(match1) / total)
	print "No match in either genome: "+str(float(nomatch) / total)
	print "Matched both genomes: "+str(float(matchboth) / total)
	print "Read matched both, seq error identical (no SNP): "+str(float(ident_errors)/total)
	print "Read has SNPs that appear to be from both BY and RM: "+str(float(both_snps)/total)
	#print "Total: "+str(total)
	#print "Different base: "+str(diffbase)
	#print "Different length reads: "+str(difflen)

	#print commonmap

	outfile = open(snpname, 'w')
	header = "CHRBY\tPOSBY\tCHRRM\tPOSRM\tMATCH\tBASE\tQUAL\tFREQ\n"
	outfile.write(header)
	bykeys = nummap.keys()
	for k in bykeys:
		rmkeys = nummap[k].keys()
		#print k
		for j in rmkeys:
			#print j
			for i in nummap[k][j]:
				#print k, j, i
				s = str(k[0])+"\t"+str(k[1])+"\t"+str(j[0])+"\t"
				s+= str(j[1])+"\t"+str(i[2])+"\t"+str(i[1])+"\t"+str(i[0])+"\t"
				s+=str(nummap[k][j][i]) + "\n"
				#print s
				outfile.write(s)  

	commonname = name+".common"	
	outfile = open(commonname, 'w')
	header = "CHRBY\tPOSBY\tCHRRM\tPOSRM\tFREQ\n"
	outfile.write(header)
	bykeys = commonmap.keys()
	for k in bykeys:
		rmkeys = commonmap[k].keys()
		#print k
		for j in rmkeys:
			s=str(k[0])+"\t"+str(k[1])+"\t"
			s+=str(j[0])+"\t"+str(j[1])+"\t"
			s+=str(commonmap[k][j])+"\n"
			outfile.write(s)
  
				


@follows(mkdir("counts"), remove_fasta_gff)
@transform(os.path.join("mapped", "*"),
           regex(os.path.join("mapped", "(.*).bam")),
           os.path.join("counts", r"\1.txt"))
def match_reads(infile, outfile):
    command = ("samtools view %s | htseq-count -i ID -s no -t gene " +
                    "- %s > %s") % (infile, "S288C.gff", outfile)
    if PARALLEL:
        command = "echo \"%s\" | qsub -cwd -b n" % command

    print command
    os.system(command)

@follows(mkdir("countsRM"), remove_fasta_gff)
@transform(os.path.join("mappedRM", "*"),
           regex(os.path.join("mappedRM", "(.*).bam")),
           os.path.join("countsRM", r"\1.txt"))
def match_readsRM(infile, outfile):
    command = ("samtools view %s | htseq-count -i ID -s no -t gene " +
                    "- %s > %s") % (infile, "RM11.gff", outfile)
    if PARALLEL:
        command = "echo \"%s\" | qsub -cwd -b n" % command

    print command
    os.system(command)


@follows(match_reads)
@merge(os.path.join("counts", "*"), "count_matrix.txt")
def combine_counts(infiles, outfile):
    """Combine all htseq-count files into one matrix"""
    d = collections.defaultdict(dict)
    for f in infiles:
        name = os.path.splitext(os.path.split(f)[1])[0]
        with open(f) as inf:
            for l in inf:
                gene, count = l[:-1].split("\t")
                d[name][gene] = count
    print d.values()
    genes = sorted(d.values()[0].keys())

    samples = sorted(d.keys())

    with open(outfile, "w") as outf:
        outf.write("\t".join(["Gene"] + samples) + "\n")
        for g in genes:
            outf.write("\t".join([g] + [d[s][g] for s in samples])
                        + "\n")


@follows(combine_counts)
@files("count_matrix.txt", "pooled_subsets.Rdata")
def pool_subsets(infile, outfile):
    """
    create an .Rdata file with a list of each subdesign of the matrix (Full,
    1 lane/1 prep, etc) that will be used for subsampling
    """
    # this is done by an R script
    os.system("../pool_subsets.R %s %s" % (infile, outfile))


SUBSET_NAMES = ["Full", "preparationA", "preparationB",
                "lane1", "lane2", "preparationA:lane1", "preparationA:lane2",
                "preparationB:lane1", "preparationB:lane2"]



@follows(pool_subsets, mkdir("subSeq_subsamples"))
@files([("pooled_subsets.Rdata",
         os.path.join("subSeq_subsamples", n + "_subsampled.Rdata"), n)
            for n in SUBSET_NAMES])
def subsample_subsets(infile, outfile, base):
    """
    For each of the subsamples in the pool_subsets file, perform subsampling
    using my R subSeq package
    """
    cmd = "../subsample_subSeq.R %s %s %s %f %f %f" % (infile, outfile,
                        base, SUBSAMPLING_START, SUBSAMPLING_END,
                        SUBSAMPLING_STEP)
    if PARALLEL:
        cmd = "qsub -cwd %s" % (cmd)
    print cmd
    os.system(cmd)


@follows(subsample_subsets)
@merge(os.path.join("subSeq_subsamples", "*"), "subsample_summary.Rdata")
def summarize_subsamples(infiles, outfile):
    """
    Combine the subsampling runs in the folder into a single table summarizing
    by depth (which is used as the input to the knitr manuscript)
    """
    infolder = os.path.split(infiles[0])[0]
    oracle = [f for f in infiles if "Full" in f][0]
    cmd = "../summarize_subsamples.R %s %s %s" % (infolder, oracle, outfile)
    print cmd
    os.system(cmd)


@follows(subsample_subsets, mkdir("GSEA_subsamples"))
@transform(os.path.join("subSeq_subsamples", "*"),
           regex(os.path.join("subSeq_subsamples", "(.*)\.Rdata")),
           os.path.join("GSEA_subsamples", r"\1.Rdata"))
def GSEA_subsamples(infile, outfile):
    """
    Perform gene set enrichment analysis using a Wilcoxon test on each
    subsample
    """
    os.system("qsub -cwd ../GSEA_subsamples.R %s %s" % (infile, outfile))


def summarize_GSEA(infile, outfile):
    """
    Summarize results from gene set enrichment analysis
    """
    os.system("qsub -cwd ../GSEA_subsamples.R %s %s" % (infile, outfile))


@merge(os.path.join("split_fastq", "[EG]*.fastq"), "length_summary.txt")
def analyze_quality(infiles, outfile):
    with open(outfile, "w") as outf:
        for f in infiles:
            from itertools import izip
            real_length = [0] * 142
            with open(f) as inf:
                # read them in four at a time, the fast way:
                for i, read in enumerate(izip(*[inf]*4)):
                    real_length[len(read[3][:-1].rstrip("#"))] += 1
                    if i > 1000000:
                        break
            sample_name = os.path.basename(f).split(".")[0]
            outf.write(sample_name + "\t" + "\t".join(map(str, real_length))
                            + "\n")


@follows(mkdir("sorted_mapped"), map_reads)
@transform("mapped/*", regex("mapped\/(.*)"), r"sorted_mapped/\1.bam")
def sort_reads(infile, outfile):
    cmd = "samtools sort %s %s" % (infile, outfile)
    if PARALLEL:
        cmd = "qsub -cwd " + cmd
    print cmd
    os.system(cmd)


@follows(sort_reads)
@transform("sorted_mapped/*", suffix(".bam"), ".bam.bai")
def index_reads(infile, outfile):
    """create .bai index of sorted read mapping files"""
    cmd = "samtools index %s" % infile
    if PARALLEL:
        cmd = "qsub -cwd " + cmd
    print cmd
    os.system(cmd)


@follows(index_reads, mkdir("selected_fastq"))
@collate("sorted_mapped/*.bam", regex("sorted_mapped\/([EG][12]).*"),
         r"selected_fastq/\1.fastq")
def select_region(infiles, outfile):
    """Select reads from a region of the genome and write them out"""
    with open(outfile, "w") as outf:
        for f in infiles:
            cmd = ["samtools", "view", f, "chrXII:270000-300000"]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            for l in p.stdout:
                spl = l.split("\t")
                i, dir = spl[:2]
                seq, phred = spl[9:11]
                if dir == "16":
                    seq = str(Seq.Seq(seq).reverse_complement())
                s = "@%s\n%s\n+\n%s\n" % (i, seq, phred)
                outf.write(s)


@follows(select_region, mkdir("subset_reads"))
@transform("selected_fastq/*", regex("selected_fastq\/(.*)"),
           r"subset_reads/\1")
def subset_reads(infile, outfile):
    """Select 1/10 of these reads, and shuffle them (they're sorted)"""
    reads = []
    for r in SeqIO.parse(infile, "fastq"):
        if random.random() < .1:
            reads.append(r)
    random.shuffle(reads)
    SeqIO.write(reads, outfile, "fastq")



### extract SNP reads ###

@files("../snps/SNPcatalogII.txt", ("../snps/SNP_seqs_BY.fastq",
                                 "../snps/SNP_seqs_RM.fastq"))
def write_fastq_SNPs(infile, outfiles):
    inf = open(infile)
    BY_outf = open(outfiles[0], "w")
    RM_outf = open(outfiles[1], "w")
    for l in inf:
        if len(l.split("\t")) == 5:
            gene, ID, change, BY, RM = l[:-1].split("\t")
            #print BY, RM
            diffs = [a != b for a, b in zip(BY, RM)]
            #print diffs.count(True)
            if diffs.count(True) != 1:
                continue

            if diffs.index(True) != 20:
                end = "I" * len(BY)
                BY_outf.write("@%s|%s|%s|BY\n%s\n+\n%s\n" %
                                    (gene, ID, change, BY, end))
                RM_outf.write("@%s|%s|%s|RM\n%s\n+\n%s\n" %
                                    (gene, ID, change, RM, end))
    BY_outf.close()
    RM_outf.close()


@follows(write_fastq_SNPs)
@files("../snps/SNP_seqs_BY.fastq", "../snps/SNP_BY.sam")
def bowtie_SNPs(infile, outfile):
    command = ["bowtie2", "genome_data/BY", infile, "-S", outfile]
    subprocess.call(command)


@follows(bowtie_SNPs)
@files("../snps/SNP_BY.sam", "../snps/SNP_positions.txt")
def SNP_positions(infile, outfile):
    import pysam
    outf = open(outfile, "w")
    samfile = pysam.Samfile(infile)
    for r in samfile:
        if r.is_unmapped:
            continue
        snp_pos = r.positions[(-20 if r.is_reverse else 19)]
        gene, num, change = r.qname.split("|")[:3]
        outf.write("\t".join(map(str, (samfile.getrname(r.tid), snp_pos,
                                        change, gene, num))) + "\n")
    outf.close()


@files((os.path.join("bowtie_index", "S288C.fsa"),
       "../snps/SNP_positions.txt"),
       os.path.join("bowtie_index", "S288C_ambiguous.fsa"))
def make_fasta_ambiguous(infiles, outfile):
    fasta_infile, snp_infile = infiles

    d = SeqIO.to_dict(SeqIO.parse(fasta_infile, "fasta"))
    chr_dict = dict((i, s.seq) for i, s in d.items())
    comp_dict = dict((i, s.seq.complement()) for i, s in d.items())

    lst_dict = dict((i, list(s.seq)) for i, s in d.items())

    equals = []
    with open(snp_infile) as inf:
        for l in inf:
            chr, pos, change, ORF, ind = l.strip().split("\t")
            pos = int(pos)

            s = comp_dict[chr] if ORF.endswith("C") else chr_dict[chr]
            match = (change[0] == s[pos])

            # turn into an ambiguous base
            lst_dict[chr][pos] = "N"

            equals.append(match)

    # write out new fasta file
    new_seqs = [SeqRecord.SeqRecord(id=i, seq=Seq.Seq("".join(l)),
                                    description="")
                    for i, l in lst_dict.iteritems()]
    SeqIO.write(new_seqs, outfile, "fasta")

    print float(equals.count(True)) / len(equals)


@follows(SNP_positions, sort_reads, mkdir("snps/snp_summary"))
@transform("sorted_mapped/*.bam.bam", regex("sorted_mapped/(.*).bam.bam"),
            r"snps/snp_summary/\1.txt")
def find_overlapping(infile, outfile):
    import pysam
    command = "python ../overlapping_snps.py %s %s" % (infile, outfile)
    os.system("qsub -cwd %s" % command)
    return
    outf = open(outfile, "w")
    s = pysam.Samfile(infile)
    for l in open("../snps/SNP_positions.txt"):
        chr, pos, change, gene, num = l[:-1].split("\t")
        pos = int(pos)
        num_overlaps = sum(1 for _ in s.fetch(chr, pos, pos + 1))
        out = l[:-1] + "\t" + str(num_overlaps)
        outf.write(out + "\n")
    outf.close()


@merge("snps/snp_summary/*", "snp_summary.txt")
def combine_snps(infiles, outfile):
    d = collections.defaultdict(list)
    for f in infiles:
        with open(f) as inf:
            name = os.path.split(f)[1].split(".")[0]
            for l in inf:
                d[name].append(l.strip().split("\t")[-1])
    with open(outfile, "w") as outf:
        outf.write("\t".join(d.keys()) + "\n")
	print d.values()
        for i in range(len(d.values()[0])):
            outf.write("\t".join(v[i] for v in d.values()))
            outf.write("\n")



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            verbose=options.verbose)

    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multiprocess    = options.jobs,
                            logger          = stderr_logger,
                            verbose         = options.verbose)
