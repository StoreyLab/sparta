#! /bin/sh

# Sample SPARTA analysis script
# This script uses SPARTA to process paired-end RNA-seq data from Skelly et al.
# In this data, each read is of ambiguous parental yeast strain origin (either BY or RM allele).
# We first map the reads to both genomes with bowtie2, and then use SPARTA to separate the reads that mapped to both genomes.

# REQUIREMENTS:
# a unix-based system
# bowtie2 2.2.3 
# SRA toolkit (specifically, version 2.3.5 of the fastq-dump tool)
# (older versions of SRA toolkit and bowtie2 MAY work)

# download the publicly available skelly data
echo 'Downloading data...'
wget -P skelly http://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR309/SRR309119/SRR309119.sra

# convert the reads from SRA to FASTQ format
echo 'Converting from SRA to FASTQ...'
fastq-dump --split-files -O skelly skelly/SRR309119.sra

# build an index of the BY and RM genomes
echo 'Building indices...'
mkdir -p skelly/BY_index
mkdir -p skelly/RM_index
bowtie2-build-s genomes/BY.fsa skelly/BY_index/index
bowtie2-build-s genomes/RM.fa skelly/RM_index/index

# align the RNA-seq reads to both genomes using bowtie2
echo 'Aligning Reads...'
bowtie2-align-s -x skelly/BY_index/index -1 skelly/SRR309119_1.fastq -2 skelly/SRR309119_2.fastq -S skelly/by.sam
bowtie2-align-s -x skelly/RM_index/index -1 skelly/SRR309119_1.fastq -2 skelly/SRR309119_2.fastq -S skelly/rm.sam

# use SPARTA to separate the reads that map to both genomes
echo 'Running SPARTA'
python sparta.py skelly/by.sam skelly/rm.sam -pe -n BY RM -o skelly/output