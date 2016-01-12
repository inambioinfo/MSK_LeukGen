*******************
Usage Documentation
*******************

:authors: Komal Rathi
:contact: rathik@mskcc.org
:organization: Memorial Sloan-Kettering
:status: This is "work in progress"
:date: |date|

.. meta::
   :keywords: rnaseq, usage
   :description: MSK_LeukGen's RNASeq analysis usage docs.

Introduction
============

RNASeq Pipeline
---------------

tools.pl
^^^^^^^^

Paths of all software/tools/reference files that we use frequently in our scripts. For e.g. python, R, perl, star, fastqc, reference fasta, gtf annotation etc.

.. note::
    This is for those who forget to update their bash_profile.

rnaseq_wrapper_v1.0.pl
^^^^^^^^^^^^^^^^^^^^^^

It the wrapper script for rnaseqFastqInput_v1.0.pl and rnaseqBamInput_v1.0.pl. It uses bsub for submitting jobs, one at a time, sample by sample.

rnaseqFastqInput_v1.0.pl
^^^^^^^^^^^^^^^^^^^^^^^^

Takes fastq files and performs fastqc, star alignment and htseq-count (using Ensembl and Gencode annotations).

rnaseqBamInput_v1.0.pl
^^^^^^^^^^^^^^^^^^^^^^

Takes bam files and performs read counting using htseq-count (using Ensembl and Gencode annotations).

rnaseq.params
^^^^^^^^^^^^^

This is the parameters file. It contains paths of project folder, sample information file, threads to use and such.

Step 1: Make sample information file

* Format: comma separated sample information
* Extension: .txt
* Required Columns: Source_name, Species
* Optional Columns: Source_type, Cell, Tissue etc
* Source_name: Sample Name
* Species: Either Human or Mouse
* Source_type: Either Tumor or Normal
* Example: `samples.txt`_.

.. include:: ../rnaseq/samples.txt

Step 2: Make parameters file

* Format: => separated key-value pairs
* Extension: .params
* Example: `rnaseq.params`_.

Step 3: Run wrapper script::

    perl rnaseq_wrapper_v1.0.pl rnaseq.params

The script will ask you if your input files are bam or fastq. Enter 1 for fastq and 2 for bam.

FastQC Summary
--------------

fastqc_parser.py
^^^^^^^^^^^^^^^^

Takes input as the name of fastqc output folder (that stores all fastqc output sub-folders) and creates a csv file in the same folder, summarizing results of all samples in one file.

Usage::

    python fastqc_parse.py <fastqc output folder>

Annotate using GTF
------------------

convertGTF2txt.pl
^^^^^^^^^^^^^^^^^

This is a very efficient & fast script that will convert a complex annotation like GTF to a simple one, extracting columns that are used frequently in genomic analysis. It will output useful information about a gene like genename, geneid, genetype, strand, chromosome, start and end. It uses an R script getgenecoord.R to summarize the start and end to gene level. (This is important because the start and end in GTF files are for exon, transcript and gene levels and could be quite confusing).

Usage::

    perl convertGTF2txt.pl <path/to/file.gtf> <path/to/outdir/>

.. references
.. _rnaseq.params: https://raw.githubusercontent.com/komalsrathi/MSK_LeukGen/master/rnaseq/rnaseq.params
.. _samples.txt: https://raw.githubusercontent.com/komalsrathi/MSK_LeukGen/master/rnaseq/samples.txt

.. attention::
gene counts output using STAR
STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which
correspond to different strandedness options:
column 1: gene ID
column 2: counts for unstranded RNA-seq
column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
