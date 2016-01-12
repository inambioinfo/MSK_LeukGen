.. |date| date::

*******************
MSK_LeukGen Scripts
*******************

:authors: Komal Rathi
:contact: rathik@mskcc.org
:organization: Memorial Sloan-Kettering
:status: This is "work in progress"
:date: |date|

.. meta::
   :keywords: rnaseq, scripts
   :description: RNAseq and Variant calling pipelines.

Welcome to **MSK_LeukGen** at luna. Here you can find scripts for RNASeq analysis and variant calling as well as other helper scripts to help you manipulate file types. These scripts were written for Leukemia Genomics Center at MSK.

* rnaseq contains scripts to *read fastq/bam files*, *perform QC*, *align* and *count* reads, and provide a *summary* file for *QC and alignment statistics*.
* rnaseqDownstreamAnalysis contains R scripts for differential gene expression as well as generating PDF reports containing heatmaps, venn diagrams, box plots, scatter plots and more.
* pileup contains scripts to generate pileups from bam files, call variants and generate PDF reports for the comparison of tumor and matched normal samples from a single individual.
* bsub_tests has a single job & array-job script to test your bsub command.
* gtf2bed has scripts that convert GTF file to an easy to use tab-delimited format that can be used to annotate your results.
* helper_scripts has scripts that help you change colors and themes of plots generated in R, get current logging time to add to your script outputs, get version of perl modules you are using (and more will be added).