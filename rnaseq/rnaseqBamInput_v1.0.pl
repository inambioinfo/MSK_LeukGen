#!/usr/bin/perl
# rnaseqBamInput_v1.0.pl
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2015-11-17_16:45:45
# Function: This script takes a bam input and performs read-counting using htseq-count
# Usage: 
# Usage1 is the stand alone usage for this program. Usage2 is when you call this script through rnaseq_wrapper.pl
# Usage1: perl rnaseqBamInput_v1.0.pl rnaseq.params <genome: either Mouse or Human> <sample_name>
# Usage2: perl rnaseqBamInput_v1.0.pl <key1 value1 key2 value2...keyN valueN> <genome: either Mouse or Human> <sample_name>

# load dependencies
# depends on tools version 1.0
require('/home/rathik/scripts/tools_v1.0.pl');
require('/home/rathik/scripts/getLoggingTime.pl');
use Path::Iterator::Rule;
use File::Find::Rule; 
use Text::Glob;
use Bio::DB::Sam;
use File::Basename;

# start time of the script
my $start_run = time();

# script version $0
($version = $0) =~ s/\.[^.]+$//;
print $version,"\n";

# timestamp getLoggingTime()

# get sample name, genome and parameters
my $sample = pop @ARGV;
my $genome = pop @ARGV;
my %param;
%param = @ARGV;

# check genome and set parameters accordingly
if($genome=="Mouse"){
	$ensembl = $mm10ensembl;
	$gencode = $mm10gencode;
	$stardb = $mm10Star;
} elsif($genome=="Human"){
	$ensembl = $hg19ensembl;
	$gencode = $hg19gencode;
	$stardb = $hg19Star;
}

# run htseq-count with Ensembl annotation
print "\nCounting reads using Htseq-count for Ensembl GTF...\n";
my $htseqe = "$python $htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample.".bam $ensembl > $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_htseq_ensembl.counts";
print $htseqe,"\n";
system($htseqe);

# run htseq-count with Gencode annotation
print "\nCounting reads using Htseq-count for Gencode GTF...\n";
my $htseqg = "$python $htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample.".bam $gencode > $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_htseq_gencode.counts";
print $htseqg,"\n";
system($htseqg);

# end time of the script
my $end_run = time();

# total time taken 
my $run_time = $end_run-$start_run;
print "Job took $run_time seconds";

# end of script
