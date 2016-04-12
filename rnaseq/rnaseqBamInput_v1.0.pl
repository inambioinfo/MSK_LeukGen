#!/usr/bin/perl
# rnaseqBamInput_v1.0.pl
# Version: 1.2
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
require('/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/perl/tools_v1.0.pl');
use Path::Iterator::Rule;
use File::Find::Rule;
use Text::Glob;
use File::Basename;

# start time of the script
my $start_run = time();

# script version $0
my $filename = basename($0);
($version = $filename) =~ s/\.[^.]+$//;
print $version,"\n";

# get sample name, genome and parameters
my $ann = pop @ARGV;
my $sample = pop @ARGV;
my $genome = pop @ARGV;
my %param;
%param = @ARGV;

# check genome and set parameters accordingly
if($genome=~/mouse/i){
    if($ann=~/ensembl/i){
        $gtf = $mm10ensembl;
    } else {
        $gtf = $mm10gencode;
    }
    $stardb = $mm10Star;
    $stargtf = $mm10gencode;
}
if($genome=~/human/i){
    if($ann=~/ensembl/i){
        $gtf = $hg19ensembl;
    } else {
        $gtf = $hg19gencode;
    }
    $stardb = $hg19Star;
    $stargtf = $hg19gencode;
}

print "Sample: ",$sample,"\n";
print "Genome: ",$genome,"\n";
print "GTF:    ",$ann,"\n";

# run htseq-count
print "\nCounting reads with Htseq-count using $ann...\n";
my $htseq = "$htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample."_*.bam $gtf > $param{'OUTDIR'}/".$sample."_".$ann."_htseq.counts";
print $htseq,"\n";
system($htseq);

# remove bam files as soon as htseq is done
print "\nRemoving bam files...\n";
my $removebam = "rm $param{'OUTDIR'}/".$sample."_*.bam";
print $removebam,"\n";
system($removebam);

# remove unwanted files and folders
print "\nCleaning up directory...\n";
my $cleanup = "rm -rf *STAR* *out* ";
print $cleanup,"\n";
system($cleanup);

# end time of the script
my $end_run = time();

# total time taken
my $run_time = $end_run-$start_run;
print "Job took $run_time seconds";

# end of script
