#!/usr/bin/perl
# rnaseqFastqInput_v1.0.pl
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2015-11-17_16:45:45
# Function: This script takes a fastq input and performs fastqc, star alignment and read-counting using htseq-count
# Usage: 
# Usage1 is the stand alone usage for this program. Usage2 is when you call this script through rnaseq_wrapper.pl
# Usage1: perl rnaseqFastqInput_v1.0.pl rnaseq.params <genome: either Mouse or Human> <sample_name>
# Usage2: perl rnaseqFastqInput_v1.0.pl <key1 value1 key2 value2...keyN valueN> <genome: either Mouse or Human> <sample_name>

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

# get R1 and R2 fastq.gz files
@R1_paths;
@R2_paths;

# check genome and set parameters accordingly
if($genome=="Mouse"){
	$ensembl = $mm10ensembl;
	$gencode = $mm10gencode;
	$stardb = $mm10Star;
} 
if($genome=="Human"){
	$ensembl = $hg19ensembl;
	$gencode = $hg19gencode;
	$stardb = $hg19Star;
}

# get sample directory paths	
my @dir_paths = File::Find::Rule->directory->name("$sample")->in("$param{'PROJECTNAME'}");
$dir = $dir_paths[0];

# retrieve fastq files from the folder
@R1_paths = File::Find::Rule->file->name("*_R1*.fastq.gz")->in($dir);
@R2_paths = File::Find::Rule->file->name("*_R2*.fastq.gz")->in($dir);

# if there are no R1 or R2 files then exit
if(scalar(@R1_paths)==0 || scalar(@R2_paths)==0){
	print "Insufficient data!!";
	exit;
}

$R1_paths = join(' ',@R1_paths);
$R2_paths = join(' ',@R2_paths);

# merge fastq files
# merge R1 files
print "\nMerging Fastq R1 files:\n";
my $mergeR1 = "cat $R1_paths > $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_R1.fastq.gz";
print $mergeR1,"\n";
system($mergeR1);

# merge R2 files
print "\nMerging Fastq R2 files:\n";
my $mergeR2 = "cat $R2_paths > $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_R2.fastq.gz";
print $mergeR2,"\n";
system($mergeR2);

# fastqc on fastq files 
print "\nQC with Fastqc:\n";
my $fastqc = "$fastqc -q --extract --threads=2 --outdir=$param{'OUTDIR'} $param{'OUTDIR'}/".$sample."*.fastq.gz";
print $fastqc,"\n";
system($fastqc);

# align with star - output is bam format
print "\nAligning with STAR:\n";
my $staralign = "$star --outSAMtype BAM Unsorted --runThreadN $param{'THREADS'} --genomeDir $stardb --outFileNamePrefix $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_ --readFilesIn $param{'OUTDIR'}/".$sample."_".$version."*_R1.fastq.gz $param{'OUTDIR'}/".$sample."_".$version."*_R2.fastq.gz --readFilesCommand zcat --quantMode GeneCounts --sjdbGTFfile $ensembl";
print $staralign,"\n";
system($staralign);

# remove fastq.gz files as soon as star is done
print "\nRemove fastq files:\n";
my $removefastq = "rm $param{'OUTDIR'}/".$sample."_".$version."*.fastq.gz";
print $removefastq,"\n";
system($removefastq);

# run htseq-count with Ensembl annotation
print "\nCounting reads using Htseq-count for Ensembl GTF...\n";
my $htseqe = "$python $htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample."_".$version."*.bam $ensembl > $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_htseq_ensembl.counts";
print $htseqe,"\n";
system($htseqe);

# run htseq-count with Gencode annotation
print "\nCounting reads using Htseq-count for Gencode GTF...\n";
my $htseqg = "$python $htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample."_".$version."*.bam $gencode > $param{'OUTDIR'}/".$sample."_".$version."_".getLoggingTime()."_htseq_gencode.counts";
print $htseqg,"\n";
system($htseqg);

# end time of the script
my $end_run = time();

# total time taken 
my $run_time = $end_run-$start_run;
print "Job took $run_time seconds";

# end of script
