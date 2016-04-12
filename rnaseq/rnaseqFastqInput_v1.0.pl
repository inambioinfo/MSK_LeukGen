#!/usr/bin/perl

# rnaseqFastqInput_v1.0.pl
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2016-02-16_14:29:49
# Function: This script takes a fastq input and performs fastqc, star alignment and read-counting using htseq-count
# Usage:
# Usage1 is the stand alone usage for this program. Usage2 is when you call this script through rnaseq_wrapper.pl
# Usage1: perl rnaseqFastqInput_v1.0.pl rnaseq.params <genome: either Mouse or Human> <sample_name>
# Usage2: perl rnaseqFastqInput_v1.0.pl <key1 value1 key2 value2...keyN valueN> <genome: either Mouse or Human> <sample_name>

# changes in v1.0
# * matches *.R1 and *.R2 and/or *_R1 and *_R2
# * does sample based or project based analysis
# * sorts merged fastq files to match read IDs
# * gets star db based on input of species and gtf

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
my $del = pop @ARGV;
my $ann = pop @ARGV;
my $sample = pop @ARGV;
my $genome = pop @ARGV;
my %param;
%param = @ARGV;

# get R1 and R2 fastq.gz files
@R1_paths;
@R2_paths;

# check genome and set parameters accordingly
if($genome=~/mouse/i){
    if($ann=~/ensembl/i){
        $gtf = $mm10ensembl;
        $stardb = $mm10Starensembl;
        $stargtf = $mm10ensembl;
    } else {
        $gtf = $mm10gencode;
        $stardb = $mm10Stargencode;
        $stargtf = $mm10gencode;
    }
}
if($genome=~/human/i){
    if($ann=~/ensembl/i){
        $gtf = $hg19ensembl;
        $stardb = $hg19Starensembl;
        $stargtf = $hg19ensembl;
    } else {
        $gtf = $hg19gencode;
        $stardb = $hg19Stargencode;
        $stargtf = $hg19gencode;
    }
}

print "Sample: ",$sample,"\n";
print "Genome: ",$genome,"\n";
print "GTF:    ",$gtf,"\n";
print "STAR db:",$stardb,"\n";

# htseq counts results stored in sample workflow folder
@dir_paths = File::Find::Rule->directory->name("$sample")->in("$param{'PROJECTNAME'}");
$param{'FINALOUTPUT'} = $dir_paths[0].'/data/htseq';
print "\nCreating htseq raw counts directory...\n";
my $createoutputdir = "mkdir -p $param{'FINALOUTPUT'}";
print $createoutputdir,"\n";
system($createoutputdir);

# retrieve fastq files from the folder
@R1_paths = File::Find::Rule->file->name("$sample*_1.fastq.gz")->in("$param{'PROJECTNAME'}");
@R2_paths = File::Find::Rule->file->name("$sample*_2.fastq.gz")->in("$param{'PROJECTNAME'}");

# if there are no R1 or R2 files then exit
if(scalar(@R1_paths)==0 || scalar(@R2_paths)==0){
    print "Insufficient data!!";
    exit;
}

$R1_paths = join(' ',@R1_paths);
$R2_paths = join(' ',@R2_paths);

# merge fastq files
# merge R1 files
print "\nMerging fastq R1 files...\n";
my $mergeR1 = "cat $R1_paths > $param{'OUTDIR'}/".$sample."_R1.fastq.gz";
print $mergeR1,"\n";
system($mergeR1);

# merge R2 files
print "\nMerging fastq R2 files...\n";
my $mergeR2 = "cat $R2_paths > $param{'OUTDIR'}/".$sample."_R2.fastq.gz";
print $mergeR2,"\n";
system($mergeR2);

# sort R1 and R2
print "\nSorting fastq files...\n";
my $sortfastq = "$pairfq makepairs -c gzip -f $param{'OUTDIR'}/".$sample."_R1.fastq.gz -r $param{'OUTDIR'}/".$sample."_R2.fastq.gz -fp $param{'OUTDIR'}/".$sample."_sorted_R1.fastq -rp $param{'OUTDIR'}/".$sample."_sorted_R2.fastq -fs $param{'OUTDIR'}/".$sample."_forwards.fastq -rs $param{'OUTDIR'}/".$sample."_reverses.fastq";
print $sortfastq,"\n";
system($sortfastq);

# fastqc on fastq files
print "\nQC with FASTQC...\n";
my $fastqc = "$fastqc -q --extract --threads=2 --outdir=$param{'OUTDIR'} $param{'OUTDIR'}/".$sample."_sorted_R*.fastq.gz";
print $fastqc,"\n";
system($fastqc);

# summarize fastqc
print "\nSummarizing fastqc output...\n";
my $fastqcparse = "python /ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/perl/fastqc_parse_v1.0.py $param{'OUTDIR'}";
print $fastqcparse,"\n";
system($fastqcparse);

# align with star - output is bam format
print "\nAligning with STAR...\n";
my $staralign = "$star --outSAMtype BAM Unsorted --runThreadN $param{'THREADS'} --genomeDir $stardb --outFileNamePrefix $param{'OUTDIR'}/".$sample."_ --readFilesIn $param{'OUTDIR'}/".$sample."_sorted_R1.fastq.gz $param{'OUTDIR'}/".$sample."_sorted_R2.fastq.gz --readFilesCommand zcat --quantMode GeneCounts --sjdbGTFfile $stargtf";
print $staralign,"\n";
system($staralign);

# summarize star
print "\nSummarizing star output...\n";
my $starparse = "perl /ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/perl/starOutputParser_v1.0.pl $param{'OUTDIR'}";
print $starparse,"\n";
system($starparse);

# run htseq-count
print "\nCounting reads with Htseq-count using $ann...\n";
# my $htseq = "$python $htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample."_*.bam $gtf > $param{'OUTDIR'}/".$sample."_".$ann."_htseq.counts";
my $htseq = "$python $htseqcount --format=bam -q -m union -s no $param{'OUTDIR'}/".$sample."_*.bam $gtf > $param{'FINALOUTPUT'}/".$sample."_".$ann."_htseq.counts";
print $htseq,"\n";
system($htseq);

# remove unwanted files and folders
if($del=~/true/i){
    print "\nCleaning up directory...\n";
    my $cleanup = "rm -rf $param{'OUTDIR'}/".$sample."*STAR* $param{'OUTDIR'}/".$sample."*.fastq.gz $param{'OUTDIR'}/".$sample."_*.bam $param{'OUTDIR'}/".$sample."*.out $param{'OUTDIR'}/".$sample."*.html $param{'OUTDIR'}/".$sample."*.zip $param{'OUTDIR'}/".$sample."*_SJ.out.tab";
    print $cleanup,"\n";
    system($cleanup);
}

# end time of the script
my $end_run = time();

# total time taken
my $run_time = $end_run-$start_run;
print "Job took $run_time seconds";

# end of script
