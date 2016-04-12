#!/usr/bin/perl -w

# rnaseq_wrapper_v1.0.pl
# Version: 1.2
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2016-02-17_14:45:11
# Function: This is a wrapper script for running rnaseq pipeline
# Usage: perl rnaseq_wrapper.pl rnaseq.params

# changes in v1.0
# * calls rnaseqFastqInput_v1.0.pl
# * runs project based or sample based analysis
# * sorts merged fastq files to match read IDs

# dependencies
use Data::Dumper;
use Path::Iterator::Rule;
use File::Find::Rule;
use Getopt::Long;
use Pod::Usage;

# default number of threads
$threads = 1;

# get commandline arguments
GetOptions(
           q(help)          => \$help,
           'f|file=s'       => \$paramfile,
           'a|analysis=s'   => \$analysis,
           't|type=s'       => \$type,
           'p|project=s'    => \$projectname,
           's|samplelist=s' => \$samplelist,
           'c|cores=i'      => \$threads,
           'o|output=s'     => \$output,
           'g|gtf=s'        => \$gtf,
           'd|delete=s'     => \$delete
           );
           # ) or pod2usage(q(-verbose) => 1);
pod2usage(q(-verbose) => 1) if $help;

# define a hash
my %param;

# check parameters
# check if type is provided
if(!$type || ($type ne "fastq" && $type ne "bam")){
    print "Please enter the type: fastq or bam","\n";
    exit;
} else {
    chomp $type;
}

# check if parameters file or alternative arguments are supplied
if(!$paramfile){
    if(!$projectname || !$samplelist){
        print "Please enter the configuration file or projectname, samplelist and threads to use","\n";
        exit;
    } else {
        %param = (
                  PROJECTNAME => $projectname,
                  SAMPLELIST => $samplelist,
                  THREADS => $threads
                  )
    }
} else {
    open PARAM, $paramfile or die print $!;
    while(<PARAM>)
    {
        chomp;
        my @r = split('=>');
        $param{$r[0]}=$r[1];
    }
}

# check if analysis level is specified
if(!$analysis){
    print "Please enter the anlaysis type: sample or project","\n";
    exit;
}

# add output directory to params if specified
if($output){
    $param{'OUTROOT'} = $output;
}

# open file with filename which is value of the key SAMPLELIST
open FILE, $param{'SAMPLELIST'} or die

# create an array samples
my @samples;

# splitting each line based on ',' and storing in an array @r
# pushing the reference of this array in another array @samples
while(<FILE>){
    chomp;
    my @r = split(',');
    push(@samples,\@r);
}

# remove header
shift @samples;

# process each sample
foreach (@samples)
{
    # get and print sample name
    my $samp = $_->[0];
    my $genome = $_->[1];

    print "\nProcessing Sample:",$samp,"\n";

    # if analysis is sample based
    if($analysis eq "sample"){
        # get the sample directory from project (root) directory
        @dir_paths = File::Find::Rule->directory->name("$samp")->in("$param{'PROJECTNAME'}");

        # check if the sample directory exists or not
        if(scalar(@dir_paths)==0){
            print "$samp doesn't exist!!";
            next;
        }

        # sample directory
        $dir = $dir_paths[0];

        # print sample directory path
        print "Sample directory:",$dir,"\n";
    }

    if($param{'OUTROOT'}){
        # make output & log folder in specified output directory
        print "Making directories in: ",$param{'OUTROOT'},"\n";
        $param{'OUTDIR'} = $param{'OUTROOT'}.'/output';
        $param{'LOGDIR'} = $param{'OUTROOT'}.'/log';
    } elsif($dir){
        # make output & log folder in sample directory
        print "Making directories in: ",$dir,"\n";
        $param{'OUTDIR'} = $dir.'/output';
        $param{'LOGDIR'} = $dir.'/log';
    } else {
        # make output & log folder in project directory
        print "Making directories in: ",$param{'PROJECTNAME'},"\n";
        $param{'OUTDIR'} = $param{'PROJECTNAME'}.'/output';
        $param{'LOGDIR'} = $param{'PROJECTNAME'}.'/log';
    }

    # if the directories does not exist, create directories
    system("mkdir -p $param{'OUTDIR'}") unless (-d $param{'OUTDIR'});
    system("mkdir -p $param{'LOGDIR'}") unless (-d $param{'LOGDIR'});

    if($type eq "fastq"){
        print "Input is fastq file...\n";
        system('bsub','-J',$samp,'-oo',"$param{'LOGDIR'}/".$samp.".out",'-eo',"$param{'LOGDIR'}/".$samp.".err",
               'perl','/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/perl/rnaseqFastqInput_v1.0.pl',%param, $genome, $samp, $gtf, $delete);
    } elsif($type eq "bam"){
        print "Input is bam file...\n";
        system('bsub','-J',$samp,'-oo',"$param{'LOGDIR'}/".$samp.".out",'-eo',"$param{'LOGDIR'}/".$samp.".err",
               'perl','/ifs/work/leukgen/local/opt/leuktools/development/leuktools/nested/run/perl/rnaseqBamInput_v1.0.pl',%param, $genome, $samp, $gtf, $delete);
    } else {
        print "Please enter a valid value: 1 for Fastq input and 2 for BAM input...\n";
        print "Exiting program...\n";
        exit;
    }
}

=head1 SYNOPSIS

    rnaseq_wrapper_v1.0.pl [options]

        required arguments:

            -file     | -f                              parameters file
            -analysis | -a                              [sample|project]
            -type     | -t                                   [fastq|bam]
            -gtf      | -g                           gtf annotation file
            -delete   | -d     delete intermediate analysis [True|False]

        required arguments when parameters file (-f or -file) is not provided:

           -project   | -p      project name
           -samplelist| -s        samplelist
           -cores     | -c   number of cores

        optional arguments:

            -output   | -o  output directory

=head1 DESCRIPTION

    This is a wrapper for a perl script to do:
    1. QC using FASTQC
    2. Alignment using STAR
    3. Counting using Htseq-count

=cut
