#!/usr/bin/perl -w

# rnaseq_wrapper_v1.2.pl
# Version: 1.2
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2016-05-01_14:45:11
# Function: This is a wrapper script for running rnaseq pipeline
# Usage: perl rnaseq_wrapper.pl rnaseq.params

# changes in v1.2
# * calls rnaseqFastqInput_v1.2.pl
# * runs project based or sample based analysis
# * sorts merged fastq files to match read IDs

# dependencies
require('/home/rathik/scripts/getLoggingTime.pl');
use Data::Dumper;
use Path::Iterator::Rule;
use File::Find::Rule;
use Getopt::Long;
use Pod::Usage;


GetOptions(
           q(help) => \$help,
           'f|file=s' => \$paramfile,
           'a|analysis=s' => \$analysis
           ) or pod2usage(q(-verbose) => 1);
pod2usage(q(-verbose) => 1) if $help;


# read parameters file
# my $paramfile = $ARGV[0];
# my $analysis = $ARGV[1];

# open the file and split by =>
# save into a hash as key & value pairs
open PARAM, $paramfile or die print $!;
my %param;
while(<PARAM>)
{
    chomp;
    my @r = split('=>');
    $param{$r[0]}=$r[1];
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

# ask user if the input is bam or fastq
# call the corresponding script
print "Do you have fastq or bam files? Enter 1 for fastq and 2 for bam:";
my $userword = <STDIN>;
chomp $userword;
exit 0 if ($userword eq "");

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

    if($dir){
        # make output & log folder in sample directory
        print "Making directories in: ",$dir,"\n";
        $param{'OUTDIR'} = $dir.'/output';
        $param{'LOGDIR'} = $dir.'/log';
    } elsif($param{'OUTROOT'}){
        # make output & log folder in specified output directory
        print "Making directories in: ",$param{'OUTROOT'},"\n";
        $param{'OUTDIR'} = $param{'OUTROOT'}.'/output';
        $param{'LOGDIR'} = $param{'OUTROOT'}.'/log';
    } else {
        # make output & log folder in project directory
        print "Making directories in: ",$param{'PROJECTNAME'},"\n";
        $param{'OUTDIR'} = $param{'PROJECTNAME'}.'/output';
        $param{'LOGDIR'} = $param{'PROJECTNAME'}.'/log';
    }

    # if the directories does not exist, create directories
    system("mkdir $param{'OUTDIR'}") unless (-d $param{'OUTDIR'});
    system("mkdir $param{'LOGDIR'}") unless (-d $param{'LOGDIR'});

    if($userword==1){
        print "Input is fastq file...\n";
        system('bsub','-J',$samp,'-oo',"$param{'LOGDIR'}/".$samp."_".getLoggingTime().".out",'-eo',"$param{'LOGDIR'}/".$samp."_".getLoggingTime().".err",
               'perl','/home/rathik/scripts/rnaseqFastqInput_v1.2.pl',%param,$genome,$samp);
    } elsif($userword==2){
        print "Input is bam file...\n";
        system('bsub','-J',$samp,'-oo',"$param{'LOGDIR'}/".$samp."_".getLoggingTime().".out",'-eo',"$param{'LOGDIR'}/".$samp."_".getLoggingTime().".err",
               'perl','/home/rathik/scripts/rnaseqBamInput_v1.2.pl',%param,$genome,$samp);
    } else {
        print "Please enter a valid value: 1 for Fastq input and 2 for BAM input...\n";
        print "Exiting program...\n";
        exit;
    }
}

=head1 SYNOPSIS

    rnaseq_wrapper_v1.2.pl [options]

        required arguments:

            -file     | -f     parameters file
            -analysis | -a    [sample|project]

=head1 DESCRIPTION

    This is a wrapper for a perl script to do:
    1. QC using FASTQC
    2. Alignment using STAR
    3. Counting using Htseq-count

=cut
