#!/usr/bin/perl -w

# rnaseq_wrapper_v1.1.pl
# Version: 1.1
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2016-05-01_14:45:11
# Function: This is a wrapper script for running rnaseq pipeline
# Usage: perl rnaseq_wrapper.pl rnaseq.params

# changes in v1.1
# * calls rnaseqFastqInput_v1.1.pl instead of rnaseqFastqInput_v1.0.pl
# * this change is to run the pipeline on a single folder of fastq files
# * sorts merged fastq files to match read IDs

# dependencies
require('/home/rathik/scripts/getLoggingTime.pl');
use Data::Dumper;
use Path::Iterator::Rule;
use File::Find::Rule;

# read parameters file
my $paramfile = $ARGV[0];

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

    # we don't need this anymore as we don't have individual sample folders
	# make output & log folder in sample folder
    # @dir_paths = File::Find::Rule->directory->name("$samp")->in("$param{'PROJECTNAME'}");

	# check if the sample directory exists or not
    # if(scalar(@dir_paths)==0){
    # 	print "$samp doesn't exit!!";
	# 	next;
    # }

	# sample directory
    # $dir = $dir_paths[0];	

	# make output folder in sample directory
	$param{'OUTDIR'} = $dir.'/output';
	$param{'LOGDIR'} = $param{'OUTDIR'}.'/log';
	
	# if the directories does not exist, create directories
	system("mkdir $param{'OUTDIR'}") unless (-d $param{'OUTDIR'});
	system("mkdir $param{'LOGDIR'}") unless (-d $param{'LOGDIR'});
	
	if($userword==1){
		print "Input is fastq file...\n";
		print "Using rnaseqFastqInput_v1.1.pl...\n";
		system('bsub','-J',$samp,'-oo',"$param{'LOGDIR'}/".$samp."_rnaseqFastqInput_v1.1_".getLoggingTime().".out",'-eo',"$param{'LOGDIR'}/".$samp."_rnaseqFastqInput_v1.1_".getLoggingTime().".err",'perl','rnaseqFastqInput_v1.1.pl',%param,$genome,$samp);
	} elsif($userword==2){
		print "Input is bam file...\n";
		print "Using rnaseqBamInput_v1.0.pl...\n";
		# system('bsub','-J',$samp,'-o',"$param{'LOGDIR'}/$samp.out",'-e',"$param{'LOGDIR'}/$samp.err",'perl','rnaseq_bamInput.pl',%param,$genome,$samp);
		system('bsub','-J',$samp,'-oo',"$param{'LOGDIR'}/".$samp."_rnaseqBamInput_v1.0_".getLoggingTime().".out",'-eo',"$param{'LOGDIR'}/".$samp."_rnaseqBamInput_v1.0_".getLoggingTime().".err",'perl','rnaseqBamInput_v1.0.pl',%param,$genome,$samp);
	} else {
		print "Please enter a valid value: 1 for Fastq input and 2 for BAM input...\n";
		print "Exiting program...\n";
        	exit;
	}
}
