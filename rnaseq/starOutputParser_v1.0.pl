#!/usr/bin/perl
# starOutputParser_v1.0.pl
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-12-01_12:10:52
# Last Modified: 2015-12-01_12:10:52
# Function: This script goes through a list of star's final output files and summarizes the required statistics in one single file
# Usage: perl starOutputParser.pl <project directory>

# load dependencies
use Path::Iterator::Rule;
use File::Find::Rule;
use File::Basename;
require('/home/rathik/scripts/getLoggingTime.pl');

# project directory
$dir = $ARGV[0];
@dir_paths = File::Find::Rule->file->name("*final.out")->in("$dir");

# open output file to write star output summary
$output = $dir."/star_alignment_summary_".getLoggingTime().".txt";
print $output,"\n";
open($oh, ">", $output) or die ("Could not open '$output'");
print $oh "Filename\tTotalReads\tReadLength\tUniqueMapped\tPercentUniqueMapped\n";

# traverse each file, get statistics and write out to summary file
foreach (@dir_paths){

	# get file basename
	$filename = basename($_);
	push @arr, $filename;

	# open the sample output file
	open($fh, "<", $_) or die ("Unable to open file");

	# start parsing
	while ($line = <$fh>){
		chomp $line;
		$line =~ s/^\s+|\s+$//g;
		
		# split each pattern matched line by space, followed by | and a tab character
		if ($line =~ "Number of input reads"){
			@tmparray = split / \|\t/,$line;
			push @arr, $tmparray[1];
		} elsif($line =~ "Average input read length"){
			@tmparray = split / \|\t/,$line;
			push @arr, $tmparray[1]; 
		} elsif($line =~ "Uniquely mapped reads number"){
			@tmparray = split / \|\t/,$line;
                        push @arr, $tmparray[1];
		} elsif($line =~ "Uniquely mapped reads %"){
			@tmparray = split / \|\t/,$line;
                        push @arr, $tmparray[1];
		} 
	}

	# join the result in tab delimited manner
	$result = join("\t",@arr);

	# print the statistics to file
	print $oh $result;		
	print $oh "\n";

	@arr = (); # set array back to null
	close $fh; # close file handle on current star output file
}

# close file handle on output file
close $oh;


