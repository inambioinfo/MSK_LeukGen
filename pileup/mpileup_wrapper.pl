# perl wrapper script for mpileup -> nucleotide frequencies
# mpileup_wrapper.pl
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified: ##------ Wed Nov 11 16:48:00 2015 ------##
# Function: This script takes in a number of bam files and generates pileups for a given list of positions. 
# The pileups are then converted to easy to read format of nucleotide frequencies.
# Usage: perl mpileup_wrapper.pl <project directory of bam files> <positions file>

#!/usr/bin/perl
# use lib "/opt/common/CentOS_6/perl/perl-5.20.1/lib";
# use lib "/home/rathik/lib";
# use required modules
use Path::Iterator::Rule;
use File::Find::Rule; 
use Text::Glob;
use Bio::DB::Sam;
use File::Basename;

# get start time
my $start_run = time();

# get path of the project directory where bam files are located
$path = $ARGV[0];

# subroutine definition
my $cb = sub {
        my ($seqid, $pos, $pileups) = @_;
        my @pileuptmp;
        my @strandtmp;
        my @reftmp;
	
	# get reference base
        my $refbase = $sam->segment($seqid,$pos,$pos)->dna;

	# for the pileups, calculate nucleotide counts and depth on forward and reverse strands
        for my $pileup (@$pileups){
                my $al = $pileup->alignment; # get alignment at that position
                next if $pileup->indel or $pileup->is_refskip; # dont deal with indels
                my $strand = $al->strand; # get strand info at that particular alignment position
                push(@strandtmp,$strand); # push it in an array where strand info is stored
                my $qBase = substr($al->qseq, $pileup->qpos, 1); # get base at that position
                
		# set per base quality score threshold
		# my $qscore = $al->qscore->[$pileup->qpos];
		# next unless $qscore > 13; # per base quality score >= 13
		
		# set mapping quality score threshold
		# my @scores    = $a->qscore
		# my $match_qual= $a->qual;
		
		next if $pileup->indel or $pileup->is_refskip; # don't deal with indels
		next if $qbase =~ /[nN]/; # don't count Ns
                
		# count nucleotides on reverse strand
		if($strand==-1){
                        $qBase = lc $qBase;
                }

                # push the bases in an array of pileup
		push(@pileuptmp,$qBase);
        }

	# get total array of pileup
        $totalpileup = join("", @pileuptmp);
	my $depth_forward = () = $totalpileup =~ m/\p{Uppercase}/g; # depth on forward strand
	my $depth_reverse = () = $totalpileup =~ m/\p{Lowercase}/g; # depth on reverse strand
        my $depth = length($totalpileup); # total depth
        my $A = $totalpileup =~ tr/A//; # count A
        my $T = $totalpileup =~ tr/T//; # count T
        my $C = $totalpileup =~ tr/C//; # count C
        my $G = $totalpileup =~ tr/G//; # count G
        my $a = $totalpileup =~ tr/a//; # count a
        my $c = $totalpileup =~ tr/c//; # count c
        my $g = $totalpileup =~ tr/g//; # count g
	my $t = $totalpileup =~ tr/t//; # count t

	# print only for the snp of interest
        if($pos==$snp){
                print $oh "$seqid\t$pos\t$refbase\t$depth\t$depth_forward\t$depth_reverse\t$A\t$T\t$G\t$C\t$a\t$t\t$g\t$c\n";
        }
        @tmp=();
};

# get all bam files from the directory
my @full_paths = File::Find::Rule->file->name('*.bam')->in($path);

for my $file (@full_paths) {
	($base, $dir, $ext) = fileparse($file,'\..*');
	print "Processing: $base\n";

	# import the bam and fasta file
	$sam = Bio::DB::Sam->new(-bam => $file,
                            	-fasta => "/ifs/work/gabow/ref/gr37/gr37.fasta"
        			);

	# this will create the name of the output file
	# save the output in the same directory as the input folder
        # $output = $dir.$base.'.out';
        # use the below command for now
        $output = '/home/rathik/data/'.$base.'_pileup.out';
	
	# open the positions file
	$filename = $ARGV[1];
	open($fh, '<:encoding(UTF-8)', $filename)
  		or die "Could not open file '$filename' $!";

	# open the output file, print header line and then proceed
	open($oh, '>:encoding(UTF-8)', $output)
  		or die "Could not open file '$output' $!";
	print $oh "CHR\tPOS\tREF\tD\tD_forward\tD_reverse\tA_forward\tT_forward\tG_forward\tC_forward\tA_reverse\tT_reverse\tG_reverse\tC_reverse\n";
	
	# for the list of positions, call pileup subroutine
	$count = 0;
	while(my $row = <$fh>){
        	chomp $row;
        	@arr = split(/\t/, $row);
        	$id = $arr[0];
        	$snp = $arr[1];
        	$sam->pileup("$id:$snp-$snp", $cb);
        	$count++;
	}
}

# time taken by script
my $end_run = time();
my $total_time = $end_run-$start_run;
print "\nTotal time to process $count snps: $total_time seconds\n";

# close file handles
close $fh;
close $oh;
