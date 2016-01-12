#!/usr/bin/perl -w

# convertGTF2txt_v1.0.pl
# Version: 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: ##------ Fri Nov  6 14:27:25 2015 ------ #
# Last Modified: ##------ Wed Nov 11 16:48:00 2015 ------##
# Usage: perl convertGTF2txt_v1.0.pl <path to GTF file> <path to output directory>
# This script also calls an R script called getgenecoord.R

require "tools_v1.0.pl";

# path to GTF file
my $gtf = $ARGV[0];

# path to output directory
my $output = $ARGV[1];

# output file name
my $outputfile = $gtf;
$outputfile =~ s{.*/}{};
$outputfile =~ s{\.[^.]+$}{};

# run command
print "\nConverting GTF to BED and TXT to make your life easier...:\n";
my $cmd = "grep -P '\\tgene\\t' $gtf | awk '{ 
	for (i = 1; i <= NF; i++){	
     		if (\$i==\$1 || \$i==\$4 || \$i==\$5 || \$i==\$7 || \$(i-1) ~ /gene_id|gene_name|gene_type|gene_biotype/) {
     		printf \"\%s \", \$i 
		}	
	}
	print \"\"
}' - | sed -e 's/;//g' -e 's/\"//g' -e 's/ /\\t/g' -e 's/\\t\$//g' > $output/$outputfile.txt";
print $cmd,"\n";
system($cmd);

print "\n Calling R script getgenecoord.R to do some cleaning up....\n";
$cmd = "$rscript --vanilla getgenecoord.R $output/$outputfile.txt";
print $cmd,"\n";
system($cmd);
