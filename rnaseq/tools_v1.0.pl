# tools_v1.0.pl
# version 1.0
# Author: Komal S Rathi
# Institute: Memorial Sloan Kettering Cancer Center
# Created: 2015-11-13_11:52:36
# Last Modified: 2015-11-13_11:52:36
# Function: This script supplies the paths to various tools and files required by rnaseq pipeline

######################
# reference datasets #
######################

# human
$hg19Starensembl="/ifs/work/leukgen/ref/homo_sapiens/37/star/ensembl/75";
$hg19Stargencode="/ifs/work/leukgen/ref/homo_sapiens/37/star/gencode/18";
$hg19ensembl="/ifs/work/leukgen/ref/homo_sapiens/37/ensembl/75/Homo_sapiens.GRCh37.75.gtf";
$hg19gencode="/ifs/work/leukgen/ref/homo_sapiens/37/gencode/18/gencode.v18.annotation.gtf";

# mouse
$mm10Starensembl="/ifs/work/leukgen/ref/mus_musculus/38/star/ensembl/80";
$mm10Stargencode="/ifs/work/leukgen/ref/mus_musculus/38/star/gencode/M7";
$mm10ensembl="/ifs/work/leukgen/ref/mus_musculus/38/ensembl/80/Mus_musculus.GRCm38.80_canonical_chromosomes_nochr.gtf";
$mm10gencode="/ifs/work/leukgen/ref/mus_musculus/38/gencode/M7/gencode.vM7.annotation.gtf";

###################
# reference tools #
###################
$star="/ifs/work/leukgen/opt/star/2.4.2a/bin/STAR";
$htseqcount="/opt/common/CentOS_6/htseq/HTSeq-0.6.1p1/bin/htseq-count";
$samtools="/ifs/work/leukgen/opt/cgp/5.18.4/alleleCount/2.1.1/bin/samtools";
$fastqc="/ifs/work/leukgen/opt/fastqc/0.11.5/bin/fastqc";
$pairfq="/ifs/work/leukgen/opt/pairfq/0.14.7/bin/pairfq";

######################
# reference binaries #
######################
$python="/opt/common/CentOS_6/python/python-2.7.8/bin/python2.7";
# $perl="/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl";
# $r="/opt/common/CentOS_6/R/R-3.2.0/bin/R";
# $rscript="/opt/common/CentOS_6/R/R-3.2.0/bin/Rscript";
