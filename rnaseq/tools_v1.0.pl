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
$hg19Star="/ifs/depot/assemblies/H.sapiens/hg19/index/star/2.4.1d/gencode/v18/overhang74";
$hg19fasta="/ifs/depot/assemblies/H.sapiens/hg19/hg19.fasta";
$hg19gencode="/ifs/work/leukgen/ref/homo_sapiens/37/gencode/18/gencode.v18.annotation.gtf";
$hg19ensembl="/ifs/work/leukgen/ref/homo_sapiens/37/ensembl/75/Homo_sapiens.GRCh37.75.gtf";
$hg19refseq="/ifs/work/leukgen/ref/homo_sapiens/37/refseq/hg19_RefSeq.bed";

# mouse
$mm10Star="/ifs/depot/assemblies/M.musculus/mm10/index/star/2.4.1d/ensembl/v80/overhang74";
$mm10fasta="/ifs/depot/assemblies/M.musculus/mm10/mm10.fasta";
$mm10gencode="/ifs/work/leukgen/ref/mus_musculus/38/gencode/M7/gencode.vM7.annotation.gtf";
$mm10ensembl="/ifs/work/leukgen/ref/mus_musculus/38/ensembl/80/Mus_musculus.GRCm38.80_canonical_chromosomes.gtf";
$mm10refseq="/ifs/work/leukgen/ref/mus_musculus/38/refseq/mm10_RefSeq.bed";

###################
# reference tools #
###################
$star="/ifs/e63data/levinelab/rapaport/software/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR";
$htseqcount="/opt/common/CentOS_6/htseq/HTSeq-0.6.1p1/bin/htseq-count";
$samtools="/opt/common/CentOS_6/samtools/samtools-1.2/samtools";
$fastqc="/opt/common/CentOS_6/fastqc/FastQC_v0.10.1/fastqc";
$pairfq="/ifs/work/leukgen/opt/Pairfq/0.14.7/bin/pairfq";
$estimatelibcomp="/opt/common/CentOS_6/picard/picard-tools-1.96/EstimateLibraryComplexity.jar";
$rseqc="/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python2.7 /opt/common/CentOS_6-dev/RSeQC/RSeQC-2.6.2/scripts/infer_experiment.py";

######################
# reference binaries #
######################
# $python="/opt/common/CentOS_6/python/python-2.7.8/bin/python2.7";
# $perl="/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl";
# $r="/opt/common/CentOS_6/R/R-3.2.0/bin/R";
# $rscript="/opt/common/CentOS_6/R/R-3.2.0/bin/Rscript";
