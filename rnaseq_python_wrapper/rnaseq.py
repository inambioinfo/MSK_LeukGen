
"""RNAseq pipeline."""

# python
import os
import os.path
import sys
import subprocess
import csv

# third party
import click

# leuktools
from leuktools.py import utils
from leuktools.py import options
from leuktools.py import constants

NAME = "RNASEQ"
VERSION = "0.1.0"

@click.command()
@options.analysis
@options.new_analysis
@click.option(
    "--project", "-p",
    required=True,
    help="Leukgen project ID.",
    metavar="<string>")
@click.option(
    "--inputtype", "-i",
    type=click.Choice(["bam", "fastq"]),
    help="Input file type. [bam|fastq]",
    required=True,
    metavar="<string>")
@click.option(
    "--analysistype", "-a",
    type=click.Choice(["sample", "project"]),
    help="Analysis type. [sample|project]",
    required=True,
    metavar="<string>")
@click.option(
    "--gtf", "-g",
    type=click.Choice(["ensembl", "gencode"]),
    help="Annotation. [ensembl|gencode]",
    required=True,
    metavar="<string>")
@click.option(
    "--outdir", "-o",
    type=click.Path(exists=False, file_okay=False, dir_okay=True, writable=True, readable=True),
    help="path to output directory",
    metavar="<path>")
@click.option(
    "--cores", "-c",
    type=click.INT,
    help="Cores. DEFAULT=1",
    default=1,
    required=False,
    metavar="<integer>")
@click.option(
    "--keep",
    is_flag=True,
    help="Add this flag to keep intermediate files.")
@options.commit
def rnaseq(project, inputtype, analysistype, cores, gtf, keep, outdir, commit, analysis, new_analysis):
    """Run Alignment and Readcount for RNASeq data."""
    data = utils.get_instance(identifier=project, model="Project")
    if data is None:
        print('Permission denied! You do not have access to the project!')
        sys.exit()
    else:
        # get primary key of samples for that project
        samples = data.get("workflow_set")
        samples = [d['pk'] for d in samples]

    # if outdir is not specified
    # create analysis instance and store results in it
    if outdir is None:
        if commit:
            analysis = utils.get_analysis(NAME, VERSION, [project], samples, new=new_analysis, identifier=analysis)
            analysesdir = utils.get_resultsdir(analysis["pk"])
            os.putenv("OUTDIR", analysesdir)
            os.system("mkdir -p $OUTDIR")
            click.secho("\nAnalysis Primary Key:\n", fg="blue")
            click.echo("\t%s\n" % analysis["pk"])
        else:
            analysesdir = "/path/to/results"
    else:
        # if outdir is specified, write output to it
        analysesdir = outdir
        if commit:
            os.putenv("OUTDIR", analysesdir)
            os.system("mkdir -p $OUTDIR")

    click.secho("\nOutput directory:\n", fg="blue")
    click.echo("\t%s\n" % analysesdir)

    # get name and path to sample files
    rows = []
    for sample in samples:
        sample = str(sample)
        tdata = utils.get_instance(identifier=sample, model="Workflow")
        twf = "W%07dF" % tdata["pk"]     # sample name
        tsp = tdata["species"]           # species
        tsource = tdata["source_type"]   # source
        textid = tdata["ext_id"]
        rows.append({'source_name': twf,
                    'species': tsp,
                    'source_type': tsource,
                    'ext_id': textid})

    # define samples file
    infile = analysesdir + '/samples.txt'
    if commit:
        with open(infile, 'w') as csvfile:
            fieldnames = ['source_name', 'species', 'source_type', 'ext_id']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    script = os.path.join(constants.LEUKTOOLS_ROOT, 'leuktools/nested/run/perl', 'rnaseq_wrapper_v1.2.pl')

    if keep:
        delete = "False"
    else:
        delete = "True"

    cmdstring = "perl {0} -p {1} -s {2} -c {3} -o {4} -a {5} -t {6} -g {7} -d {8}".format(script, constants.WORKFLOWS_ROOT, infile, cores, analysesdir, analysistype, inputtype, gtf, delete)
    click.secho("\nCommand:\n", fg="blue")
    click.echo("\t%s\n" % cmdstring)
    if commit:
        os.system(cmdstring)

    if not commit:
        click.secho("Add --commit to submit.\n", fg="green", blink=True)
        raise click.Abort

if __name__ == "__main__":
    rnaseq()
