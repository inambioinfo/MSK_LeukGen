
"""RNAseq downstream analysis pipeline."""

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

@click.command()
@options.commit
@click.option(
    "--analysis", "-a",
    type=click.INT,
    required=True,
    help="Analysis primary key.",
    metavar="<integer>")
@click.option(
    "--sampcond", "-sc",
    help="Sample condition file. Should contain fields: source_name and source_condition.",
    required=True,
    metavar="<string>")
@click.option(
    "--gtf", "-g",
    type=click.Choice(["ensembl", "gencode"]),
    help="Annotation to use.",
    required=True,
    metavar="<string>")
@click.option(
    "--group", "-gr",
    help="Groups file path.",
    required=True,
    metavar="<string>")
@click.option(
    "--npca", "-np",
    help="Number of genes in PCA. DEFAULT=10",
    required=False,
    default=10,
    metavar="<integer>")
@click.option(
    "--nheatmap", "-nh",
    help="Number of genes in Heatmap. DEFAULT=50",
    required=False,
    default=50,
    metavar="<integer>")
def rnaseqdownstream(analysis, sampcond, gtf, group, npca, nheatmap, commit):
    """Run Differential Test analysis on readcounts."""
    pk = "A%07dL" % analysis
    analysesdir = constants.LEUKDC_ROOT + '/analyses'
    samples = analysesdir + "/%s" % pk + "/results/samples.txt"
    analysesdir = analysesdir + "/%s" % pk + "/results/output"
    click.secho("\nOutput directory:\n", fg="blue")
    click.echo("\t%s\n" % analysesdir)
    script = os.path.join(constants.LEUKTOOLS_ROOT, 'leuktools/nested/run/r', 'DGE_analysis.R')

    cmdstring = "Rscript {0} {1} {2} {3} {4} {5} {6} {7}".format(script, analysesdir, sampcond, gtf, group, npca, nheatmap, samples)
    click.secho("\nCommand:\n", fg="blue")
    click.echo("\t%s\n" % cmdstring)

    if commit:
        os.system(cmdstring)

    if not commit:
        click.secho("Add --commit to submit.\n", fg="green", blink=True)
        raise click.Abort

if __name__ == "__main__":
    rnaseqdownstream()
