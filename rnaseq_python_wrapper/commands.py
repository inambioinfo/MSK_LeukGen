
""":mod:`~leuktools.nested.run` module."""

# third-party
import click

# local
from .cgp import cgp
from .cgptest import cgptest
from .itd import itd
from .facets import facets
from .excavator import excavator
from .varcorrect import varcorrect
from .rnaseq import rnaseq
from .rnaseqdownstream import rnaseqdownstream

@click.group()
def run():
    """Run leukgen pipelines."""
    pass

run.add_command(cgp)
run.add_command(cgptest)
run.add_command(itd)
run.add_command(facets)
run.add_command(excavator)
run.add_command(varcorrect)
run.add_command(rnaseq)
run.add_command(rnaseqdownstream)

if __name__ == "__main__":
    run()
