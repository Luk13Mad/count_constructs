import click
from .screen_2D.loop_2D import setup_2D
from .screen_3D.loop_3D import setup_3D

@click.group(help = "Command line tool for counting FASTQ files from 2D/3D CRISPR screen.")
def cli():
    pass

cli.add_command(setup_2D, name = "2D")
cli.add_command(setup_3D, name = "3D")


if __name__ == "__main__":
    cli()