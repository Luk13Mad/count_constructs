import click


@click.command()
@click.option("--forward","-f",type = click.File(mode = "r"), required = True, help = "FASTQ file with forward reads.")
@click.option("--reverse","-r",type = click.File(mode = "r"), required = True, help = "FASTQ file with reverse reads.")
@click.option("--expected","-e",type = click.File(mode = "r"), required = True, help = "TSV file with expected constructs.")
@click.option("--samplefolder",type = click.Path(exists = True), required = True, help = "Absolute path to your sample folder.")
@click.option("--startspacer1","-s1",type = click.INT, required = True, help = "Start default extraction window first spacer.")
@click.option("--startspacer2","-s2",type = click.INT, required = True, help = "Start default extraction window second spacer.")
@click.option("--startspacer3","-s3",type = click.INT, required = True, help = "Start default extraction window third spacer.")
@click.option("--endspacer1","-e1",type = click.INT, required = True, help = "End default extraction window first spacer.")
@click.option("--endspacer2","-e2",type = click.INT, required = True, help = "End default extraction window second spacer.")
@click.option("--endspacer3","-e3",type = click.INT, required = True, help = "End default extraction window third spacer.")
@click.option("--resfile", default = "count_results.tsv", help = "Name of results file. STDOUT prints results to standard out.")
@click.option("--workers","-w",type = click.INT,default = 1,help = "Max number of workers to be used.")
@click.option("--batchsize",type = click.INT, default = 5000, help = "Batchsize for processing FASTQ pairs.")
def setup_3D(expected,samplefolder,forward,reverse,batchsize,workers,resfile):
    print("Called 3D subcommand")
    print("Not implemented yet.")