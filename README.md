Aylward21
=========

This is a command-line tool to create a taxonomically annotated fasta and blast database of NCLDV genomes based on [Aylward21](https://doi.org/10.1371/journal.pbio.3001430).

You need [`snakemake`](https://snakemake.readthedocs.io/en/stable/) to run the pipeline, [GVMAGs database](https://genome.jgi.doe.gov/portal/GVMAGs/GVMAGs.home.html) unpacked in `sources/` and your email specified as a config parameter. The workflow makes use of a resource called `ncbi` that limits the number of the parallel requests to NCBI.

There is a convenience script `run.sh` which takes care of the parameters. The email can be specified in a config file (`config.yaml` by default) as: `email: <my@email.com>). Any arguments to `run.sh` are passed to `snakemake`.
