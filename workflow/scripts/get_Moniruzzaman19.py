from Bio import SeqIO
from os import path

input_dir = str(snakemake.input)
output_file = str(snakemake.output)
genome = snakemake.params['genome']

input_file = path.join(input_dir, genome['search'] + '.fa')
records = SeqIO.parse(input_file, 'fasta')
with open(output_file, 'w') as fh_out:
    for record in records:
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write(record, fh_out, 'genbank')
