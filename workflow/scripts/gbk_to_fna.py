from Bio import SeqIO

input_file = str(snakemake.input)
output_file = str(snakemake.output)

fmt_from = str(snakemake.params['from'])
fmt_to = str(snakemake.params['to'])

records = SeqIO.parse(input_file, fmt_from)
with open(output_file, 'w') as fh_out:
    for record in records:
        SeqIO.write(record, fh_out, fmt_to)
