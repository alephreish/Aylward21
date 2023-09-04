from Bio import SeqIO
from os import path, listdir

output_file = str(snakemake.output)
genome = snakemake.params['genome']

input_file = None
found = False
for dirname1 in snakemake.input:
    for dirname2 in listdir(path.join(dirname1)):
        filename = path.join(dirname1, dirname2, genome['search'] + '.fna')
        if path.exists(filename):
            input_file = filename
            found = True
            break
    if found: break
assert found, "File not found: {filename}".format(filename = genome['search'] + '.fna')
records = SeqIO.parse(input_file, 'fasta')
with open(output_file, 'w') as fh_out:
    for record in records:
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write(record, fh_out, 'genbank')
