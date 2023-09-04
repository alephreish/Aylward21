
from Bio import SeqIO
import re

input_file = str(snakemake.input)
output_file = str(snakemake.output)

genome = snakemake.params['genome']
search = genome['search']
search_re = re.compile(r'\b' + search + r'\b')

records = SeqIO.parse(input_file, 'genbank')
with open(output_file, 'w') as out:
    for record in records:
        if search_re.search(record.description):
            SeqIO.write(record, out, 'genbank')
