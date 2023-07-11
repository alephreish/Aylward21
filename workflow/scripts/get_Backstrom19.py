
from Bio import GenBank
import re

input_file = str(snakemake.input)
output_file = str(snakemake.output)

genome = snakemake.params['genome']
search = genome['search']
search_re = re.compile(search + r'\b')

with open(input_file) as fh:
    with open(output_file, 'w') as fh_out:
        records = GenBank.parse(fh)
        for record in records:
            if search_re.match(record.source):
                fh_out.write(str(record) + '\n')
