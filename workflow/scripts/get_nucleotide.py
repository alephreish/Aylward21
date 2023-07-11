from Bio import Entrez
import re
from shutil import copyfileobj

Entrez.email = str(snakemake.config['email'])
output_file = str(snakemake.output)
genome = snakemake.params['genome']

def expand_interval(interval):
    match = re.match(r'(.+?)(\d+)\.(\d+)- ?(.+?)(\d+)\.(\d+)', interval)
    if not match: return []
    prefix1, start, suffix1, prefix2, end, suffix2 = match.groups()
    num_pos = len(start)
    expanded_list = []
    for i in range(int(start), int(end) + 1):
        expanded_list.append(f"{prefix1}{str(i).zfill(num_pos)}.{suffix1}")
    return expanded_list

ids = expand_interval(genome['search'])
assert len(ids) < 100, "Way too many ids: {ids}".format(ids = ', '.join(ids))
handle = Entrez.efetch(db = "nucleotide", id = ids, rettype = "gb", retmode = "text")

with open(output_file, 'w') as fh:
    copyfileobj(handle, fh)
