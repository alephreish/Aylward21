from Bio import Entrez
import urllib.request
from os.path import basename
import gzip
import shutil

Entrez.email = str(snakemake.config['email'])
output_file = str(snakemake.output)
genome = snakemake.params['genome']
acc = genome['search']

def decompress_gzip_file(url, local_file_path):
    file = basename(url) + '_genomic.gbff.gz'
    response = urllib.request.urlopen(f"{url}/{file}".replace("ftp://", "https://"))
    with gzip.GzipFile(fileobj = response) as gzipped_file:
        with open(local_file_path, 'wb') as decompressed_file:
            shutil.copyfileobj(gzipped_file, decompressed_file)

handle = Entrez.esearch(db = "assembly", term = f"{acc}[ASAC]")
record = Entrez.read(handle)
assert len(record['IdList']) > 0, f"Nothing found for {acc}"
id = record['IdList'][0]
esummary_handle = Entrez.esummary(db = "assembly", id = id, report = "full")
esummary_record = Entrez.read(esummary_handle)
url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
decompress_gzip_file(url, output_file)
