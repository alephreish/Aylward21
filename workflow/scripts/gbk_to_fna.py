from Bio import SeqIO

input_file = str(snakemake.input)
output_file = str(snakemake.output)
g = snakemake.params['genome']

g["Fam"] =  g["Family"] if g["Family_AKA"] == "-" else g["Family_AKA"]
g_id = g["genome_id"]

records = SeqIO.parse(input_file, 'genbank')
with open(output_file, 'w') as fh_out:
    seq_num = 0
    for record in records:
        seq_num += 1
        taxonomy = f'p__{g["Phylum"]};c__{g["Class"]};o__{g["Order"]};f__{g["Fam"]};g__{g["Genus"]}'
        record.description = f'{record.id} genome:{g_id} common_name:{g["common_name"]} taxonomy:{taxonomy} fam_num:{g["Family"]} approach:{g["Sequencing-approach"]} source:{g["Source"]} source_detail:{g["Source-Detail"]}'
        record.id = f'{g_id}_{seq_num}'
        SeqIO.write(record, fh_out, 'fasta')
