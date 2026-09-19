import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

de_dir = REPO_ROOT / 'data' / 'differential_expression' / 'gene_level'
trans_map_file = Path('/home/angsuman/rna_pipeline/trinity_output/Trinity.fasta.gene_trans_map')
fasta_file = Path('/home/angsuman/rna_pipeline/trinity_output/Trinity.fasta')
output_fasta = Path('/home/angsuman/rna_pipeline/trinity_output/Trinity_DE_subset.fasta')

# 1. Collect DE gene IDs with FDR < 0.05
de_genes = set()
for f in os.listdir(de_dir):
    if f.endswith('.DE_results'):
        path = os.path.join(de_dir, f)
        with open(path, 'r') as fh:
            header = fh.readline()
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    gene_id = parts[0]
                    try:
                        fdr = float(parts[6])
                        if fdr < 0.05:
                            de_genes.add(gene_id)
                    except ValueError:
                        continue

print(f'Collected {len(de_genes)} unique DE genes.')

# 2. Map gene IDs to transcript IDs
gene_to_trans = {}
with open(trans_map_file, 'r') as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            gene, trans = parts[0], parts[1]
            if gene in de_genes:
                if gene not in gene_to_trans:
                    gene_to_trans[gene] = []
                gene_to_trans[gene].append(trans)

de_transcripts = set()
for trans_list in gene_to_trans.values():
    de_transcripts.update(trans_list)

print(f'Mapped to {len(de_transcripts)} unique transcript sequences.')

# 3. Extract FASTA records without Biopython
with open(fasta_file, 'r') as infile, open(output_fasta, 'w') as outfile:
    current_seq = []
    current_id = None
    for line in infile:
        if line.startswith('>'):
            if current_id and current_id in de_transcripts:
                outfile.write('>' + current_id + '\n' + ''.join(current_seq))
            current_id = line[1:].strip().split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    # Write the last sequence
    if current_id and current_id in de_transcripts:
        outfile.write('>' + current_id + '\n' + ''.join(current_seq))

print('Subset FASTA file written successfully.')
