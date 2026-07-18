import os
import re
import pandas as pd
import numpy as np

# Set directories
base_dir = '/home/angsuman/rna_pipeline/trinity_output'
de_dir = os.path.join(base_dir, 'edgeR_gene_dir')
report_file = os.path.join(base_dir, 'trinotate_DE_annotation_report.xls')
gene_tmm_file = '/home/angsuman/rna_pipeline/trinity_output/trinity_matrix.gene.TMM.EXPR.matrix'
busco_summary = '/home/angsuman/rna_pipeline/trinity_output/busco_results/short_summary.specific.eukaryota_odb12.2.busco_results.txt'

tables_dir = os.path.join(base_dir, 'tables')
os.makedirs(tables_dir, exist_ok=True)

# 1. Venn Diagram overlap counts csv
print("Saving Venn overlaps...")
de_files = {
    'AZ1_vs_AZ2': 'trinity_matrix.gene.counts.matrix.AZ1_vs_AZ2.edgeR.DE_results',
    'AZ1_vs_AZ3': 'trinity_matrix.gene.counts.matrix.AZ1_vs_AZ3.edgeR.DE_results',
    'AZ2_vs_AZ3': 'trinity_matrix.gene.counts.matrix.AZ2_vs_AZ3.edgeR.DE_results'
}
de_sets = {}
for comp, fname in de_files.items():
    fpath = os.path.join(de_dir, fname)
    if os.path.exists(fpath):
        df = pd.read_csv(fpath, sep='\t')
        de_sets[comp] = set(df[df['FDR'] < 0.05].index)
    else:
        de_sets[comp] = set()

# Calculate specific overlaps
a12 = de_sets['AZ1_vs_AZ2']
a13 = de_sets['AZ1_vs_AZ3']
a23 = de_sets['AZ2_vs_AZ3']

venn_data = {
    'Comparison_Group': [
        'AZ1_vs_AZ2_Only', 'AZ1_vs_AZ3_Only', 'AZ2_vs_AZ3_Only',
        'Shared_AZ1_vs_AZ2_AND_AZ1_vs_AZ3', 'Shared_AZ1_vs_AZ2_AND_AZ2_vs_AZ3', 'Shared_AZ1_vs_AZ3_AND_AZ2_vs_AZ3',
        'Shared_All_Three', 'Total_Unique_DE_Genes'
    ],
    'Gene_Count': [
        len(a12 - a13 - a23), len(a13 - a12 - a23), len(a23 - a12 - a13),
        len((a12 & a13) - a23), len((a12 & a23) - a13), len((a13 & a23) - a12),
        len(a12 & a13 & a23), len(a12 | a13 | a23)
    ]
}
pd.DataFrame(venn_data).to_csv(os.path.join(tables_dir, 'de_overlap_counts.tsv'), sep='\t', index=False)

# 2. Clustered Heatmap raw TMM values csv
print("Saving Heatmap raw data...")
all_de_genes = a12.union(a13).union(a23)
if os.path.exists(gene_tmm_file) and len(all_de_genes) > 0:
    tmm_df = pd.read_csv(gene_tmm_file, sep='\t', index_col=0)
    de_tmm = tmm_df.loc[tmm_df.index.intersection(all_de_genes)]
    de_tmm.to_csv(os.path.join(tables_dir, 'de_genes_expression_matrix.tsv'), sep='\t')

# 3. Pfam and GO top categories csvs
print("Saving Pfam & GO Top lists...")
if os.path.exists(report_file):
    rep_df = pd.read_csv(report_file, sep='\t', header=0, low_memory=False)
    
    # Pfam domains
    pfams = []
    pfam_col = [col for col in rep_df.columns if 'pfam' in col.lower() or 'hmmer' in col.lower()]
    if pfam_col:
        for val in rep_df[pfam_col[0]].dropna():
            for line in str(val).split('`'):
                if '^' in line:
                    parts = line.split('^')
                    if len(parts) >= 3:
                        pfams.append((parts[0], parts[2])) # ID, name
                        
    if pfams:
        pfam_df = pd.DataFrame(pfams, columns=['Pfam_ID', 'Domain_Description'])
        pfam_counts = pfam_df.value_counts().head(30).reset_index(name='Counts')
        pfam_counts.to_csv(os.path.join(tables_dir, 'top_pfam_domains.tsv'), sep='\t', index=False)

    # GO categories
    go_terms = []
    go_col = [col for col in rep_df.columns if 'go' in col.lower() or 'ontology' in col.lower()]
    if go_col:
        for val in rep_df[go_col[0]].dropna():
            for term in str(val).split('`'):
                if '^' in term:
                    parts = term.split('^')
                    if len(parts) >= 3:
                        go_terms.append((parts[0], parts[1], parts[2])) # ID, Ontology, Term
                        
    if go_terms:
        go_df = pd.DataFrame(go_terms, columns=['GO_ID', 'Ontology', 'Term_Description'])
        go_counts = go_df.value_counts().head(30).reset_index(name='Counts')
        go_counts.to_csv(os.path.join(tables_dir, 'top_go_categories.tsv'), sep='\t', index=False)

# 4. BUSCO completeness csv
print("Saving BUSCO percentages...")
if os.path.exists(busco_summary):
    with open(busco_summary, 'r') as f:
        content = f.read()
    match = re.search(r'C:([\d\.]+)%\[S:([\d\.]+)%,D:([\d\.]+)%\],F:([\d\.]+)%,M:([\d\.]+)%,n:(\d+)', content)
    if match:
        c_pct, s_pct, d_pct, f_pct, m_pct, n_val = map(float, match.groups())
        busco_data = {
            'Category': ['Complete_Single_Copy', 'Complete_Duplicated', 'Fragmented', 'Missing'],
            'Orthologs_Count': [int(round(s_pct * n_val / 100)), int(round(d_pct * n_val / 100)), int(round(f_pct * n_val / 100)), int(round(m_pct * n_val / 100))],
            'Percentage': [s_pct, d_pct, f_pct, m_pct]
        }
        pd.DataFrame(busco_data).to_csv(os.path.join(tables_dir, 'busco_completeness_percentages.tsv'), sep='\t', index=False)

# 5. Top 6 Candidate genes expressions csv
print("Saving top candidate genes expression profiles...")
de_files = [
    'trinity_matrix.gene.counts.matrix.AZ1_vs_AZ2.edgeR.DE_results',
    'trinity_matrix.gene.counts.matrix.AZ1_vs_AZ3.edgeR.DE_results',
    'trinity_matrix.gene.counts.matrix.AZ2_vs_AZ3.edgeR.DE_results'
]
lowest_fdr_genes = []
for fname in de_files:
    fpath = os.path.join(de_dir, fname)
    if os.path.exists(fpath):
        df = pd.read_csv(fpath, sep='\t')
        lowest_fdr_genes.extend(df.sort_values('FDR').head(10).index.tolist())
top_6_genes = []
for g in lowest_fdr_genes:
    if g not in top_6_genes:
        top_6_genes.append(g)
    if len(top_6_genes) == 6:
        break

if len(top_6_genes) == 6 and os.path.exists(gene_tmm_file):
    tmm_df = pd.read_csv(gene_tmm_file, sep='\t', index_col=0)
    top_6_tmm = tmm_df.loc[top_6_genes]
    top_6_tmm.to_csv(os.path.join(tables_dir, 'top_candidate_genes_expression.tsv'), sep='\t')

print("All CSV raw data tables generated successfully!")
