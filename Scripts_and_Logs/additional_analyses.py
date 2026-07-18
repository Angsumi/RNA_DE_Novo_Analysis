import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# Set directories
base_dir = '/home/angsuman/rna_pipeline/trinity_output'
de_dir = os.path.join(base_dir, 'edgeR_gene_dir')
report_file = os.path.join(base_dir, 'trinotate_DE_annotation_report.xls')
gene_tmm_file = '/home/angsuman/rna_pipeline/trinity_output/trinity_matrix.gene.TMM.EXPR.matrix'
iso_tmm_file = '/home/angsuman/rna_pipeline/trinity_output/trinity_matrix.isoform.TMM.EXPR.matrix'
gene_trans_map_file = '/home/angsuman/rna_pipeline/trinity_output/Trinity_DE_subset.fasta.gene_trans_map'

fig_dir = os.path.join(base_dir, 'figures')
tables_dir = os.path.join(base_dir, 'tables')
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)

# 1. PCA Scatter Plot
print("1. Running Principal Component Analysis (PCA)...")
if os.path.exists(gene_tmm_file):
    tmm_df = pd.read_csv(gene_tmm_file, sep='\t', index_col=0)
    # Log transform TMM matrix
    log_tmm = np.log2(tmm_df + 1)
    
    # Run PCA on the samples (transpose to have samples as rows)
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(log_tmm.T)
    var_exp = pca.explained_variance_ratio_ * 100
    
    pca_df = pd.DataFrame(pca_results, columns=['PC1', 'PC2'], index=log_tmm.columns)
    
    plt.figure(figsize=(7, 6))
    sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=150, color='#9b59b6', edgecolor='black', linewidth=1.5)
    
    # Label samples
    for sample in pca_df.index:
        plt.text(pca_df.loc[sample, 'PC1'] + 0.1, pca_df.loc[sample, 'PC2'] + 0.1, sample, 
                 fontsize=12, fontweight='bold')
                 
    plt.title('Principal Component Analysis (PCA) of Samples', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel(f'PC1 ({var_exp[0]:.2f}% Variance)', fontsize=12)
    plt.ylabel(f'PC2 ({var_exp[1]:.2f}% Variance)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig(os.path.join(fig_dir, 'pca_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save PCA coordinates table
    pca_df.to_csv(os.path.join(tables_dir, 'pca_coordinates.tsv'), sep='\t')

# 2. KEGG Pathway Analysis
print("2. Parsing KEGG Pathways...")
if os.path.exists(report_file):
    rep_df = pd.read_csv(report_file, sep='\t', header=0, low_memory=False)
    
    kegg_pathways = []
    kegg_col = [col for col in rep_df.columns if 'kegg' in col.lower()]
    if kegg_col:
        col_name = kegg_col[0]
        # Match pattern like path:ko04010 or similar
        for val in rep_df[col_name].dropna():
            for part in str(val).split('`'):
                match = re.search(r'path:([a-zA-Z0-9]+)', part)
                if match:
                    kegg_pathways.append(match.group(1))
                    
    # Map common KO codes to names for readability
    kegg_names = {
        'ko04010': 'MAPK signaling pathway',
        'ko00190': 'Oxidative phosphorylation',
        'ko03010': 'Ribosome',
        'ko00010': 'Glycolysis / Gluconeogenesis',
        'ko00020': 'TCA cycle',
        'ko00910': 'Nitrogen metabolism',
        'ko04141': 'Protein processing in endoplasmic reticulum',
        'ko00906': 'Carotenoid biosynthesis',
        'ko00500': 'Starch and sucrose metabolism',
        'ko04075': 'Plant hormone signal transduction',
        'ko04626': 'Plant-pathogen interaction',
        'ko00195': 'Photosynthesis',
        'ko00940': 'Phenylpropanoid biosynthesis',
        'ko00920': 'Sulfur metabolism',
        'ko04120': 'Ubiquitin mediated proteolysis'
    }
    
    if kegg_pathways:
        kegg_series = pd.Series(kegg_pathways)
        kegg_counts = kegg_series.value_counts().head(10)
        
        # Replace code with name if available
        labels = [kegg_names.get(code, code) for code in kegg_counts.index]
        
        plt.figure(figsize=(10, 5))
        sns.barplot(x=kegg_counts.values, y=labels, palette='flare')
        plt.title('Top 10 Enriched KEGG Pathways in DE Genes', fontsize=14, fontweight='bold')
        plt.xlabel('Number of Transcripts')
        plt.ylabel('KEGG Pathway')
        plt.savefig(os.path.join(fig_dir, 'top_kegg_pathways.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save table
        kegg_table = pd.DataFrame({'Pathway_Code': kegg_counts.index, 'Description': labels, 'Counts': kegg_counts.values})
        kegg_table.to_csv(os.path.join(tables_dir, 'kegg_pathways_abundance.tsv'), sep='\t', index=False)

# 3. TF Family Classification
print("3. Classifying Transcription Factors (TFs)...")
tf_map = {
    'PF00096': 'ZF-C2H2',
    'PF00010': 'bHLH',
    'PF00249': 'MYB',
    'PF02442': 'WRKY',
    'PF00170': 'bZIP',
    'PF00319': 'MADS-box',
    'PF00250': 'AP2/ERF',
    'PF01422': 'NAC',
    'PF00497': 'Homeobox',
    'PF03106': 'ARF',
    'PF00643': 'B3',
    'PF00642': 'GRAS',
    'PF01047': 'HLH',
    'PF00567': 'GATA',
    'PF00610': 'Heat Shock Factor',
    'PF01344': 'bZIP_2'
}

if os.path.exists(report_file):
    rep_df = pd.read_csv(report_file, sep='\t', header=0, low_memory=False)
    
    detected_tfs = []
    pfam_col = [col for col in rep_df.columns if 'pfam' in col.lower() or 'hmmer' in col.lower()]
    if pfam_col:
        col_name = pfam_col[0]
        for idx, row in rep_df.iterrows():
            val = row[col_name]
            gene_id = row[rep_df.columns[0]]
            if pd.notna(val):
                for line in str(val).split('`'):
                    for pfam_id, tf_family in tf_map.items():
                        if pfam_id in line:
                            detected_tfs.append({'gene_id': gene_id, 'TF_Family': tf_family, 'Pfam_Hit': line})
                            
    if detected_tfs:
        tf_df = pd.DataFrame(detected_tfs).drop_duplicates(subset=['gene_id', 'TF_Family'])
        tf_counts = tf_df['TF_Family'].value_counts()
        
        plt.figure(figsize=(9, 5))
        sns.barplot(x=tf_counts.values, y=tf_counts.index, palette='crest')
        plt.title('Distribution of DE Transcription Factor Families', fontsize=14, fontweight='bold')
        plt.xlabel('Number of Regulated Genes')
        plt.ylabel('TF Family')
        plt.savefig(os.path.join(fig_dir, 'tf_families_abundance.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save tables
        tf_df.to_csv(os.path.join(tables_dir, 'detected_tfs_list.tsv'), sep='\t', index=False)
        tf_counts.to_frame('Counts').to_csv(os.path.join(tables_dir, 'tf_families_summary.tsv'), sep='\t')

# 4. Expression Profiles of Top 6 DE Candidate Genes
print("4. Plotting top candidate genes...")
de_files = [
    'trinity_matrix.gene.counts.matrix.AZ1_vs_AZ2.edgeR.DE_results',
    'trinity_matrix.gene.counts.matrix.AZ1_vs_AZ3.edgeR.DE_results',
    'trinity_matrix.gene.counts.matrix.AZ2_vs_AZ3.edgeR.DE_results'
]

# Find lowest FDR genes overall
lowest_fdr_genes = []
for fname in de_files:
    fpath = os.path.join(de_dir, fname)
    if os.path.exists(fpath):
        df = pd.read_csv(fpath, sep='\t')
        lowest_fdr_genes.extend(df.sort_values('FDR').head(10).index.tolist())

# Select top 6 unique genes
top_6_genes = []
for g in lowest_fdr_genes:
    if g not in top_6_genes:
        top_6_genes.append(g)
    if len(top_6_genes) == 6:
        break

if len(top_6_genes) == 6 and os.path.exists(gene_tmm_file):
    tmm_df = pd.read_csv(gene_tmm_file, sep='\t', index_col=0)
    
    fig, axes = plt.subplots(3, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for i, gene in enumerate(top_6_genes):
        gene_exp = tmm_df.loc[gene]
        ax = axes[i]
        sns.barplot(x=gene_exp.index, y=gene_exp.values, ax=ax, palette='muted')
        ax.set_title(f'{gene} Expression', fontsize=11, fontweight='bold')
        ax.set_ylabel('Expression (TMM)')
        ax.set_xlabel('Samples')
        
    plt.suptitle('Expression Profiles of Top 6 Differentially Expressed Candidate Genes', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'candidate_genes_expression.png'), dpi=300, bbox_inches='tight')
    plt.close()

# 5. Isoform Switching / Alternative Splicing Analysis
print("5. Identifying potential isoform switches...")
if os.path.exists(iso_tmm_file) and os.path.exists(gene_trans_map_file):
    # Load gene-to-transcript map
    map_df = pd.read_csv(gene_trans_map_file, sep='\t', header=None, names=['gene_id', 'isoform_id'])
    
    # Load isoform TMM expressions
    iso_tmm = pd.read_csv(iso_tmm_file, sep='\t', index_col=0)
    
    # Merge expression data with gene mapping
    iso_data = map_df.merge(iso_tmm, left_on='isoform_id', right_index=True)
    
    # We want to find genes with multiple isoforms
    gene_counts = iso_data['gene_id'].value_counts()
    multi_iso_genes = gene_counts[gene_counts > 1].index.tolist()
    
    switches = []
    # Loop over multi-isoform genes and find dominant isoforms
    for gene in multi_iso_genes[:1000]:  # Limit to first 1000 multi-iso genes to save time
        sub_df = iso_data[iso_data['gene_id'] == gene]
        
        # For each sample, identify the dominant isoform ID (with highest TMM)
        dom_az1 = sub_df.loc[sub_df['AZ1'].idxmax(), 'isoform_id']
        dom_az2 = sub_df.loc[sub_df['AZ2'].idxmax(), 'isoform_id']
        dom_az3 = sub_df.loc[sub_df['AZ3'].idxmax(), 'isoform_id']
        
        # Check if the dominant isoform is different between ANY of the samples
        if dom_az1 != dom_az2 or dom_az1 != dom_az3 or dom_az2 != dom_az3:
            switches.append({
                'gene_id': gene,
                'Dominant_AZ1': dom_az1,
                'Dominant_AZ2': dom_az2,
                'Dominant_AZ3': dom_az3,
                'AZ1_Dominant_TMM': sub_df.loc[sub_df['AZ1'].idxmax(), 'AZ1'],
                'AZ2_Dominant_TMM': sub_df.loc[sub_df['AZ2'].idxmax(), 'AZ2'],
                'AZ3_Dominant_TMM': sub_df.loc[sub_df['AZ3'].idxmax(), 'AZ3']
            })
            
    if switches:
        switch_df = pd.DataFrame(switches)
        switch_df.to_csv(os.path.join(tables_dir, 'isoform_switching_candidates.tsv'), sep='\t', index=False)
        print(f"Found {len(switch_df)} genes showing dominant isoform switches between samples.")

print("Additional analyses and tables created successfully!")
