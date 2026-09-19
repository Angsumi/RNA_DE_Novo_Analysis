import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from matplotlib_venn import venn3
    VENN_AVAILABLE = True
except ImportError:
    VENN_AVAILABLE = False

from pathlib import Path

# Resolve repository root directory
REPO_ROOT = Path(__file__).resolve().parents[2]

de_dir = REPO_ROOT / 'data' / 'differential_expression' / 'gene_level'
report_file = REPO_ROOT / 'data' / 'annotation' / 'trinotate_DE_annotation_report.xls'
fig_dir = REPO_ROOT / 'results' / 'figures'
os.makedirs(fig_dir, exist_ok=True)

# Expression matrix: check external matrix first, fallback to repository table
tmm_file = Path('/home/angsuman/rna_pipeline/trinity_output/trinity_matrix.gene.TMM.EXPR.matrix')
if not tmm_file.exists():
    tmm_file = REPO_ROOT / 'results' / 'tables' / 'tsv' / 'de_genes_expression_matrix.tsv'

# 1. Generate Venn Diagram
print("Generating Venn Diagram...")
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
        # filter FDR < 0.05 (the gene IDs are the index)
        sig_genes = set(df[df['FDR'] < 0.05].index)
        de_sets[comp] = sig_genes
    else:
        de_sets[comp] = set()

if VENN_AVAILABLE:
    plt.figure(figsize=(8, 8))
    venn3(subsets=[de_sets['AZ1_vs_AZ2'], de_sets['AZ1_vs_AZ3'], de_sets['AZ2_vs_AZ3']],
          set_labels=('AZ1 vs AZ2', 'AZ1 vs AZ3', 'AZ2 vs AZ3'))
    plt.title('Differentially Expressed Genes Overlap (FDR < 0.05)', fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(fig_dir, 'de_overlap_venn.png'), dpi=300, bbox_inches='tight')
    plt.close()
else:
    print("Warning: 'matplotlib_venn' not found. Skipping Venn diagram generation.")

# 2. Generate Heatmap
print("Generating Heatmap...")
all_de_genes = de_sets['AZ1_vs_AZ2'].union(de_sets['AZ1_vs_AZ3']).union(de_sets['AZ2_vs_AZ3'])
if os.path.exists(tmm_file) and len(all_de_genes) > 0:
    tmm_df = pd.read_csv(tmm_file, sep='\t', index_col=0)
    
    # filter for DE genes
    de_tmm = tmm_df.loc[tmm_df.index.intersection(all_de_genes)]
    
    if len(de_tmm) > 0:
        # log2 transform
        log_de_tmm = np.log2(de_tmm + 1)
        
        # Plot clustermap
        cg = sns.clustermap(log_de_tmm, cmap='viridis', metric='euclidean', method='average',
                            col_cluster=False, figsize=(8, 10), 
                            yticklabels=False, cbar_kws={'label': 'log2(TMM + 1)'})
        cg.ax_heatmap.set_title(f'Expression Heatmap of DE Genes (n={len(de_tmm)})', fontsize=14, fontweight='bold', pad=20)
        plt.savefig(os.path.join(fig_dir, 'de_heatmap.png'), dpi=300, bbox_inches='tight')
        plt.close()

# 3. Parse Trinotate Report for Pfam and GO
print("Generating Pfam & GO figures...")
if os.path.exists(report_file):
    try:
        rep_df = pd.read_csv(report_file, sep='\t', header=0, low_memory=False)
    except Exception as e:
        print("Error reading xls:", e)
        rep_df = None
        
    if rep_df is not None:
        # Parse Pfam
        pfams = []
        pfam_col = [col for col in rep_df.columns if 'pfam' in col.lower() or 'hmmer' in col.lower()]
        if pfam_col:
            col_name = pfam_col[0]
            for val in rep_df[col_name].dropna():
                for line in str(val).split('`'):
                    if '^' in line:
                        parts = line.split('^')
                        if len(parts) >= 3:
                            pfams.append(parts[2])
                            
        if pfams:
            pfam_counts = pd.Series(pfams).value_counts().head(15)
            plt.figure(figsize=(10, 6))
            sns.barplot(x=pfam_counts.values, y=pfam_counts.index, palette='crest')
            plt.title('Top 15 Most Abundant Pfam Domains in DE Transcripts', fontsize=14, fontweight='bold')
            plt.xlabel('Number of Transcripts')
            plt.ylabel('Pfam Domain')
            plt.savefig(os.path.join(fig_dir, 'top_pfam_domains.png'), dpi=300, bbox_inches='tight')
            plt.close()

        # Parse GO terms
        go_terms = []
        go_col = [col for col in rep_df.columns if 'go' in col.lower() or 'ontology' in col.lower()]
        if go_col:
            col_name = go_col[0]
            for val in rep_df[col_name].dropna():
                for term in str(val).split('`'):
                    if '^' in term:
                        parts = term.split('^')
                        if len(parts) >= 3:
                            go_terms.append(parts[2])
                            
        if go_terms:
            go_counts = pd.Series(go_terms).value_counts().head(15)
            plt.figure(figsize=(10, 6))
            sns.barplot(x=go_counts.values, y=go_counts.index, palette='mako')
            plt.title('Top 15 Most Abundant GO Terms in DE Transcripts', fontsize=14, fontweight='bold')
            plt.xlabel('Number of Annotations')
            plt.ylabel('GO Term')
            plt.savefig(os.path.join(fig_dir, 'top_go_categories.png'), dpi=300, bbox_inches='tight')
            plt.close()

# 4. Generate BUSCO Plot
busco_summary = REPO_ROOT / 'data' / 'qc' / 'busco_results' / 'short_summary.specific.eukaryota_odb12.2.busco_results.txt'
if os.path.exists(busco_summary):
    with open(busco_summary, 'r') as f:
        content = f.read()
    
    # Extract using regex (Format: C:68.0%[S:41.6%,D:26.4%],F:28.8%,M:3.2%,n:125)
    match = re.search(r'C:([\d\.]+)%\[S:([\d\.]+)%,D:([\d\.]+)%\],F:([\d\.]+)%,M:([\d\.]+)%,n:(\d+)', content)
    if match:
        c_pct, s_pct, d_pct, f_pct, m_pct, n_val = map(float, match.groups())
        
        # Plot stacked bar chart
        fig, ax = plt.subplots(figsize=(10, 2))
        
        # Stacked bar values
        ax.barh([0], [s_pct], label=f'Complete Single-copy (S: {s_pct}%)', color='#3498db')
        ax.barh([0], [d_pct], left=[s_pct], label=f'Complete Duplicated (D: {d_pct}%)', color='#2c3e50')
        ax.barh([0], [f_pct], left=[s_pct + d_pct], label=f'Fragmented (F: {f_pct}%)', color='#e67e22')
        ax.barh([0], [m_pct], left=[s_pct + d_pct + f_pct], label=f'Missing (M: {m_pct}%)', color='#e74c3c')
        
        ax.set_yticks([])
        ax.set_xlim(0, 100)
        ax.set_xlabel('Percentage (%)', fontsize=12)
        ax.set_title(f'BUSCO Completeness Assessment (n={int(n_val)} eukaryota_odb10)', fontsize=14, fontweight='bold', pad=15)
        
        # Legend below
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=4, frameon=False, fontsize=10)
        
        plt.savefig(os.path.join(fig_dir, 'busco_completeness.png'), dpi=300, bbox_inches='tight')
        plt.close()

print("All figures generated successfully!")
