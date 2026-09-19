# De Novo RNA-Seq Transcriptome Assembly and Downstream Analysis

This repository contains the complete processed results, publication-ready scientific figures, parsed CSV/TSV tables, and reproducible scripts from the de novo transcriptome assembly, differential expression, and functional annotation workflow across three conditions/timepoints (**AZ1**, **AZ2**, **AZ3**).

---

## 📂 Repository Structure

```
RNA_DE_Novo_Analysis/
├── README.md                          # Master project documentation
├── requirements.txt                   # Python package dependencies
├── environment.yml                    # Conda environment definition
│
├── scripts/                           # Modular workflow scripts organized by stage
│   ├── 01_assembly/                   # Trinity de novo assembly scripts
│   │   ├── run_trinity_norm.sh
│   │   ├── run_trinity_norm_12cpu.sh
│   │   └── run_trinity_norm_12cpu_fix.sh
│   ├── 02_quantification/            # Salmon pseudoalignment & abundance estimation
│   │   └── run_salmon_manual.sh
│   ├── 03_diff_expression/          # edgeR DE analysis & FASTA subsetting
│   │   ├── run_de_analysis_no_replicates.sh
│   │   └── extract_de_sequences.py
│   ├── 04_annotation/               # Functional annotation (Trinotate, TransDecoder, DIAMOND, HMMER)
│   │   ├── run_trinotate.sh
│   │   └── run_trinotate_subset.sh
│   ├── 05_downstream_analysis/      # Data mining, publication plotting & table export
│   │   ├── additional_analyses.py    # PCA, KEGG pathways, TF classification, isoform switches
│   │   ├── generate_figures.py       # High-res publication figure generation
│   │   └── generate_csv_data.py      # Dual TSV and CSV tabular compilation
│   └── pipeline_orchestration/       # Monitoring and multi-step pipeline runners
│       ├── run_downstream.sh
│       └── run_downstream_all.sh
│
├── data/                             # Processed data, statistics & annotations
│   ├── qc/                           # Assembly quality & completeness assessment
│   │   ├── TrinityStats.txt          # N50, GC content, length distribution
│   │   └── busco_results/            # BUSCO v6.1.0 (eukaryota_odb12.2) benchmarking summaries
│   │       ├── short_summary.specific.eukaryota_odb12.2.busco_results.json
│   │       └── short_summary.specific.eukaryota_odb12.2.busco_results.txt
│   ├── annotation/                   # Comprehensive functional annotation outputs
│   │   └── trinotate_DE_annotation_report.xls
│   └── differential_expression/      # Count matrices, DE results & volcano plots
│       ├── gene_level/               # Gene-level edgeR comparisons (.DE_results, .pdf, .Rscript)
│       └── isoform_level/            # Transcript-level edgeR comparisons
│
└── results/                          # Final deliverables
    ├── figures/                      # High-resolution (300 DPI) publication figures (.png)
    │   ├── busco_completeness.png    # BUSCO ortholog completeness assessment
    │   ├── de_overlap_venn.png       # 3-way Venn diagram of DE gene overlaps
    │   ├── de_heatmap.png            # Hierarchical clustering heatmap of 432 DE genes
    │   ├── pca_plot.png              # Principal Component Analysis (PCA) of samples
    │   ├── candidate_genes_expression.png # Expression profiles of top 6 candidate genes
    │   ├── top_go_categories.png     # Top 15 most frequent Gene Ontology (GO) terms
    │   └── top_pfam_domains.png      # Top 15 most frequent Pfam domains
    └── tables/                       # Processed tabular deliverables
        ├── csv/                      # Comma-Separated Values (for Excel, R, Python)
        └── tsv/                      # Tab-Separated Values (for UNIX pipelines, Prism)
```

---

## 🧬 Summary of Completed Analyses

### 1. De Novo Assembly (Trinity)
* **Total Genes**: 2,682,805
* **Total Transcripts**: 2,944,151
* **Contig N50**: 363 bp
* **Total Assembled Bases**: 1.07 Gbp
* **Assembly Completeness (BUSCO)**: 68.0% Complete (85/125 orthologs), 28.8% Fragmented (36/125), 3.2% Missing (4/125).

### 2. Abundance Estimation (Salmon)
* Abundance quantification mapped raw reads back to assembled contigs:
  * **AZ1**: 80.37% alignment rate (10.75M mapped reads)
  * **AZ2**: 79.47% alignment rate (10.78M mapped reads)
  * **AZ3**: 80.28% alignment rate (6.73M mapped reads)

### 3. Differential Expression (edgeR without replicates)
* Pairwise DE genes identified at FDR < 0.05 (dispersion fixed at 0.1):
  * **AZ1 vs AZ2**: 993 genes (441 with |log2FC| > 2)
  * **AZ1 vs AZ3**: 1,103 genes (417 with |log2FC| > 2)
  * **AZ2 vs AZ3**: 839 genes (453 with |log2FC| > 2)
* A union of **432 unique differentially expressed genes** (mapping to 2,840 transcripts) was extracted for downstream annotation.

### 4. Functional Annotation (Trinotate)
* Functional profiling of the 2,840 DE transcripts:
  * Open Reading Frames (ORFs) predicted by **TransDecoder**.
  * Homologous proteins annotated by **DIAMOND** blastp (against SwissProt).
  * Functional domains annotated by **HMMER** (against Pfam).
  * Non-coding RNAs scanned by **Infernal (cmscan)** (against Rfam).
* Compiled into the comprehensive report sheet `data/annotation/trinotate_DE_annotation_report.xls`.

---

## 🚀 Quickstart & Reproducibility

### 1. Environment Setup

Using conda:
```bash
conda env create -f environment.yml
conda activate rna_downstream
```

Or using pip:
```bash
pip install -r requirements.txt
```

### 2. Regenerate Figures and Tables

All downstream analysis scripts use repository-relative paths and run standalone:

```bash
# 1. Regenerate all publication figures (.png)
python scripts/05_downstream_analysis/generate_figures.py

# 2. Re-export all TSV and CSV tables
python scripts/05_downstream_analysis/generate_csv_data.py

# 3. Re-run PCA, TF classification, KEGG pathways, and isoform switching analysis
python scripts/05_downstream_analysis/additional_analyses.py
```

### 3. Accessing Results
- Figures are stored in [`results/figures/`](results/figures/).
- Tables are stored in [`results/tables/csv/`](results/tables/csv/) and [`results/tables/tsv/`](results/tables/tsv/).
