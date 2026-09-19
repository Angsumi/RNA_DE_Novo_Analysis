#!/bin/bash
export PATH=/home/angsuman/miniconda3/envs/downstream/bin:/home/angsuman/miniconda3/bin:$PATH
cd /home/angsuman/rna_pipeline/trinity_output

echo "Subset Trinotate pipeline started at $(date)" > trinotate_subset.log

# 1. Download and create databases
echo "[1/5] Preparing Trinotate boilerplate SQLite database (downloading databases)..." >> trinotate_subset.log
mkdir -p /home/angsuman/trinotate_data
Trinotate --db /home/angsuman/rna_pipeline/trinity_output/trinotate_subset.sqlite --create --trinotate_data_dir /home/angsuman/trinotate_data >> trinotate_subset.log 2>&1

# 2. Run TransDecoder
echo "[2/5] Running TransDecoder LongOrfs on DE subset..." >> trinotate_subset.log
TransDecoder.LongOrfs -t Trinity_DE_subset.fasta >> trinotate_subset.log 2>&1

echo "[2/5] Running TransDecoder Predict on DE subset..." >> trinotate_subset.log
TransDecoder.Predict -t Trinity_DE_subset.fasta >> trinotate_subset.log 2>&1

# 3. Initialize database
echo "[3/5] Creating subset gene-to-transcript map and initializing database..." >> trinotate_subset.log
# Filter gene trans map to only include DE transcripts
grep -F -f <(grep '>' Trinity_DE_subset.fasta | tr -d '>') Trinity.fasta.gene_trans_map > Trinity_DE_subset.fasta.gene_trans_map

Trinotate --db trinotate_subset.sqlite --init   --gene_trans_map Trinity_DE_subset.fasta.gene_trans_map   --transcript_fasta Trinity_DE_subset.fasta   --transdecoder_pep Trinity_DE_subset.fasta.transdecoder.pep >> trinotate_subset.log 2>&1

# 4. Run computes with DIAMOND
echo "[4/5] Running DIAMOND and Pfam alignments on DE subset (using 12 CPUs)..." >> trinotate_subset.log
Trinotate --db trinotate_subset.sqlite --run ALL --CPU 12 --use_diamond   --transcript_fasta Trinity_DE_subset.fasta   --transdecoder_pep Trinity_DE_subset.fasta.transdecoder.pep   --trinotate_data_dir /home/angsuman/trinotate_data >> trinotate_subset.log 2>&1

# 5. Report generation
echo "[5/5] Generating final functional annotation report..." >> trinotate_subset.log
Trinotate --db trinotate_subset.sqlite --report > trinotate_DE_annotation_report.xls 2>> trinotate_subset.log

echo "Subset Trinotate pipeline completed successfully at $(date)" >> trinotate_subset.log
