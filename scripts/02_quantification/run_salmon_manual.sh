#!/bin/bash
export PATH=/home/angsuman/miniconda3/envs/downstream/bin:$PATH
cd /home/angsuman/rna_pipeline

echo "Manual Salmon pipeline started at $(date)" > salmon_manual.log

# 1. Indexing
echo "[1/4] Indexing transcriptome..." >> salmon_manual.log
salmon index   -t /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta   -i /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta.salmon.idx   -k 31 -p 12 >> salmon_manual.log 2>&1

# 2. Quantifying AZ1
echo "[2/4] Quantifying AZ1..." >> salmon_manual.log
salmon quant   -i /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta.salmon.idx   -l A   -1 /home/angsuman/rna_pipeline/trimmed/AZ1_1.trimmed.fq.gz   -2 /home/angsuman/rna_pipeline/trimmed/AZ1_2.trimmed.fq.gz   -p 12   -o /home/angsuman/rna_pipeline/trinity_output/salmon_out/AZ1 >> salmon_manual.log 2>&1

# 3. Quantifying AZ2
echo "[3/4] Quantifying AZ2..." >> salmon_manual.log
salmon quant   -i /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta.salmon.idx   -l A   -1 /home/angsuman/rna_pipeline/trimmed/AZ2_1.trimmed.fq.gz   -2 /home/angsuman/rna_pipeline/trimmed/AZ2_2.trimmed.fq.gz   -p 12   -o /home/angsuman/rna_pipeline/trinity_output/salmon_out/AZ2 >> salmon_manual.log 2>&1

# 4. Quantifying AZ3
echo "[4/4] Quantifying AZ3..." >> salmon_manual.log
salmon quant   -i /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta.salmon.idx   -l A   -1 /home/angsuman/rna_pipeline/trimmed/AZ3_1.trimmed.fq.gz   -2 /home/angsuman/rna_pipeline/trimmed/AZ3_2.trimmed.fq.gz   -p 12   -o /home/angsuman/rna_pipeline/trinity_output/salmon_out/AZ3 >> salmon_manual.log 2>&1

echo "Manual Salmon pipeline finished at $(date)" >> salmon_manual.log
