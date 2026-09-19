#!/bin/bash
export PATH=/home/angsuman/miniconda3/envs/downstream/bin:/home/angsuman/miniconda3/bin:$PATH
cd /home/angsuman/rna_pipeline/trinity_output

echo "Trinotate pipeline started at $(date)" > trinotate_pipeline.log

# 1. Download and create databases
echo "[1/5] Downloading and creating Trinotate database..." >> trinotate_pipeline.log
mkdir -p /home/angsuman/trinotate_data
Trinotate --db /home/angsuman/rna_pipeline/trinity_output/trinotate.sqlite --create --trinotate_data_dir /home/angsuman/trinotate_data >> trinotate_pipeline.log 2>&1

# 2. Run TransDecoder
echo "[2/5] Running TransDecoder LongOrfs..." >> trinotate_pipeline.log
TransDecoder.LongOrfs -t Trinity.fasta >> trinotate_pipeline.log 2>&1

echo "[2/5] Running TransDecoder Predict..." >> trinotate_pipeline.log
TransDecoder.Predict -t Trinity.fasta >> trinotate_pipeline.log 2>&1

# 3. Initialize database
echo "[3/5] Initializing Trinotate database..." >> trinotate_pipeline.log
Trinotate --db trinotate.sqlite --init   --gene_trans_map Trinity.fasta.gene_trans_map   --transcript_fasta Trinity.fasta   --transdecoder_pep Trinity.fasta.transdecoder.pep >> trinotate_pipeline.log 2>&1

# 4. Run computes with DIAMOND
echo "[4/5] Running annotation computes (using DIAMOND)..." >> trinotate_pipeline.log
Trinotate --db trinotate.sqlite --run ALL --CPU 12 --use_diamond   --transcript_fasta Trinity.fasta   --transdecoder_pep Trinity.fasta.transdecoder.pep   --trinotate_data_dir /home/angsuman/trinotate_data >> trinotate_pipeline.log 2>&1

# 5. Report generation
echo "[5/5] Generating functional annotation report..." >> trinotate_pipeline.log
Trinotate --db trinotate.sqlite --report > trinotate_annotation_report.xls 2>> trinotate_pipeline.log

echo "Trinotate pipeline finished successfully at $(date)" >> trinotate_pipeline.log
