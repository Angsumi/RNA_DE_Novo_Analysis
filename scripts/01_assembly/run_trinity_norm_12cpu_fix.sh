#!/bin/bash
export PATH=/home/angsuman/miniconda3/bin:$PATH
cd /home/angsuman/rna_pipeline

# Remove the broken symlinks just in case
rm -f trinity_output_old/insilico_read_normalization/left.norm.fq
rm -f trinity_output_old/insilico_read_normalization/right.norm.fq

# Point directly to the actual files
LEFT_FQ="/home/angsuman/rna_pipeline/trinity_output_old/insilico_read_normalization/AZ1_1.trimmed.fq.gz_ext_all_reads.normalized_K25_C50_pctSD10000.fq"
RIGHT_FQ="/home/angsuman/rna_pipeline/trinity_output_old/insilico_read_normalization/AZ1_2.trimmed.fq.gz_ext_all_reads.normalized_K25_C50_pctSD10000.fq"

nohup perl /home/angsuman/miniconda3/bin/Trinity --seqType fq \
  --left $LEFT_FQ \
  --right $RIGHT_FQ \
  --max_memory 120G \
  --CPU 12 \
  --no_normalize_reads \
  --output /home/angsuman/rna_pipeline/trinity_output > trinity_output_norm.log 2>&1 &
echo "Trinity started with fixed paths!"
