#!/bin/bash
export PATH=/home/angsuman/miniconda3/bin:$PATH
cd /home/angsuman/rna_pipeline
nohup perl /home/angsuman/miniconda3/bin/Trinity --seqType fq \
  --left /home/angsuman/rna_pipeline/trinity_output_old/insilico_read_normalization/left.norm.fq \
  --right /home/angsuman/rna_pipeline/trinity_output_old/insilico_read_normalization/right.norm.fq \
  --max_memory 90G \
  --CPU 12 \
  --no_normalize_reads \
  --output /home/angsuman/rna_pipeline/trinity_output > trinity_output_norm.log 2>&1 &
echo "Trinity started!"
