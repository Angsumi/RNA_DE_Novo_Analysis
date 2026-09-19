#!/bin/bash
export PATH=/home/angsuman/miniconda3/envs/downstream/bin:/home/angsuman/miniconda3/bin:$PATH
cd /home/angsuman/rna_pipeline/trinity_output

# Create sample relationships
echo -e "AZ1\tAZ1\nAZ2\tAZ2\nAZ3\tAZ3" > de_samples.txt

echo "DE Analysis started at $(date)" > de_analysis.log

# 1. Gene-level DE
echo "Running gene-level DE analysis..." >> de_analysis.log
perl /home/angsuman/miniconda3/opt/trinity-2.5.1/Analysis/DifferentialExpression/run_DE_analysis.pl   --matrix /home/angsuman/rna_pipeline/trinity_output/trinity_matrix.gene.counts.matrix   --method edgeR   --samples_file /home/angsuman/rna_pipeline/trinity_output/de_samples.txt   --dispersion 0.1   --output /home/angsuman/rna_pipeline/trinity_output/edgeR_gene_dir >> de_analysis.log 2>&1

# 2. Isoform-level DE
echo "Running isoform-level DE analysis..." >> de_analysis.log
perl /home/angsuman/miniconda3/opt/trinity-2.5.1/Analysis/DifferentialExpression/run_DE_analysis.pl   --matrix /home/angsuman/rna_pipeline/trinity_output/trinity_matrix.isoform.counts.matrix   --method edgeR   --samples_file /home/angsuman/rna_pipeline/trinity_output/de_samples.txt   --dispersion 0.1   --output /home/angsuman/rna_pipeline/trinity_output/edgeR_isoform_dir >> de_analysis.log 2>&1

echo "DE Analysis completed successfully at $(date)" >> de_analysis.log
