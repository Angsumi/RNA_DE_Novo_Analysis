#!/bin/bash
export PATH=/home/angsuman/miniconda3/envs/downstream/bin:/home/angsuman/miniconda3/bin:$PATH
cd /home/angsuman/rna_pipeline

# 1. Create samples.txt
cat << 'EOT' > samples.txt
AZ1    AZ1_rep1    /home/angsuman/rna_pipeline/trimmed/AZ1_1.trimmed.fq.gz    /home/angsuman/rna_pipeline/trimmed/AZ1_2.trimmed.fq.gz
AZ2    AZ2_rep1    /home/angsuman/rna_pipeline/trimmed/AZ2_1.trimmed.fq.gz    /home/angsuman/rna_pipeline/trimmed/AZ2_2.trimmed.fq.gz
AZ3    AZ3_rep1    /home/angsuman/rna_pipeline/trimmed/AZ3_1.trimmed.fq.gz    /home/angsuman/rna_pipeline/trimmed/AZ3_2.trimmed.fq.gz
EOT

echo "Downstream analyses started at $(date)" > downstream_run.log

# 2. Run Salmon abundance estimation for all samples
echo "Starting Salmon Abundance Estimation..." >> downstream_run.log
perl /home/angsuman/miniconda3/opt/trinity-2.5.1/util/align_and_estimate_abundance.pl   --transcripts /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta   --seqType fq   --samples_file /home/angsuman/rna_pipeline/samples.txt   --est_method salmon   --prep_reference   --trinity_mode   --thread_count 12   --output_dir /home/angsuman/rna_pipeline/trinity_output/salmon_out >> downstream_run.log 2>&1

echo "Salmon Abundance Estimation completed at $(date)" >> downstream_run.log

# 3. Run BUSCO completeness
echo "Starting BUSCO completeness analysis..." >> downstream_run.log
cd /home/angsuman/rna_pipeline/trinity_output
/home/angsuman/miniconda3/envs/downstream/bin/python /home/angsuman/miniconda3/envs/downstream/bin/busco   -i /home/angsuman/rna_pipeline/trinity_output/Trinity.fasta   -m transcriptome   --auto-lineage-euk   -o busco_results   -c 12 >> /home/angsuman/rna_pipeline/downstream_run.log 2>&1

echo "BUSCO completeness analysis completed at $(date)" >> /home/angsuman/rna_pipeline/downstream_run.log
echo "All downstream analyses finished successfully!" >> /home/angsuman/rna_pipeline/downstream_run.log
