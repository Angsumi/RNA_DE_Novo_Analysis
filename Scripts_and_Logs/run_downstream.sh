#!/bin/bash
export PATH=/home/angsuman/miniconda3/bin:$PATH

# PIDs and Paths
TRINITY_PID=$(pgrep -f "perl /home/angsuman/miniconda3/bin/Trinity" | head -n 1)
TRINITY_FASTA="/home/angsuman/rna_pipeline/trinity_output/Trinity.fasta"
OUT_STATS="/home/angsuman/rna_pipeline/trinity_output/TrinityStats.txt"
LOG_FILE="/home/angsuman/rna_pipeline/downstream_analysis.log"
TRINITY_HOME="/home/angsuman/miniconda3/opt/trinity-2.5.1"

# The read files used for Trinity. For best abundance estimation, 
# raw trimmed reads should be used, but we use the normalized ones here 
# as they are guaranteed to exist and match what was assembled.
LEFT_FQ="/home/angsuman/rna_pipeline/trinity_output_old/insilico_read_normalization/AZ1_1.trimmed.fq.gz_ext_all_reads.normalized_K25_C50_pctSD10000.fq"
RIGHT_FQ="/home/angsuman/rna_pipeline/trinity_output_old/insilico_read_normalization/AZ1_2.trimmed.fq.gz_ext_all_reads.normalized_K25_C50_pctSD10000.fq"

echo "Starting downstream monitoring at $(date)" > $LOG_FILE
echo "Waiting for Trinity (PID $TRINITY_PID) to finish..." >> $LOG_FILE

while kill -0 $TRINITY_PID 2>/dev/null; do
    sleep 60
done

echo "Trinity process finished at $(date). Checking for output fasta..." >> $LOG_FILE

if [ -f "$TRINITY_FASTA" ]; then
    echo "Found $TRINITY_FASTA." >> $LOG_FILE
    
    # 1. Basic Stats
    echo "Running TrinityStats.pl..." >> $LOG_FILE
    perl $TRINITY_HOME/util/TrinityStats.pl $TRINITY_FASTA > $OUT_STATS
    echo "TrinityStats finished. Results saved to $OUT_STATS." >> $LOG_FILE

    # 2. Abundance Estimation (Salmon)
    echo "Starting Abundance Estimation with Salmon..." >> $LOG_FILE
    cd /home/angsuman/rna_pipeline/trinity_output
    if command -v salmon &> /dev/null; then
        perl $TRINITY_HOME/util/align_and_estimate_abundance.pl \
            --transcripts $TRINITY_FASTA \
            --seqType fq \
            --left $LEFT_FQ \
            --right $RIGHT_FQ \
            --est_method salmon \
            --prep_reference \
            --thread_count 12 \
            --output_dir salmon_outdir >> $LOG_FILE 2>&1
        echo "Abundance estimation finished." >> $LOG_FILE
    else
        echo "Salmon is not in PATH. Skipping abundance estimation." >> $LOG_FILE
    fi

    # 3. BUSCO completeness
    echo "Starting BUSCO completeness check..." >> $LOG_FILE
    if command -v busco &> /dev/null; then
        # Running busco in auto-lineage mode (eukaryota)
        busco -i $TRINITY_FASTA -m transcriptome --auto-lineage-euk -o busco_results -c 12 >> $LOG_FILE 2>&1
        echo "BUSCO finished." >> $LOG_FILE
    else
        echo "BUSCO is not in PATH. Skipping BUSCO." >> $LOG_FILE
    fi

    echo "All downstream analyses completed at $(date)." >> $LOG_FILE

else
    echo "Error: $TRINITY_FASTA not found. The assembly might have encountered an error during the final steps." >> $LOG_FILE
fi
