#!/bin/bash

# Script: preprocess_and_run_blaze.sh
# Usage : ./preprocess_and_run_blaze.sh <file_path> <file_type: bam|fastq> <path_to_ref_whitelist> <sample_name>

set -euo pipefail  # Safer error handling: exit on error, undefined variable, or failed pipe

# === Input arguments ===
file_path=$1        # Path to directory containing .bam or .fastq files
file_type=$2        # Either "bam" or "fastq"
path_to_ref=$3      # Path to barcode whitelist (e.g. 737K-august-2016.txt)
sample_name=$4      # Sample name for output file prefixing

# === Input validation ===
if [[ ! -d "$file_path" ]]; then
    echo "Error: Directory '$file_path' does not exist." >&2
    exit 1
fi

if [[ "$file_type" != "bam" && "$file_type" != "fastq" ]]; then
    echo "Error: file_type must be either 'bam' or 'fastq'." >&2
    exit 1
fi

if [[ ! -f "$path_to_ref" ]]; then
    echo "Error: Whitelist file not found at $path_to_ref." >&2
    exit 1
fi

# === Create necessary directories ===
mkdir -p "$file_path/fastq"
mkdir -p "$file_path/parts"
mkdir -p "$file_path/blaze/merged"

# === Step 1: Convert BAM to FASTQ if needed ===
if [[ "$file_type" == "bam" ]]; then
    echo "Converting BAM files to FASTQ..."
    for bam in "$file_path"/*.bam; do
        base=$(basename "$bam" .bam)
        samtools fastq "$bam" > "$file_path/fastq/${base}.fastq"
    done
fi

# === Step 2: Merge FASTQ files by barcode groups ===
echo "Merging FASTQ files by barcode groups..."
cat "$file_path"/fastq/*_1*.fastq > "$file_path/parts/merged1.fastq" || true
cat "$file_path"/fastq/*_2*.fastq > "$file_path/parts/merged2.fastq" || true
cat "$file_path"/fastq/*_3*.fastq > "$file_path/parts/merged3.fastq" || true
cat "$file_path"/fastq/*_[4567890]*.fastq > "$file_path/parts/merged4.fastq" || true

# === Step 3: Activate conda environment for Blaze ===
echo "Activating Conda environment 'blaze_env'..."
source ~/.bashrc
conda activate blaze_env

# === Step 4: Run Blaze on each merged FASTQ group ===
echo "Running Blaze for each merged group..."
for i in 1 2 3 4; do
    input_fastq="$file_path/parts/merged${i}.fastq"
    if [[ -f "$input_fastq" ]]; then
        blaze --expect-cells=2000 --threads=20 \
            --full-bc-whitelist="$path_to_ref" \
            --output-prefix "${sample_name}_part${i}" \
            "$input_fastq"
    fi
done

# === Step 5: Final merge of Blaze output ===
echo "Merging all matched reads from Blaze..."
cat "$file_path"/blaze/*.fastq.gz > "$file_path/blaze/merged/merged_matched_reads_${sample_name}.fastq.gz"

echo "Process complete. Merged output saved to:"
echo "$file_path/blaze/merged/merged_matched_reads_${sample_name}.fastq.gz"

