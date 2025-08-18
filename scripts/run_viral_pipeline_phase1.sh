#!/bin/bash
# Viral Genomics Pipeline - Phase 1: Initial Analysis
# Runs assembly, depth analysis, and diagnostics for QC review
# Author: Kathie Mihindu
# Date: 2025-08-18

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_status() {
    echo -e "${GREEN}[STATUS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# Usage function
usage() {
    cat << USAGE
========================================================================
VIRAL GENOMICS PIPELINE - PHASE 1: Initial Analysis & QC
========================================================================

This script runs the initial analysis modules for quality assessment:
  - Reference-guided assembly and SNP calling
  - Depth calculation and visualization
  - Diagnostic analysis (optional but recommended)

After Phase 1 completes, review the outputs to determine appropriate
quality, depth, and frequency cutoffs for Phase 2.

Usage: $0 -s SAMPLE_NAME -r READ1 -f READ2 -a ACCESSION -t THREADS [OPTIONS]

Required arguments:
    -s SAMPLE_NAME   Sample name (e.g., SMS_5_DENV1)
    -r READ1         Path to R1 fastq.gz file
    -f READ2         Path to R2 fastq.gz file  
    -a ACCESSION     Reference genome accession (e.g., NC_001477.1)
    -t THREADS       Number of threads to use

Optional arguments:
    -o OUTPUT_DIR    Output directory (default: current directory)
    -D               Skip diagnostic module
    -h               Show this help message

Example:
    $0 -s SMS_5_DENV1 \\
       -r ./sample_R1.fastq.gz \\
       -f ./sample_R2.fastq.gz \\
       -a NC_001477.1 \\
       -t 4
USAGE
    exit 1
}

# Default parameters
OUTPUT_DIR="."
SKIP_DIAGNOSTIC=0

# Parse command line arguments
while getopts "s:r:f:a:t:o:Dh" opt; do
    case $opt in
        s) SAMPLE_NAME="$OPTARG" ;;
        r) READ1="$OPTARG" ;;
        f) READ2="$OPTARG" ;;
        a) ACCESSION="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        D) SKIP_DIAGNOSTIC=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$ACCESSION" ] || [ -z "$THREADS" ]; then
    print_error "Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$READ1" ]; then
    print_error "Read 1 file not found: $READ1"
    exit 1
fi

if [ ! -f "$READ2" ]; then
    print_error "Read 2 file not found: $READ2"
    exit 1
fi

# Set up paths
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"

# Create output directory structure
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}_results"
RESULTS_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}_results"

echo ""
echo "======================================================================="
echo "         VIRAL GENOMICS PIPELINE - PHASE 1: Initial Analysis"
echo "======================================================================="
print_status "Sample: $SAMPLE_NAME"
print_status "Reference: $ACCESSION"
print_status "Output directory: $RESULTS_DIR"
echo ""

# Submit both jobs in parallel
print_status "Submitting parallel analysis jobs..."

# Module 1: Reference-guided assembly and SNP calling
print_info "Module 1: Submitting reference-guided assembly..."
ASSEMBLY_JOB=$(sbatch --parsable ${PIPELINE_BASE}/viral_pipeline/analysis/submit_viral_pipeline.sh \\
    "$READ1" "$READ2" "$ACCESSION" "$THREADS")
print_status "Assembly job submitted: $ASSEMBLY_JOB"

# Module 6: Diagnostic report (run in parallel)
if [ $SKIP_DIAGNOSTIC -eq 0 ]; then
    print_info "Module 6: Submitting diagnostic analysis..."
    DIAG_JOB=$(sbatch --parsable ${PIPELINE_BASE}/viral_pipeline/analysis/submit_viral_diagnostic.sh \\
        "$READ1" "$READ2" "$ACCESSION" "diagnostic_${SAMPLE_NAME}" "$THREADS")
    print_status "Diagnostic job submitted: $DIAG_JOB"
else
    print_info "Skipping diagnostic module (use -D flag to enable)"
fi

# Wait for assembly job to complete
print_status "Waiting for assembly to complete (Job: $ASSEMBLY_JOB)..."
while squeue -j $ASSEMBLY_JOB 2>/dev/null | grep -q $ASSEMBLY_JOB; do
    sleep 30
done

# Check if assembly job completed successfully
if ! sacct -j $ASSEMBLY_JOB --format=State --noheader | grep -q COMPLETED; then
    print_error "Assembly job failed. Check logs for job $ASSEMBLY_JOB"
    exit 1
fi

print_status "Assembly completed successfully"

# Module 2: Generate depth file
print_info "Module 2: Generating depth file..."
$MAMBA_CMD bash -c "
    echo -e 'chrom\tposition\tdepth' > ${RESULTS_DIR}/${SAMPLE_NAME}_depth.txt &&
    samtools depth cleaned_seqs/mapping/*.lofreq.final.bam >> ${RESULTS_DIR}/${SAMPLE_NAME}_depth.txt
"
print_status "Depth file generated"

# Module 5: Visualize depth (run early for QC)
print_info "Module 5: Generating depth visualization..."
$MAMBA_CMD python3 ${PIPELINE_BASE}/viral_pipeline/visualization/visualize_depth.py \\
    --depth ${RESULTS_DIR}/${SAMPLE_NAME}_depth.txt \\
    --output ${RESULTS_DIR}/${SAMPLE_NAME}_depth.png \\
    --output-html ${RESULTS_DIR}/${SAMPLE_NAME}_depth.html \\
    --accession $ACCESSION
print_status "Depth visualization completed"

# Generate initial SNP statistics for QC review
print_info "Generating SNP statistics for quality assessment..."
$MAMBA_CMD python3 -c "
import pandas as pd
import numpy as np

# Read the SNP file
snp_file = 'cleaned_seqs/variants/' + [f for f in os.listdir('cleaned_seqs/variants/') if f.endswith('.snpEFF.ann.tsv')][0]
df = pd.read_csv(snp_file, sep='\t')

# Calculate statistics
print('\n========== SNP QUALITY STATISTICS ==========')
print(f'Total variants: {len(df)}')
print(f'\nQuality score distribution:')
print(f'  Min: {df[\"QUAL\"].min():.1f}')
print(f'  25%%: {df[\"QUAL\"].quantile(0.25):.1f}')
print(f'  Median: {df[\"QUAL\"].quantile(0.5):.1f}')
print(f'  75%%: {df[\"QUAL\"].quantile(0.75):.1f}')
print(f'  Max: {df[\"QUAL\"].max():.1f}')
print(f'\nDepth distribution:')
print(f'  Min: {df[\"DP\"].min()}')
print(f'  25%%: {df[\"DP\"].quantile(0.25):.0f}')
print(f'  Median: {df[\"DP\"].quantile(0.5):.0f}')
print(f'  75%%: {df[\"DP\"].quantile(0.75):.0f}')
print(f'  Max: {df[\"DP\"].max()}')
print(f'\nFrequency distribution (AF):')
print(f'  Variants >= 1%%: {len(df[df[\"AF\"] >= 0.01])}')
print(f'  Variants >= 5%%: {len(df[df[\"AF\"] >= 0.05])}')
print(f'  Variants >= 10%%: {len(df[df[\"AF\"] >= 0.10])}')
print(f'  Variants >= 50%%: {len(df[df[\"AF\"] >= 0.50])}')
print('============================================\n')
" 2>/dev/null || echo "Could not generate SNP statistics"

# Wait for diagnostic job if it was submitted
if [ $SKIP_DIAGNOSTIC -eq 0 ]; then
    print_status "Waiting for diagnostic analysis to complete (Job: $DIAG_JOB)..."
    while squeue -j $DIAG_JOB 2>/dev/null | grep -q $DIAG_JOB; do
        sleep 30
    done
    if sacct -j $DIAG_JOB --format=State --noheader | grep -q COMPLETED; then
        print_status "Diagnostic analysis completed"
    else
        print_error "Diagnostic job failed - continuing with main pipeline"
    fi
fi

echo ""
echo "======================================================================="
echo "                    PHASE 1 COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
print_status "Initial analysis complete. Review the following outputs:"
echo ""
echo "  1. Depth visualization: ${RESULTS_DIR}/${SAMPLE_NAME}_depth.png"
echo "  2. Interactive depth: ${RESULTS_DIR}/${SAMPLE_NAME}_depth.html"
echo "  3. Raw variants: cleaned_seqs/variants/*.snpEFF.ann.tsv"
if [ $SKIP_DIAGNOSTIC -eq 0 ]; then
    echo "  4. Diagnostic report: diagnostic_${SAMPLE_NAME}/*"
fi
echo ""
print_info "Based on the QC metrics above, determine appropriate cutoffs for:"
echo "  - Quality score (QUAL) - suggested: median or higher"
echo "  - Minimum depth (DP) - suggested: based on depth plot"
echo "  - Minimum frequency (AF) - suggested: 0.01 for minor variants"
echo ""
print_status "Next step: Run Phase 2 with your chosen parameters"
echo ""
echo "Example command for Phase 2:"
echo "  ./run_viral_pipeline_phase2.sh -s $SAMPLE_NAME \\"
echo "     -a $ACCESSION \\"
echo "     -q 41645 \\"
echo "     -d 200 \\"
echo "     -m 0.01"
echo ""
echo "======================================================================="

# Save phase 1 completion marker
echo "$(date): Phase 1 completed for $SAMPLE_NAME" > ${RESULTS_DIR}/.phase1_complete

