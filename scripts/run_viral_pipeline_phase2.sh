#!/bin/bash
# Viral Genomics Pipeline - Phase 2: Filtering and Analysis
# Runs mutation filtering, visualization, and haplotype reconstruction
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

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Usage function
usage() {
    cat << USAGE
========================================================================
VIRAL GENOMICS PIPELINE - PHASE 2: Filtering & Analysis
========================================================================

This script runs the filtering and analysis modules using QC-informed
parameters from Phase 1:
  - Parse and filter mutations based on quality cutoffs
  - Generate mutation visualizations
  - Perform haplotype reconstruction and quasispecies analysis

Prerequisites: Phase 1 must be completed first

Usage: $0 -s SAMPLE_NAME -a ACCESSION -q QUALITY -d DEPTH -m MIN_FREQ [OPTIONS]

Required arguments:
    -s SAMPLE_NAME   Sample name (must match Phase 1)
    -a ACCESSION     Reference genome accession (e.g., NC_001477.1)
    -q QUALITY       Quality threshold for filtering (e.g., 41645)
    -d DEPTH         Minimum depth for filtering (e.g., 200)
    -m MIN_FREQ      Minimum frequency for mutations (e.g., 0.01)

Optional arguments:
    -o OUTPUT_DIR    Output directory (default: current directory)
    -h               Show this help message

Example:
    $0 -s SMS_5_DENV1 \\
       -a NC_001477.1 \\
       -q 41645 \\
       -d 200 \\
       -m 0.01
USAGE
    exit 1
}

# Default parameters
OUTPUT_DIR="."

# Parse command line arguments
while getopts "s:a:q:d:m:o:h" opt; do
    case $opt in
        s) SAMPLE_NAME="$OPTARG" ;;
        a) ACCESSION="$OPTARG" ;;
        q) QUALITY="$OPTARG" ;;
        d) DEPTH="$OPTARG" ;;
        m) MIN_FREQ="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$ACCESSION" ] || [ -z "$QUALITY" ] || [ -z "$DEPTH" ] || [ -z "$MIN_FREQ" ]; then
    print_error "Missing required arguments"
    usage
fi

# Set up paths
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"
RESULTS_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}_results"

# Check if Phase 1 was completed
if [ ! -f "${RESULTS_DIR}/.phase1_complete" ]; then
    print_error "Phase 1 has not been completed for this sample"
    print_info "Please run run_viral_pipeline_phase1.sh first"
    exit 1
fi

# Check for required input files from Phase 1
if [ ! -f "${RESULTS_DIR}/${SAMPLE_NAME}_depth.txt" ]; then
    print_error "Depth file not found. Please ensure Phase 1 completed successfully"
    exit 1
fi

echo ""
echo "======================================================================="
echo "         VIRAL GENOMICS PIPELINE - PHASE 2: Filtering & Analysis"
echo "======================================================================="
print_status "Sample: $SAMPLE_NAME"
print_status "Reference: $ACCESSION"
print_status "Parameters:"
echo "  - Quality threshold: $QUALITY"
echo "  - Minimum depth: $DEPTH"
echo "  - Minimum frequency: $MIN_FREQ"
echo ""

# Module 3: Parse and filter mutations
print_info "Module 3: Parsing and filtering mutations..."

# Find the SNP file
SNP_FILE=$(ls cleaned_seqs/variants/*.snpEFF.ann.tsv 2>/dev/null | head -1)
if [ -z "$SNP_FILE" ]; then
    print_error "SNP file not found in cleaned_seqs/variants/"
    exit 1
fi

$MAMBA_CMD python3 ${PIPELINE_BASE}/viral_pipeline/visualization/parse_snpeff_tsv.py \\
    "$SNP_FILE" \\
    ${RESULTS_DIR}/${SAMPLE_NAME}_filtered_mutations.tsv \\
    --quality $QUALITY \\
    --depth $DEPTH \\
    --freq $MIN_FREQ

# Count filtered mutations
MUTATION_COUNT=$(tail -n +2 ${RESULTS_DIR}/${SAMPLE_NAME}_filtered_mutations.tsv | wc -l)
print_status "Filtered mutations: $MUTATION_COUNT variants pass filters"

# Module 4: Visualize mutations
print_info "Module 4: Generating mutation visualization..."
$MAMBA_CMD python3 ${PIPELINE_BASE}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \\
    --input ${RESULTS_DIR}/${SAMPLE_NAME}_filtered_mutations.tsv \\
    --output ${SAMPLE_NAME}_mutations.png \\
    --accession $ACCESSION \\
    --cutoff $MIN_FREQ \\
    --outdir ${RESULTS_DIR}

print_status "Mutation visualization completed"

# Module 7: Generate haplotype consensus
print_info "Module 7: Generating haplotype consensus and proteins..."

# Find reference file
REF_FILE="cleaned_seqs/${ACCESSION}.fasta"
if [ ! -f "$REF_FILE" ]; then
    print_warning "Reference file not found at $REF_FILE, searching..."
    REF_FILE=$(find cleaned_seqs -name "*.fasta" -o -name "*.fa" | head -1)
    if [ -z "$REF_FILE" ]; then
        print_error "No reference file found"
        exit 1
    fi
fi

$MAMBA_CMD python3 ${PIPELINE_BASE}/viral_pipeline/utils/generate_realistic_haplotype_consensus.py \\
    --vcf ${RESULTS_DIR}/${SAMPLE_NAME}_filtered_mutations.tsv \\
    --reference "$REF_FILE" \\
    --accession $ACCESSION \\
    --quality $QUALITY \\
    --depth $DEPTH \\
    --freq $MIN_FREQ \\
    --output-prefix ${RESULTS_DIR}/${SAMPLE_NAME}

print_status "Haplotype reconstruction completed"

# Generate summary report
print_info "Generating analysis summary..."
cat > ${RESULTS_DIR}/${SAMPLE_NAME}_phase2_summary.txt << SUMMARY
========================================
PHASE 2 ANALYSIS SUMMARY
========================================
Sample: $SAMPLE_NAME
Reference: $ACCESSION
Date: $(date)

FILTERING PARAMETERS:
- Quality threshold: $QUALITY
- Minimum depth: $DEPTH  
- Minimum frequency: $MIN_FREQ

RESULTS:
- Total filtered mutations: $MUTATION_COUNT
- Output files generated:
  * Filtered mutations: ${SAMPLE_NAME}_filtered_mutations.tsv
  * Mutation visualization: ${SAMPLE_NAME}_mutations.png
  * Mutation table: ${SAMPLE_NAME}_mutations_mutations_table.tsv
  * Haplotype sequences: ${SAMPLE_NAME}_realistic_haplotypes.fasta
  * Protein sequences: ${SAMPLE_NAME}_realistic_proteins.fasta
  * Haplotype report: ${SAMPLE_NAME}_realistic_haplotype_report.txt
========================================
SUMMARY

echo ""
echo "======================================================================="
echo "                    PHASE 2 COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
print_status "Analysis complete! Results available in: $RESULTS_DIR"
echo ""
print_info "Key outputs:"
echo "  1. Filtered mutations: ${SAMPLE_NAME}_filtered_mutations.tsv"
echo "  2. Mutation visualization: ${SAMPLE_NAME}_mutations.png"
echo "  3. Haplotype consensus: ${SAMPLE_NAME}_realistic_haplotypes.fasta"
echo "  4. Translated proteins: ${SAMPLE_NAME}_realistic_proteins.fasta"
echo "  5. Analysis summary: ${SAMPLE_NAME}_phase2_summary.txt"
echo ""

# Quality check suggestions
if [ $MUTATION_COUNT -eq 0 ]; then
    print_warning "No mutations passed the filters. Consider:"
    echo "  - Lowering the quality threshold"
    echo "  - Reducing the minimum depth requirement"
    echo "  - Checking the depth plot for coverage issues"
elif [ $MUTATION_COUNT -gt 100 ]; then
    print_warning "High number of mutations ($MUTATION_COUNT). Consider:"
    echo "  - Increasing the quality threshold"
    echo "  - Raising the minimum frequency cutoff"
    echo "  - Reviewing for potential contamination"
fi

echo ""
echo "======================================================================="

# Mark phase 2 completion
echo "$(date): Phase 2 completed with Q=$QUALITY, D=$DEPTH, F=$MIN_FREQ" >> ${RESULTS_DIR}/.phase2_complete

