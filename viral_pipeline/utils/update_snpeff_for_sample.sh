#!/bin/bash
# Wrapper script to update SnpEff database for a sample
# Usage: bash update_snpeff_for_sample.sh KU955591.1

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <accession>"
    echo "Example: $0 KU955591.1"
    echo ""
    echo "This script will:"
    echo "  1. Find the BLAST annotation results"
    echo "  2. Update the SnpEff database"
    echo "  3. Test the new database"
    exit 1
fi

ACCESSION=$1
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"

echo "========================================="
echo "Updating SnpEff Database"
echo "========================================="
echo "Accession: $ACCESSION"
echo "Working directory: $(pwd)"
echo ""

# Find required files
GENOME_FILE="cleaned_seqs/${ACCESSION}.fasta"
BLAST_DIR="blast_annotation_${ACCESSION}"
GFF_FILE="${BLAST_DIR}/${ACCESSION}.gff3"

# Check if files exist
if [ ! -f "$GENOME_FILE" ]; then
    echo "❌ Error: Genome file not found: $GENOME_FILE"
    echo "Please run Module 1 first to generate the reference"
    exit 1
fi

if [ ! -f "$GFF_FILE" ]; then
    echo "❌ Error: GFF file not found: $GFF_FILE"
    echo "Please run BLAST annotation first:"
    echo "  sbatch ${PIPELINE_BASE}/viral_pipeline/analysis/submit_blast_annotation.sh $ACCESSION <refseq_accession>"
    exit 1
fi

echo "Found required files:"
echo "  Genome: $GENOME_FILE"
echo "  GFF: $GFF_FILE"
echo ""

# First regenerate the GFF with simplified structure
echo "Regenerating simplified GFF..."
$MAMBA_CMD python3 "${PIPELINE_BASE}/viral_pipeline/utils/annotate_virus_by_blast.py" \
    "$GENOME_FILE" \
    --refseq NC_012532.1 \
    --output-dir "$BLAST_DIR" \
    --update-snpeff \
    --min-coverage 0.8

if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to regenerate GFF"
    exit 1
fi

# Run the simplified database updater
echo "Updating SnpEff database..."
$MAMBA_CMD python3 "${PIPELINE_BASE}/viral_pipeline/utils/update_snpeff_database_simple.py" \
    --accession "$ACCESSION" \
    --genome "$GENOME_FILE" \
    --gff "$GFF_FILE" \
    --backup \
    --test

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ SnpEff database updated successfully!"
    echo ""
    echo "Next steps:"
    echo "  1. Re-run variant calling to get proper annotations"
    echo "  2. Check that Gene_106_10377 is replaced with real protein names"
    echo "  3. Run Module 8 for consensus generation with individual proteins"
else
    echo ""
    echo "❌ SnpEff database update failed!"
    exit 1
fi