#!/bin/bash
# Script to fix polyprotein-only annotations using BLAST
# Integrates with viral genomics pipeline

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <accession> <refseq_accession>"
    echo "Example: $0 KU955591.1 NC_012532.1"
    echo ""
    echo "Common RefSeq viruses with good annotations:"
    echo "  Zika: NC_012532.1 or NC_035889.1"
    echo "  Dengue 1: NC_001477.1"
    echo "  West Nile: NC_009942.1"
    echo "  VEEV: NC_001449.1"
    exit 1
fi

ACCESSION=$1
REFSEQ=$2

# Set paths
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
SCRIPT_PATH="${PIPELINE_BASE}/viral_pipeline/utils/annotate_virus_by_blast.py"

echo "========================================="
echo "Fixing annotation for $ACCESSION"
echo "Using RefSeq: $REFSEQ"
echo "========================================="

# Find genome file
if [ -f "cleaned_seqs/${ACCESSION}.fasta" ]; then
    GENOME_FILE="cleaned_seqs/${ACCESSION}.fasta"
elif [ -f "${ACCESSION}.fasta" ]; then
    GENOME_FILE="${ACCESSION}.fasta"
else
    echo "❌ Error: Cannot find genome file ${ACCESSION}.fasta"
    echo "Please run Module 1 first to generate the reference"
    exit 1
fi

echo "Using genome file: $GENOME_FILE"

# Create output directory
OUTPUT_DIR="blast_annotation_${ACCESSION}"
mkdir -p "$OUTPUT_DIR"

# Run BLAST annotation
echo ""
echo "Running BLAST annotation..."
$MAMBA_CMD python3 "$SCRIPT_PATH" \
    "$GENOME_FILE" \
    --refseq "$REFSEQ" \
    --output-dir "$OUTPUT_DIR" \
    --update-snpeff \
    --min-coverage 0.8

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Annotation successful!"
    echo ""
    echo "Next steps:"
    echo "1. Review gene coordinates in: $OUTPUT_DIR/${ACCESSION}_annotation_summary.json"
    echo "2. Update SnpEff database with: $OUTPUT_DIR/${ACCESSION}.gff3"
    echo "3. Re-run Module 1 with updated annotation"
    echo "4. Run Module 8 for proper protein translation"
else
    echo ""
    echo "❌ Annotation failed!"
    exit 1
fi