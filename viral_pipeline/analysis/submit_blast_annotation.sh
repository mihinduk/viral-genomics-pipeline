#!/bin/bash
#SBATCH --job-name=blast_annotation
#SBATCH --output=blast_annotation_%j.out
#SBATCH --error=blast_annotation_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=general

# SLURM submission script for BLAST-based gene annotation
# Usage: sbatch submit_blast_annotation.sh <accession> <refseq_accession>

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

# Change to submission directory
cd ${SLURM_SUBMIT_DIR}

echo "========================================="
echo "BLAST-based Gene Annotation"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "========================================="
echo "Target: $ACCESSION"
echo "Reference: $REFSEQ"
echo "========================================="

# Set paths
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
SCRIPT_PATH="${PIPELINE_BASE}/viral_pipeline/utils/annotate_virus_by_blast.py"

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

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed at: $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Annotation successful!"
    echo ""
    echo "Results:"
    echo "  - Gene coordinates: $OUTPUT_DIR/${ACCESSION}_annotation_summary.json"
    echo "  - GFF3 for SnpEff: $OUTPUT_DIR/${ACCESSION}.gff3"
    echo "  - known_viruses.json has been updated"
    echo ""
    echo "Next steps:"
    echo "  1. Review the annotation summary"
    echo "  2. Re-run variant calling if needed"
    echo "  3. Run Module 8 for proper protein translation"
else
    echo "❌ Annotation failed with exit code: $EXIT_CODE"
fi
echo "========================================="

exit $EXIT_CODE