#!/bin/bash
#SBATCH --job-name=genbank_parse
#SBATCH --output=genbank_parse_%j.out
#SBATCH --error=genbank_parse_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Script to parse mat_peptide annotations from GenBank file
# Usage: sbatch submit_genbank_annotation.sh ACCESSION

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: sbatch submit_genbank_annotation.sh ACCESSION"
    echo "Example: sbatch submit_genbank_annotation.sh NC_001477.1"
    exit 1
fi

ACCESSION=$1

# Change to submission directory
cd ${SLURM_SUBMIT_DIR}

echo "GenBank mat_peptide Parser"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "=========================================="
echo "Target: $ACCESSION"
echo "=========================================="

# Set paths
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"
SCRIPT_DIR="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/analysis"
OUTPUT_DIR="blast_annotation_${ACCESSION}"
GENBANK_FILE="${ACCESSION}.gb"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Download GenBank file if not present
if [ \! -f "$GENBANK_FILE" ]; then
    echo "Downloading GenBank file..."
    $MAMBA_CMD python3 -c "
from Bio import Entrez
Entrez.email = 'your.email@example.com'
handle = Entrez.efetch(db='nucleotide', id='$ACCESSION', rettype='gb', retmode='text')
with open('$GENBANK_FILE', 'w') as f:
    f.write(handle.read())
handle.close()
print('Downloaded $GENBANK_FILE')
"
fi

# Parse GenBank annotations
echo ""
echo "Parsing mat_peptide annotations..."
$MAMBA_CMD python3 "${SCRIPT_DIR}/parse_genbank_annotations.py" "$GENBANK_FILE" --output-dir .

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✅ Annotation successful\!"
    echo ""
    echo "Results:"
    echo "  - Gene coordinates: $OUTPUT_DIR/${ACCESSION}_gene_coords.json"
    echo "  - GFF3 for SnpEff: $OUTPUT_DIR/${ACCESSION}_genes.gff"
    echo "  - Summary: $OUTPUT_DIR/${ACCESSION}_annotation_summary.json"
    echo ""
    echo "Next steps:"
    echo "  1. Review the annotation summary"
    echo "  2. Update SnpEff database if needed"
else
    echo "❌ Error: Annotation failed"
    exit 1
fi

echo ""
echo "Job completed at: $(date)"
