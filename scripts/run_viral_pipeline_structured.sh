#!/bin/bash
# Structured Viral Genomics Pipeline - Phase 1
# Uses standardized directory structure with predictable relative paths
# Author: Kathie Mihindu

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_status() { echo -e "${GREEN}[✓]${NC} $1"; }
print_error() { echo -e "${RED}[✗]${NC} $1" >&2; }
print_info() { echo -e "${BLUE}[i]${NC} $1"; }
print_step() { echo -e "\n${YELLOW}Step $1:${NC} $2"; }

usage() {
    cat << USAGE
========================================================================
STRUCTURED VIRAL GENOMICS PIPELINE
========================================================================

Creates standardized directory structure for seamless module integration.
All modules use predictable relative paths - no hardcoding required.

Usage: $0 -s SAMPLE_NAME -r READ1 -f READ2 -a ACCESSION -o OUTPUT_DIR -t THREADS

Required:
    -s SAMPLE_NAME   Sample identifier
    -r READ1         Path to R1 fastq.gz file
    -f READ2         Path to R2 fastq.gz file  
    -a ACCESSION     Reference genome accession
    -o OUTPUT_DIR    Output directory (will be created)
    -t THREADS       Number of threads

Example:
    $0 -s SMS_4_WNV \\
       -r /data/sample_R1.fastq.gz \\
       -f /data/sample_R2.fastq.gz \\
       -a AY532665.1 \\
       -o /scratch/results/SMS_4_analysis \\
       -t 4

Directory Structure Created:
    OUTPUT_DIR/
    ├── inputs/      # Raw inputs (linked)
    ├── qc/          # Quality control
    ├── mapping/     # Read alignment
    ├── variants/    # SNP calling
    ├── assembly/    # Diagnostic assembly
    ├── analysis/    # Phase 2 outputs
    ├── logs/        # All log files
    └── metadata/    # Pipeline status
USAGE
    exit 1
}

# Parse arguments
while getopts "s:r:f:a:o:t:h" opt; do
    case $opt in
        s) SAMPLE_NAME="$OPTARG" ;;
        r) READ1="$OPTARG" ;;
        f) READ2="$OPTARG" ;;
        a) ACCESSION="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$ACCESSION" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$THREADS" ]; then
    print_error "Missing required arguments"
    usage
fi

# Validate input files
for file in "$READ1" "$READ2"; do
    if [ ! -f "$file" ]; then
        print_error "Input file not found: $file"
        exit 1
    fi
done

# Pipeline configuration
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"

# Header
echo ""
echo "======================================================================="
echo "         STRUCTURED VIRAL GENOMICS PIPELINE - PHASE 1"
echo "======================================================================="
print_info "Sample: $SAMPLE_NAME"
print_info "Reference: $ACCESSION"
print_info "Output: $OUTPUT_DIR"
print_info "Threads: $THREADS"
echo ""

# Step 1: Database validation
print_step 1 "Validating virus databases"
${PIPELINE_BASE}/scripts/validate_virus_databases.sh "$ACCESSION"
if [ $? -ne 0 ]; then
    print_error "Database validation failed"
    exit 1
fi
print_status "Database validation passed"

# Step 2: Create directory structure
print_step 2 "Creating standardized directory structure"
mkdir -p "$OUTPUT_DIR"/{inputs,qc,mapping,variants,assembly,analysis,logs,metadata}

# Save pipeline metadata
cat > "$OUTPUT_DIR/metadata/sample_info.json" << METADATA
{
    "sample_name": "$SAMPLE_NAME",
    "accession": "$ACCESSION",
    "read1_original": "$READ1",
    "read2_original": "$READ2",
    "threads": $THREADS,
    "pipeline_version": "structured_v1.0",
    "created": "$(date -Iseconds)",
    "pipeline_base": "$PIPELINE_BASE"
}
METADATA

print_status "Directory structure created"

# Step 3: Setup inputs
print_step 3 "Setting up input files"
cd "$OUTPUT_DIR"

# Link input files with consistent naming
ln -sf "$READ1" "inputs/${SAMPLE_NAME}_R1.fastq.gz"
ln -sf "$READ2" "inputs/${SAMPLE_NAME}_R2.fastq.gz"

# Download reference sequence
$MAMBA_CMD python3 -c "
from Bio import Entrez, SeqIO
Entrez.email = 'pipeline@wustl.edu'
try:
    handle = Entrez.efetch(db='nucleotide', id='$ACCESSION', rettype='fasta', retmode='text')
    with open('inputs/reference.fasta', 'w') as f:
        f.write(handle.read())
    handle.close()
    print('Reference sequence downloaded')
except Exception as e:
    print(f'Error: {e}')
    exit(1)
"

print_status "Input files configured"

# Step 4: Submit assembly job
print_step 4 "Submitting structured assembly job"

# Create wrapper script that uses structured paths
cat > "scripts/structured_assembly_wrapper.sh" << WRAPPER
#!/bin/bash
#SBATCH --job-name=structured_assembly
#SBATCH --output=logs/assembly_%j.out
#SBATCH --error=logs/assembly_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=$THREADS

# Change to analysis directory
cd $OUTPUT_DIR

# Activate environment
source /ref/sahlab/software/anaconda3/bin/activate
$MAMBA_CMD bash -c "
    # Quality control
    fastp -i inputs/${SAMPLE_NAME}_R1.fastq.gz \
          -I inputs/${SAMPLE_NAME}_R2.fastq.gz \
          -o qc/${SAMPLE_NAME}_R1.qc.fastq.gz \
          -O qc/${SAMPLE_NAME}_R2.qc.fastq.gz \
          -h qc/fastp_report.html \
          -j qc/fastp_report.json \
          --thread $THREADS

    # Index reference
    bwa index inputs/reference.fasta

    # Map reads
    bwa mem -t $THREADS inputs/reference.fasta \
            qc/${SAMPLE_NAME}_R1.qc.fastq.gz \
            qc/${SAMPLE_NAME}_R2.qc.fastq.gz | \
    samtools view -bS - | \
    samtools sort -@ $THREADS -o mapping/${SAMPLE_NAME}.bam

    # Index BAM
    samtools index mapping/${SAMPLE_NAME}.bam

    # Variant calling with LoFreq
    lofreq call-parallel --pp-threads $THREADS \
                         -f inputs/reference.fasta \
                         -o variants/${SAMPLE_NAME}.vcf \
                         mapping/${SAMPLE_NAME}.bam

    # Annotate with SnpEff
    snpEff $ACCESSION variants/${SAMPLE_NAME}.vcf > variants/${SAMPLE_NAME}.snpEFF.ann.vcf

    # Convert to TSV
    java -jar \$SNPEFF_JAR/SnpSift.jar extractFields \
         variants/${SAMPLE_NAME}.snpEFF.ann.vcf \
         CHROM POS REF ALT QUAL DP AF \
         'ANN[0].EFFECT' 'ANN[0].IMPACT' 'ANN[0].GENE' 'ANN[0].FEATURE' \
         > variants/${SAMPLE_NAME}.snpEFF.ann.tsv
"

# Mark Phase 1 complete
echo "$(date -Iseconds): Phase 1 assembly completed" > metadata/.phase1_assembly_complete
WRAPPER

chmod +x "scripts/structured_assembly_wrapper.sh"

# Step 5: Submit diagnostic job
print_step 5 "Submitting structured diagnostic job"

cat > "scripts/structured_diagnostic_wrapper.sh" << DIAGNOSTIC
#!/bin/bash
#SBATCH --job-name=structured_diagnostic
#SBATCH --output=logs/diagnostic_%j.out
#SBATCH --error=logs/diagnostic_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=$THREADS

# Change to analysis directory
cd $OUTPUT_DIR

# Activate environment
source /ref/sahlab/software/anaconda3/bin/activate
$MAMBA_CMD bash -c "
    # Assembly with MEGAHIT
    megahit -1 qc/${SAMPLE_NAME}_R1.qc.fastq.gz \
            -2 qc/${SAMPLE_NAME}_R2.qc.fastq.gz \
            -o assembly/megahit_output \
            --presets meta-sensitive \
            --min-contig-len 500 \
            -t $THREADS

    # Copy final contigs
    if [ -f assembly/megahit_output/final.contigs.fa ]; then
        cp assembly/megahit_output/final.contigs.fa assembly/contigs.fasta
        
        # BLAST against viral database
        blastn -query assembly/contigs.fasta \
               -db /ref/sahlab/databases/viral_genomes/viral_genomes \
               -out assembly/blast_results.tsv \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
               -max_target_seqs 10 \
               -num_threads $THREADS
               
        # Generate assembly report
        echo 'Assembly Report' > assembly/assembly_report.txt
        echo 'Date: $(date)' >> assembly/assembly_report.txt
        echo 'Contigs: $(grep -c '^>' assembly/contigs.fasta)' >> assembly/assembly_report.txt
        echo '' >> assembly/assembly_report.txt
        echo 'Top BLAST hits:' >> assembly/assembly_report.txt
        head -10 assembly/blast_results.tsv >> assembly/assembly_report.txt
    else
        echo 'Assembly failed - no contigs produced' > assembly/assembly_report.txt
    fi
"

# Mark diagnostic complete
echo "$(date -Iseconds): Phase 1 diagnostic completed" > metadata/.phase1_diagnostic_complete
DIAGNOSTIC

chmod +x "scripts/structured_diagnostic_wrapper.sh"

# Submit both jobs
mkdir -p scripts
ASSEMBLY_JOB=$(sbatch --parsable scripts/structured_assembly_wrapper.sh)
DIAGNOSTIC_JOB=$(sbatch --parsable scripts/structured_diagnostic_wrapper.sh)

print_status "Jobs submitted:"
echo "  Assembly: $ASSEMBLY_JOB"
echo "  Diagnostic: $DIAGNOSTIC_JOB"

# Final setup
print_step 6 "Finalizing setup"

# Create Phase 2 parameters template
cat > "metadata/phase2_parameters_template.json" << PARAMS
{
    "quality_threshold": "TBD - Review assembly results",
    "depth_threshold": "TBD - Review depth plots",
    "frequency_threshold": "TBD - Review mutation distribution",
    "recommended_values": {
        "quality": "Use median+ from SNP statistics",
        "depth": "200-500 based on coverage plot",
        "frequency": "0.01 for minor variants, 0.05+ for major"
    }
}
PARAMS

echo ""
echo "======================================================================="
echo "                    STRUCTURED PIPELINE INITIATED"
echo "======================================================================="
echo ""
print_status "Analysis directory: $OUTPUT_DIR"
print_status "All modules will use standardized relative paths"
print_status "Jobs running: Assembly ($ASSEMBLY_JOB), Diagnostic ($DIAGNOSTIC_JOB)"
echo ""
print_info "Monitor progress:"
echo "  tail -f $OUTPUT_DIR/logs/assembly_$ASSEMBLY_JOB.out"
echo "  tail -f $OUTPUT_DIR/logs/diagnostic_$DIAGNOSTIC_JOB.out"
echo ""
print_info "When complete, run Phase 2:"
echo "  ${PIPELINE_BASE}/scripts/run_structured_phase2.sh -d $OUTPUT_DIR"
echo ""
