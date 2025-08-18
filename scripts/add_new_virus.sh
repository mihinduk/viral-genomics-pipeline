#!/bin/bash
# Streamlined Virus Addition Helper
# Assists with adding new viruses to the pipeline database
# Author: Kathie Mihindu
# Date: 2025-08-18

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
NC='\033[0m'

print_header() {
    echo -e "${MAGENTA}========================================${NC}"
    echo -e "${MAGENTA}    VIRUS ADDITION HELPER TOOL${NC}"
    echo -e "${MAGENTA}========================================${NC}"
}

print_status() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_error() {
    echo -e "${RED}[✗]${NC} $1" >&2
}

print_info() {
    echo -e "${BLUE}[i]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

print_step() {
    echo -e "\n${YELLOW}Step $1:${NC} $2"
}

usage() {
    cat << USAGE

This tool streamlines the process of adding new viruses to the pipeline.

Usage: $0 <accession> [options]

Arguments:
    accession        GenBank/RefSeq accession (e.g., KU955591.1)

Options:
    -r REF_ACCESSION Reference virus for BLAST annotation
    -n NAME          Virus name (e.g., "Zika virus strain XYZ")
    -f FAMILY        Virus family (e.g., Flaviviridae, Togaviridae, Peribunyaviridae)
    -g GENUS         Virus genus (e.g., Orthoflavivirus, Alphavirus, Orthobunyavirus)
    -G               Use GenBank annotation (if available)
    -h               Show this help

Taxonomy Reference:
    Family: Flaviviridae
        Genus: Orthoflavivirus
            - Zika virus (NC_012532.1, NC_035889.1)
            - Dengue virus 1-4 (NC_001477.1)
            - West Nile virus (NC_009942.1)
            - Yellow fever virus (NC_002031.1)
            - Japanese encephalitis virus (NC_001437.1)
    
    Family: Togaviridae
        Genus: Alphavirus
            - Venezuelan equine encephalitis virus (NC_001449.1, NC_075022.1)
            - Chikungunya virus (NC_004162.2)
            - Eastern equine encephalitis virus (NC_003899.1)
            - Western equine encephalitis virus (NC_003908.1)
    
    Family: Peribunyaviridae
        Genus: Orthobunyavirus
            - Powassan virus (HM440560.1)
            - La Crosse virus (NC_004108.1)

Example workflows:

1. Add Flavivirus with GenBank annotation:
   $0 NC_035889.1 -G -n "Zika virus isolate PRVABC59" \\
      -f Flaviviridae -g Orthoflavivirus

2. Add Alphavirus using BLAST annotation against same species:
   $0 AY532665.1 -r NC_001563.2 -n "WNV strain NY99" \\
      -f Flaviviridae -g Orthoflavivirus

3. Add virus using cross-species BLAST (requires careful review):
   $0 NEW_VIRUS.1 -r NC_001477.1 -n "Novel flavivirus" \\
      -f Flaviviridae -g Orthoflavivirus

USAGE
    exit 1
}

# Parse arguments
ACCESSION=$1
shift

if [ -z "$ACCESSION" ]; then
    print_error "No accession provided"
    usage
fi

# Default values
USE_GENBANK=0
REF_ACCESSION=""
VIRUS_NAME=""
VIRUS_FAMILY=""
VIRUS_GENUS=""

# Parse options
while getopts "r:n:f:g:Gh" opt; do
    case $opt in
        r) REF_ACCESSION="$OPTARG" ;;
        n) VIRUS_NAME="$OPTARG" ;;
        f) VIRUS_FAMILY="$OPTARG" ;;
        g) VIRUS_GENUS="$OPTARG" ;;
        G) USE_GENBANK=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Pipeline paths
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"

print_header
echo ""
print_info "Adding virus: $ACCESSION"
[ -n "$VIRUS_NAME" ] && print_info "Name: $VIRUS_NAME"
[ -n "$VIRUS_FAMILY" ] && print_info "Family: $VIRUS_FAMILY"
[ -n "$VIRUS_GENUS" ] && print_info "Genus: $VIRUS_GENUS"
echo ""

# Validate taxonomy if provided
if [ -n "$VIRUS_FAMILY" ] && [ -n "$VIRUS_GENUS" ]; then
    print_info "Validating taxonomy..."
    case "$VIRUS_FAMILY" in
        Flaviviridae)
            if [ "$VIRUS_GENUS" != "Orthoflavivirus" ]; then
                print_error "For family Flaviviridae, genus should be Orthoflavivirus"
                print_info "Common error: 'flavivirus' is not a genus, it's a common name"
                exit 1
            fi
            ;;
        Togaviridae)
            if [ "$VIRUS_GENUS" != "Alphavirus" ]; then
                print_error "For family Togaviridae, genus should be Alphavirus"
                exit 1
            fi
            ;;
        Peribunyaviridae)
            if [ "$VIRUS_GENUS" != "Orthobunyavirus" ]; then
                print_error "For family Peribunyaviridae, expected genus Orthobunyavirus"
                print_info "Note: Powassan virus is in genus Orthobunyavirus"
            fi
            ;;
        *)
            print_info "Unrecognized family: $VIRUS_FAMILY"
            print_info "Ensure taxonomy is correct before proceeding"
            ;;
    esac
    print_status "Taxonomy validation passed"
fi

# Step 1: Check if virus already exists
print_step 1 "Checking if virus is already in database"

if $MAMBA_CMD python3 -c "
import json
with open('${PIPELINE_BASE}/viral_pipeline/visualization/known_viruses.json', 'r') as f:
    data = json.load(f)
    if '$ACCESSION' in data:
        print('Virus already exists in database')
        exit(1)
" 2>/dev/null; then
    print_status "Virus not in database, proceeding"
else
    print_error "Virus $ACCESSION already exists in known_viruses.json"
    print_info "To update existing virus, edit known_viruses.json directly"
    exit 1
fi

# Step 2: Download sequence
print_step 2 "Downloading sequence from NCBI"
cd $PIPELINE_BASE

$MAMBA_CMD python3 -c "
from Bio import Entrez, SeqIO
Entrez.email = 'your.email@example.com'

try:
    handle = Entrez.efetch(db='nucleotide', id='$ACCESSION', rettype='fasta', retmode='text')
    with open('${ACCESSION}.fasta', 'w') as f:
        f.write(handle.read())
    handle.close()
    print('Sequence downloaded successfully')
except Exception as e:
    print(f'Error downloading sequence: {e}')
    exit(1)
"

if [ -f "${ACCESSION}.fasta" ]; then
    print_status "Sequence downloaded: ${ACCESSION}.fasta"
    SEQ_LENGTH=$(grep -v '>' ${ACCESSION}.fasta | tr -d '\n' | wc -c)
    print_info "Sequence length: $SEQ_LENGTH bp"
else
    print_error "Failed to download sequence"
    exit 1
fi

# Step 3: Run annotation
if [ $USE_GENBANK -eq 1 ]; then
    print_step 3 "Running GenBank annotation"
    print_info "Submitting GenBank annotation job..."
    
    JOB_ID=$(sbatch --parsable ${PIPELINE_BASE}/viral_pipeline/analysis/submit_genbank_annotation.sh $ACCESSION)
    print_status "Job submitted: $JOB_ID"
    
    ANNOTATION_TYPE="GenBank"
else
    if [ -z "$REF_ACCESSION" ]; then
        print_error "Reference accession required for BLAST annotation"
        print_info "Use -r option to specify reference virus"
        usage
    fi
    
    print_step 3 "Running BLAST annotation"
    print_info "Using reference: $REF_ACCESSION"
    
    JOB_ID=$(sbatch --parsable ${PIPELINE_BASE}/viral_pipeline/analysis/submit_blast_annotation.sh $ACCESSION $REF_ACCESSION)
    print_status "Job submitted: $JOB_ID"
    
    ANNOTATION_TYPE="BLAST vs $REF_ACCESSION"
fi

# Step 4: Monitor job
print_step 4 "Monitoring annotation job"
print_info "Waiting for job $JOB_ID to complete..."

while squeue -j $JOB_ID 2>/dev/null | grep -q $JOB_ID; do
    sleep 10
done

if sacct -j $JOB_ID --format=State --noheader | grep -q COMPLETED; then
    print_status "Annotation completed successfully"
else
    print_error "Annotation job failed. Check logs:"
    echo "  sbatch output: slurm-${JOB_ID}.out"
    exit 1
fi

# Step 5: Check annotation results and compare with reference
print_step 5 "Reviewing annotation results"

if [ $USE_GENBANK -eq 1 ]; then
    ANNOTATION_DIR="genbank_annotation_${ACCESSION}"
else
    ANNOTATION_DIR="blast_annotation_${ACCESSION}"
fi

if [ -d "$ANNOTATION_DIR" ]; then
    print_status "Annotation directory found: $ANNOTATION_DIR"
    echo ""
    print_info "Contents:"
    ls -la $ANNOTATION_DIR/ | tail -n +2
    echo ""
    
    # Display annotation summary if available
    if [ -f "${ANNOTATION_DIR}/${ACCESSION}_annotation_summary.json" ]; then
        print_info "Annotation summary:"
        $MAMBA_CMD python3 -c "
import json
with open('${ANNOTATION_DIR}/${ACCESSION}_annotation_summary.json', 'r') as f:
    data = json.load(f)
    for key, value in data.items():
        print(f'  {key}: {value}')
"
    fi
else
    print_error "Annotation directory not found"
    exit 1
fi

# Step 6: Compare with reference virus if using BLAST
if [ ! $USE_GENBANK -eq 1 ] && [ -n "$REF_ACCESSION" ]; then
    print_step 6 "Comparing annotations with reference virus"
    
    $MAMBA_CMD python3 -c "
import json
import os

ref_accession = '$REF_ACCESSION'
target_accession = '$ACCESSION'
annotation_dir = '$ANNOTATION_DIR'

# Load known_viruses.json to get reference structure
with open('${PIPELINE_BASE}/viral_pipeline/visualization/known_viruses.json', 'r') as f:
    known_viruses = json.load(f)

if ref_accession in known_viruses:
    ref_data = known_viruses[ref_accession]
    print(f'\nReference virus: {ref_data.get(\"name\", ref_accession)}')
    print(f'Family: {ref_data.get(\"family\", \"unknown\")}')
    print(f'Genus: {ref_data.get(\"genus\", \"unknown\")}')
    print(f'Genome length: {ref_data.get(\"genome_length\", \"unknown\")} bp')
    
    # Show reference gene structure
    print('\nReference gene structure:')
    ref_genes = ref_data.get('gene_coords', {})
    for gene, coords in ref_genes.items():
        print(f'  {gene}: {coords[0]}-{coords[1]} ({coords[1]-coords[0]+1} bp)')
    
    # Try to load annotation results
    annotation_file = os.path.join(annotation_dir, f'{target_accession}_annotation_summary.json')
    if os.path.exists(annotation_file):
        with open(annotation_file, 'r') as f:
            target_data = json.load(f)
        
        print('\n⚠️  ANNOTATION COMPARISON:')
        print('Reference genes:', list(ref_genes.keys()))
        target_genes = target_data.get('genes', [])
        print('Annotated genes:', target_genes)
        
        # Check for missing or extra genes
        ref_gene_set = set(ref_genes.keys())
        target_gene_set = set(target_genes) if target_genes else set()
        
        missing = ref_gene_set - target_gene_set
        extra = target_gene_set - ref_gene_set
        
        if missing:
            print(f'\n⚠️  Missing genes in annotation: {list(missing)}')
            print('   These may need manual adjustment')
        if extra:
            print(f'\n⚠️  Extra genes in annotation: {list(extra)}')
            print('   Verify these are correct for your strain')
        
        if not missing and not extra:
            print('\n✓ Gene names match reference')
        
        # Check genome length difference
        target_length = int('$SEQ_LENGTH')
        ref_length = ref_data.get('genome_length', 0)
        length_diff = abs(target_length - ref_length)
        if length_diff > 100:
            print(f'\n⚠️  Significant genome length difference: {length_diff} bp')
            print(f'   Reference: {ref_length} bp, Target: {target_length} bp')
            print('   Gene coordinates may need adjustment')
    else:
        print('\nNo detailed annotation summary found')
else:
    print(f'\n⚠️  Reference {ref_accession} not in known_viruses.json')
    print('   Cannot perform detailed comparison')
    print('   Manual review of annotations strongly recommended')
" || print_warning "Could not perform annotation comparison"
fi

# Step 7: Manual curation reminder
print_step 7 "Manual Curation Required"
echo ""
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}  MANUAL CURATION STEPS:${NC}"
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "1. Review the annotation files in: $ANNOTATION_DIR/"
echo "   - Check gene coordinates in *.gff3 file"
echo "   - Verify gene names and boundaries"
echo "   - Confirm polyprotein cleavage sites (if applicable)"

if [ ! $USE_GENBANK -eq 1 ] && [ -n "$REF_ACCESSION" ]; then
    echo ""
    echo "2. Address annotation differences identified above:"
    echo "   - Missing genes may indicate annotation transfer issues"
    echo "   - Extra genes may be strain-specific features"
    echo "   - Length differences may require coordinate adjustments"
fi

echo ""
echo "3. Update known_viruses.json with:"
echo "   - Virus name, family, AND genus"
echo "   - Gene coordinates and names"
echo "   - Gene colors for visualization"
echo "   - Structural vs non-structural gene classification"
echo ""
echo "4. Create SnpEff database:"
echo "   bash ${PIPELINE_BASE}/viral_pipeline/utils/update_snpeff_for_sample.sh $ACCESSION"
echo ""
echo "5. Test the new configuration:"
echo "   - Run a test sample through the pipeline"
echo "   - Verify visualizations show correct genes"
echo "   - Check that mutations are annotated correctly"
echo ""

# Step 8: Generate template entry
print_step 8 "Generating template for known_viruses.json"

# Generate template based on annotation results or reference
$MAMBA_CMD python3 -c "
import json
import os

accession = '$ACCESSION'
ref_accession = '$REF_ACCESSION'
annotation_dir = '$ANNOTATION_DIR'
virus_name = '$VIRUS_NAME' or f'{accession} (annotated by $ANNOTATION_TYPE)'
virus_family = '$VIRUS_FAMILY' or 'unknown'
virus_genus = '$VIRUS_GENUS' or 'unknown'
seq_length = int('$SEQ_LENGTH')

template = {
    accession: {
        'name': virus_name,
        'family': virus_family,
        'genus': virus_genus,
        'genome_length': seq_length,
        'gene_coords': {},
        'structural_genes': [],
        'nonstructural_genes': [],
        'colors': {}
    }
}

# Try to populate from reference if available
if ref_accession:
    with open('${PIPELINE_BASE}/viral_pipeline/visualization/known_viruses.json', 'r') as f:
        known_viruses = json.load(f)
    
    if ref_accession in known_viruses:
        ref_data = known_viruses[ref_accession]
        
        # Copy structure from reference as starting point
        template[accession]['gene_coords'] = ref_data.get('gene_coords', {}).copy()
        template[accession]['structural_genes'] = ref_data.get('structural_genes', []).copy()
        template[accession]['nonstructural_genes'] = ref_data.get('nonstructural_genes', []).copy()
        template[accession]['colors'] = ref_data.get('colors', {}).copy()
        
        # Add TODO comments
        template[accession]['_TODO'] = [
            'Verify gene coordinates match your strain',
            'Adjust for genome length differences',
            'Check polyprotein cleavage sites',
            'Update colors if desired'
        ]

# Write template
with open(f'{accession}_template.json', 'w') as f:
    json.dump(template, f, indent=4)

print(f'Template saved to: {accession}_template.json')
" || echo "Could not generate template automatically"

if [ -f "${ACCESSION}_template.json" ]; then
    print_status "Template saved to: ${ACCESSION}_template.json"
    echo ""
    print_info "Review and edit this template before adding to known_viruses.json"
    
    # Show template content
    echo ""
    echo "Template preview:"
    head -20 ${ACCESSION}_template.json
    echo "..."
fi

# Final summary
echo ""
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}  VIRUS ADDITION INITIATED${NC}"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
print_status "Annotation completed using $ANNOTATION_TYPE"
print_status "Files ready for manual curation in: $ANNOTATION_DIR/"
print_status "Template configuration in: ${ACCESSION}_template.json"
[ -n "$VIRUS_FAMILY" ] && print_status "Family: $VIRUS_FAMILY"
[ -n "$VIRUS_GENUS" ] && print_status "Genus: $VIRUS_GENUS"
echo ""

if [ ! $USE_GENBANK -eq 1 ] && [ -n "$REF_ACCESSION" ]; then
    print_warning "BLAST annotation requires careful review of gene transfers"
    print_info "Pay special attention to any differences noted above"
fi

print_info "Next: Complete manual curation steps above"
echo ""

