#!/bin/bash
# Pre-Pipeline Virus Database Validation
# Checks SnpEff and visualization databases before pipeline execution
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
    echo -e "${MAGENTA}   VIRUS DATABASE VALIDATION${NC}"
    echo -e "${MAGENTA}========================================${NC}"
}

print_status() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_error() {
    echo -e "${RED}[✗]${NC} $1" >&2
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

print_info() {
    echo -e "${BLUE}[i]${NC} $1"
}

usage() {
    cat << USAGE
========================================================================
PRE-PIPELINE VIRUS DATABASE VALIDATION
========================================================================

This script validates that a virus accession is properly configured in
both SnpEff and visualization databases before pipeline execution.

CRITICAL: Run this BEFORE any pipeline modules to avoid failures.

Usage: $0 <accession> [options]

Arguments:
    accession        Virus accession to validate (e.g., NC_001477.1)

Options:
    -a               Auto-fix: Add to databases if missing (requires manual curation)
    -v               Verbose output
    -h               Show this help

Exit codes:
    0 - All databases ready, pipeline can proceed
    1 - Databases missing, need updates
    2 - Error in validation process

Example:
    # Check if DENV1 is ready for pipeline
    $0 NC_001477.1
    
    # Check and auto-add if missing
    $0 NEW_VIRUS.1 -a

USAGE
    exit 1
}

# Parse arguments
ACCESSION=$1
shift || { print_error "No accession provided"; usage; }

AUTO_FIX=0
VERBOSE=0

while getopts "avh" opt; do
    case $opt in
        a) AUTO_FIX=1 ;;
        v) VERBOSE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Pipeline paths
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
SNPEFF_CONFIG="/home/mihindu/software/snpEff/snpEff.config"
KNOWN_VIRUSES="${PIPELINE_BASE}/viral_pipeline/visualization/known_viruses.json"

print_header
echo ""
print_info "Validating databases for accession: $ACCESSION"
echo ""

# Validation results
SNPEFF_OK=0
VISUALIZATION_OK=0
NEEDS_SETUP=0

# Check 1: SnpEff Database
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1. SNPEFF DATABASE CHECK"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ ! -f "$SNPEFF_CONFIG" ]; then
    print_error "SnpEff config not found: $SNPEFF_CONFIG"
    exit 2
fi

if grep -q "${ACCESSION}.genome" "$SNPEFF_CONFIG"; then
    print_status "SnpEff database entry exists"
    SNPEFF_OK=1
    
    if [ $VERBOSE -eq 1 ]; then
        echo "   Entry: $(grep "${ACCESSION}.genome" "$SNPEFF_CONFIG")"
    fi
else
    print_warning "SnpEff database entry MISSING"
    SNPEFF_OK=0
    NEEDS_SETUP=1
fi

# Check 2: Visualization Database
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2. VISUALIZATION DATABASE CHECK"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ ! -f "$KNOWN_VIRUSES" ]; then
    print_error "Visualization database not found: $KNOWN_VIRUSES"
    exit 2
fi

# Check if accession exists in known_viruses.json
if python3 -c "
import json
import sys
with open('$KNOWN_VIRUSES', 'r') as f:
    data = json.load(f)
sys.exit(0 if '$ACCESSION' in data else 1)
" 2>/dev/null; then
    print_status "Visualization database entry exists"
    VISUALIZATION_OK=1
    
    if [ $VERBOSE -eq 1 ]; then
        python3 -c "
import json
with open('$KNOWN_VIRUSES', 'r') as f:
    data = json.load(f)
virus = data.get('$ACCESSION', {})
print(f'   Name: {virus.get(\"name\", \"Unknown\")}')
print(f'   Family: {virus.get(\"family\", \"Unknown\")}')
print(f'   Genus: {virus.get(\"genus\", \"Unknown\")}')
print(f'   Genes: {len(virus.get(\"gene_coords\", {}))} configured')
"
    fi
else
    print_warning "Visualization database entry MISSING"
    VISUALIZATION_OK=0
    NEEDS_SETUP=1
fi

# Check 3: Cross-reference validation
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3. CROSS-REFERENCE CHECK"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Look for related accessions (same virus, different strains)
python3 -c "
import json
accession = '$ACCESSION'

with open('$KNOWN_VIRUSES', 'r') as f:
    data = json.load(f)

# Find exact match
if accession in data:
    print('✓ Exact accession match found')
else:
    # Look for related accessions
    related = []
    for acc, info in data.items():
        name = info.get('name', '').lower()
        family = info.get('family', '')
        
        # Basic matching (can be enhanced)
        if any(virus in name for virus in ['dengue', 'west nile', 'zika', 'veev', 'powassan']):
            if any(virus in accession.lower() for virus in ['nc_001477', 'nc_009942', 'nc_012532']):
                related.append((acc, info.get('name', 'Unknown')))
    
    if related:
        print('Related accessions found:')
        for acc, name in related:
            print(f'  {acc}: {name}')
        print('⚠️  Consider if these share gene structure')
    else:
        print('No related accessions found')
" 2>/dev/null || echo "Could not perform cross-reference check"

# Summary and recommendations
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "4. VALIDATION SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ $SNPEFF_OK -eq 1 ] && [ $VISUALIZATION_OK -eq 1 ]; then
    echo ""
    print_status "ALL DATABASES READY - Pipeline can proceed"
    echo ""
    echo "✓ SnpEff database configured"
    echo "✓ Visualization database configured"
    echo ""
    print_info "Safe to run pipeline modules:"
    echo "  ./scripts/run_viral_pipeline_phase1.sh -s SAMPLE -a $ACCESSION ..."
    echo ""
    exit 0
else
    echo ""
    print_error "DATABASES INCOMPLETE - Setup required before pipeline"
    echo ""
    
    if [ $SNPEFF_OK -eq 0 ]; then
        echo "✗ SnpEff database missing"
    else
        echo "✓ SnpEff database ready"
    fi
    
    if [ $VISUALIZATION_OK -eq 0 ]; then
        echo "✗ Visualization database missing"
    else
        echo "✓ Visualization database ready"
    fi
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "REQUIRED ACTIONS:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    if [ $AUTO_FIX -eq 1 ]; then
        print_info "Auto-fix requested - initiating database setup"
        echo ""
        
        # Launch virus addition workflow
        print_info "Launching virus addition workflow..."
        exec ${PIPELINE_BASE}/scripts/add_new_virus.sh $ACCESSION
        
    else
        echo "1. Add virus to databases:"
        echo "   ./scripts/add_new_virus.sh $ACCESSION [options]"
        echo ""
        echo "2. For well-annotated viruses:"
        echo "   ./scripts/add_new_virus.sh $ACCESSION -G \\"
        echo "     -n \"Virus name\" -f Family -g Genus"
        echo ""
        echo "3. For BLAST annotation:"
        echo "   ./scripts/add_new_virus.sh $ACCESSION -r REF_ACCESSION \\"
        echo "     -n \"Virus name\" -f Family -g Genus"
        echo ""
        echo "4. After setup, re-run this validation:"
        echo "   $0 $ACCESSION"
        echo ""
    fi
    
    exit 1
fi
