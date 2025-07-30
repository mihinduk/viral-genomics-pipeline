#\!/bin/bash
# Coordinate validation workflow
# Run after BLAST annotation to verify coordinate conservation

QUERY_ACC=$1
REF_ACC=$2
BLAST_OUTPUT=$3

echo "🔍 Validating coordinates for ${QUERY_ACC} against ${REF_ACC}"

# Run validation
python3 viral_pipeline/analysis/validate_blast_coordinates.py \
    ${QUERY_ACC} \
    ${REF_ACC} \
    --blast-output ${BLAST_OUTPUT} \
    --json-config viral_pipeline/visualization/known_viruses.json \
    --output-json validation_${QUERY_ACC}_vs_${REF_ACC}.json

# Check exit code
if [ $? -eq 0 ]; then
    echo "✅ Validation passed - coordinates are conserved"
    echo "📝 You can safely use ${REF_ACC} coordinates for ${QUERY_ACC}"
else
    echo "❌ Validation failed - please review mismatches"
    echo "📧 Contact Kathie to update JSON configuration if needed"
fi
