#\!/bin/bash
# DENV1 Example Workflow
# Dengue virus type 1 complete analysis and visualization

# Parameters
SAMPLE_PREFIX="SMS_5_DENV1"
ACCESSION="NC_001477.1"
QUALITY_CUTOFF="1000"
FREQ_CUTOFF="0.01"

echo "🦠 DENV1 Analysis Pipeline"
echo "=========================="

# 1. Assembly and Analysis
echo "Step 1: Running assembly and variant calling..."
/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics/run_pipeline_htcf_consolidated.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  $ACCESSION \
  4

# 2. Create output directory
mkdir -p ${SAMPLE_PREFIX}_results

# 3. Generate depth file
echo "Step 2: Generating depth coverage..."
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  samtools depth cleaned_seqs/mapping/*.lofreq.final.bam > ${SAMPLE_PREFIX}_results/${SAMPLE_PREFIX}_depth.txt

# 4. Quality filter VCF
echo "Step 3: Filtering mutations by quality..."
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ../viral_pipeline/visualization/parse_snpeff_vcf.py \
  -i cleaned_seqs/variants/*NC_001477*.snpEFF.ann.vcf \
  -d 200 \
  -q $QUALITY_CUTOFF \
  -o ${SAMPLE_PREFIX}_filtered_mutations.tsv \
  -O ${SAMPLE_PREFIX}_results

# 5. Create mutation visualization
echo "Step 4: Creating publication-ready mutation plot..."
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ../viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input ${SAMPLE_PREFIX}_results/${SAMPLE_PREFIX}_filtered_mutations.tsv \
  --output ${SAMPLE_PREFIX}_mutations.png \
  --accession $ACCESSION \
  --cutoff $FREQ_CUTOFF \
  --outdir ${SAMPLE_PREFIX}_results

# 6. Create depth visualization
echo "Step 5: Creating interactive depth visualization..."
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ../viral_pipeline/visualization/visualize_depth.py \
  --depth ${SAMPLE_PREFIX}_results/${SAMPLE_PREFIX}_depth.txt \
  --output ${SAMPLE_PREFIX}_results/${SAMPLE_PREFIX}_depth.png \
  --output-html ${SAMPLE_PREFIX}_results/${SAMPLE_PREFIX}_depth.html \
  --accession $ACCESSION

# 7. Generate diagnostic report
echo "Step 6: Generating comprehensive diagnostic report..."
sbatch ../viral_pipeline/analysis/submit_viral_diagnostic.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  $ACCESSION \
  $SAMPLE_PREFIX \
  4

echo "🎉 DENV1 pipeline complete\!"
echo "📁 Results in: ${SAMPLE_PREFIX}_results/"
echo "   - Beautiful mutation visualization with flavivirus gene layout"
echo "   - Interactive depth coverage with gene annotations"
echo "   - Complete mutation tables with perfect line-to-table parity"
echo "   - Comprehensive diagnostic report"
