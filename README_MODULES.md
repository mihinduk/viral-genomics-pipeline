# Viral Genomics Pipeline - Complete Module Guide

## 🚀 Quick Start

This pipeline provides end-to-end viral genome analysis with publication-ready visualizations. Each module can be run independently or as part of the complete workflow.

## 📋 Prerequisites

```bash
# Activate conda environment
source /ref/sahlab/software/anaconda3/bin/activate
conda activate viral_genomics

# Set pipeline path for easy reference
PIPELINE_DIR="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
```

## 🎯 Complete Workflow Example: ZIKV Analysis (SMS_9)

This real-world example shows the complete analysis of a Zika virus sample with high-quality results:

### Module 1: Reference-guided Assembly and SNP Calling
```bash
# Submit the main pipeline job
sbatch /scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/analysis/submit_viral_pipeline_enhanced.sh \
      "./NovaSeq_N1027_I14068_Cell_Culture_RNA_Diamond_Scheaffer_SMS_9_ZIKV_Dakar_41525_MA_ZIKV_25.256_Vero_R1.fastq.gz" \
      "./NovaSeq_N1027_I14068_Cell_Culture_RNA_Diamond_Scheaffer_SMS_9_ZIKV_Dakar_41525_MA_ZIKV_25.256_Vero_R2.fastq.gz" \
      KU955591.1 \
      4 \
      --large-files

# Monitor job status
squeue -u $USER
# Example output: Submitted batch job 25780118
```

**What this does:**
- Quality control and trimming
- Reference-guided assembly
- SNP calling with LoFreq
- SnpEff annotation
- Generates: BAM files, VCF files, annotated TSV

### Module 2: Create Output Directory & Generate Depth File
```bash
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics bash -c '
  mkdir -p SMS_9_ZIKV_results &&
  echo -e "chrom\tposition\tdepth" > SMS_9_ZIKV_results/SMS_9_ZIKV_depth.txt &&
  samtools depth cleaned_seqs/mapping/NovaSeq_N1027_I14068_Cell_Culture_RNA_Diamond_Scheaffer_SMS_9_ZIKV_Dakar_41525_MA_ZIKV_25.256_Vero.lofreq.final.bam >> SMS_9_ZIKV_results/SMS_9_ZIKV_depth.txt'
```

**What this does:**
- Creates organized output directory
- Extracts read depth at every position
- Prepares data for visualization

### Module 3: Parse & Filter Mutations
```bash
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 /scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/visualization/parse_snpeff_tsv.py \
  "cleaned_seqs/variants/NovaSeq_N1027_I14068_Cell_Culture_RNA_Diamond_Scheaffer_SMS_9_ZIKV_Dakar_41525_MA_ZIKV_25.256_Vero.snpEFF.ann.tsv" \
  "SMS_9_ZIKV_results/SMS_9_ZIKV_filtered_mutations.tsv" \
  --quality 42814 \
  --depth 200 \
  --freq 0.01
```

**Parameters explained:**
- `--quality 42814`: Minimum variant quality score
- `--depth 200`: Minimum read depth required
- `--freq 0.01`: Minimum allele frequency (1%)

### Module 4: Mutation Visualization
```bash
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 /scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input SMS_9_ZIKV_results/SMS_9_ZIKV_filtered_mutations.tsv \
  --output SMS_9_ZIKV_mutations.png \
  --accession KU955591.1 \
  --cutoff 0.01 \
  --outdir SMS_9_ZIKV_results
```

**Output:** Publication-ready mutation visualization with:
- Color-coded genes (structural vs non-structural)
- Clean gene labels (C, ancC, pr, prM, M, E, NS1, etc.)
- Mutation tables by gene
- Consistent styling

### Module 5: Depth Visualization
```bash
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 /scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/visualization/visualize_depth.py \
  --depth SMS_9_ZIKV_results/SMS_9_ZIKV_depth.txt \
  --output SMS_9_ZIKV_results/SMS_9_ZIKV_depth.png \
  --output-html SMS_9_ZIKV_results/SMS_9_ZIKV_depth.html \
  --accession KU955591.1
```

**Output:** Professional depth visualization with:
- Color-coded gene regions
- Coverage statistics
- Genomically ordered legend
- Interactive HTML version

### Module 6: Diagnostic Report (Optional)
```bash
sbatch /scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/analysis/submit_viral_diagnostic.sh \
  "./NovaSeq_N1027_I14068_Cell_Culture_RNA_Diamond_Scheaffer_SMS_9_ZIKV_Dakar_41525_MA_ZIKV_25.256_Vero_R1.fastq.gz" \
  "./NovaSeq_N1027_I14068_Cell_Culture_RNA_Diamond_Scheaffer_SMS_9_ZIKV_Dakar_41525_MA_ZIKV_25.256_Vero_R2.fastq.gz" \
  KU955591.1 \
  diagnostic_SMS_9 \
  4
```

**Output:** Comprehensive diagnostic report with assembly metrics, coverage analysis, and quality assessments.

### Module 7: Generate Filtered Consensus
```bash
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
    python3 /scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/utils/generate_filtered_consensus.py \
    --vcf SMS_9_ZIKV_results/SMS_9_ZIKV_filtered_mutations.tsv \
    --reference cleaned_seqs/KU955591.1.fasta \
    --accession KU955591.1 \
    --quality 42814 \
    --depth 200 \
    --freq 0.50 \
    --output-prefix SMS_9_ZIKV_results/SMS_9_ZIKV
```

**Output:** High-quality consensus genome with mutations incorporated.

### Module 8: Download All Results
```bash
# Run from your local machine
cd "/Users/handley_lab/Handley Lab Dropbox/virome/Diamond_lab_isolate_seq/2025_07_03_Diamond_NovaSeq_N1027/ZIKV_Dakar"

# Download main analysis results
rsync -avh htcf:/scratch/sahlab/kathie/NovaSeq_N1027_ZIKV_Dakar/SMS_9_ZIKV_results/* .

# Download diagnostic results (if Module 6 was run)
for ext in fa tsv txt html; do rsync -avh htcf:/scratch/sahlab/kathie/NovaSeq_N1027_ZIKV_Dakar/diagnostic_SMS_9/*.$ext . ; done
```

**Logical Flow:** Generate all files on HTCF first (Modules 6-7), then download everything at once (Module 8).

---

## 🎯 Quick Workflow Templates

### Template 1: Standard Analysis
```bash
# Replace with your values:
SAMPLE_NAME="your_sample"
R1_FASTQ="./your_sample_R1.fastq.gz"  
R2_FASTQ="./your_sample_R2.fastq.gz"
ACCESSION="NC_001477.1"  # or KU955591.1 for ZIKV, HM440560.1 for POWV
THREADS=4

# Module 1: Main pipeline
sbatch ${PIPELINE_DIR}/viral_pipeline/analysis/submit_viral_pipeline_enhanced.sh \
  "${R1_FASTQ}" "${R2_FASTQ}" "${ACCESSION}" "${THREADS}"

# Module 2: Depth extraction (after job completes)
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics bash -c "
  mkdir -p ${SAMPLE_NAME}_results &&
  echo -e 'chrom\tposition\tdepth' > ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt &&
  samtools depth cleaned_seqs/mapping/*.lofreq.final.bam >> ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt"

# Module 3: Parse mutations
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/parse_snpeff_tsv.py \
  "cleaned_seqs/variants/*.snpEFF.ann.tsv" \
  "${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv" \
  --quality 1000 --depth 200 --freq 0.01

# Module 4: Mutation visualization
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \
  --output ${SAMPLE_NAME}_mutations.png \
  --accession ${ACCESSION} --cutoff 0.01 \
  --outdir ${SAMPLE_NAME}_results

# Module 5: Depth visualization  
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_depth.py \
  --depth ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt \
  --output ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.png \
  --output-html ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.html \
  --accession ${ACCESSION}

# Module 6: Diagnostic report (optional)
sbatch ${PIPELINE_DIR}/viral_pipeline/analysis/submit_viral_diagnostic.sh \
  "${R1_FASTQ}" "${R2_FASTQ}" "${ACCESSION}" "diagnostic_${SAMPLE_NAME}" "${THREADS}"

# Module 7: Generate consensus (optional)
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/utils/generate_filtered_consensus.py \
  --vcf ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \
  --reference cleaned_seqs/${ACCESSION}.fasta \
  --accession ${ACCESSION} --quality 1000 --depth 200 --freq 0.50 \
  --output-prefix ${SAMPLE_NAME}_results/${SAMPLE_NAME}

# Module 8: Download all results (run from local machine)
# rsync -avh htcf:/path/to/your/analysis/${SAMPLE_NAME}_results/* .
# for ext in fa tsv txt html; do rsync -avh htcf:/path/to/your/analysis/diagnostic_${SAMPLE_NAME}/*.$ext . ; done
```

### Template 2: High-Depth Analysis (for samples >10,000x coverage)
```bash
# Use higher quality thresholds for high-depth samples
QUALITY_THRESHOLD=50000
DEPTH_THRESHOLD=1000
FREQ_THRESHOLD=0.001  # 0.1% for high sensitivity

# Same pipeline, different filtering parameters
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/parse_snpeff_tsv.py \
  "cleaned_seqs/variants/*.snpEFF.ann.tsv" \
  "${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv" \
  --quality ${QUALITY_THRESHOLD} \
  --depth ${DEPTH_THRESHOLD} \
  --freq ${FREQ_THRESHOLD}
```

---

## 📊 Supported Virus References

Our curated reference database provides optimal visualization for:

| Virus Family | Accession | Virus | Quality Score | Status |
|-------------|-----------|-------|---------------|---------|
| **Flaviviridae** | KU955591.1 | Zika virus (Dakar) | ⭐ 10/10 | Perfect |
| | NC_012532.1 | Zika virus (MR766) | 8/10 | Good |
| | NC_075022.1 | West Nile virus | 9/10 | Excellent |
| **Flaviviridae (Tick-borne)** | HM440560.1 | Powassan virus | 7/10 | Adequate |
| **Togaviridae** | NC_075022.1 | Venezuelan EEV | 8/10 | Good |

**Quality 10/10** references provide perfect gene label positioning with no overlapping text.

---

## 🔧 Common Parameter Adjustments

### Quality Filtering
```bash
# Strict filtering (high-quality samples)
--quality 50000 --depth 1000 --freq 0.001

# Standard filtering (most samples)  
--quality 1000 --depth 200 --freq 0.01

# Permissive filtering (low-depth samples)
--quality 500 --depth 50 --freq 0.05
```

### Visualization Options
```bash
# Focus on specific genes
--mutation-genes "structural"           # Only structural proteins
--mutation-genes "non-structural"       # Only NS proteins  
--mutation-genes "NS3,NS5"             # Specific genes

# Adjust frequency cutoffs
--cutoff 0.01    # Show variants ≥1%
--cutoff 0.10    # Show variants ≥10% 
--cutoff 0.50    # Show variants ≥50%
```

---

## 🚨 Troubleshooting

### Common Issues:

**1. File not found errors:**
```bash
# Check file paths
ls -la ./your_sample_R*.fastq.gz
ls -la cleaned_seqs/variants/*.tsv
```

**2. Low mutation counts:**
```bash
# Reduce filtering stringency
--quality 500 --depth 50 --freq 0.01
```

**3. Pipeline job failures:**
```bash
# Check SLURM logs
cat slurm-*.out
tail -50 slurm-*.out
```

**4. Missing gene annotations:**
```bash
# Update SnpEff database for new viruses
bash ${PIPELINE_DIR}/viral_pipeline/utils/update_snpeff_for_sample.sh YOUR_ACCESSION
```

---

## 📚 Output Files

After completing all modules, you'll have:

```
sample_results/
├── sample_mutations.png          # Publication-ready mutation plot
├── sample_depth.png              # Coverage visualization  
├── sample_depth.html             # Interactive depth plot
├── sample_filtered_mutations.tsv # Filtered variant calls
├── sample_mutations_table.tsv    # Mutation summary table
├── sample_consensus.fasta         # Consensus genome (Module 8)
└── sample_proteins.fasta          # Protein sequences (Module 8)
```

---

## 🎯 Next Steps

1. **Customize for your virus**: Add new references using our database system
2. **Batch processing**: Scale up for multiple samples  
3. **Publication**: Use high-quality visualizations in manuscripts
4. **Collaboration**: Share standardized results with colleagues

For questions or new virus requests, open a GitHub issue at: `https://github.com/mihinduk/viral-genomics-pipeline/issues`
## 🎯 Latest Enhancements (January 2025)

### Enhanced Mutation Visualization with Comprehensive Classification

Our mutation visualization script now includes advanced features for publication-ready figures:

#### **Comprehensive Mutation Type Support**
- **High Impact**: Stop gained/lost, start lost, frameshift, splice variants
- **Moderate Impact**: Missense, in-frame indels  
- **Low Impact**: Synonymous variants
- **Modifier**: Regulatory (UTR), intronic, intergenic variants

#### **Visual Priority System**
When multiple mutations occur at nearby positions (e.g., affecting the same codon), the visualization automatically prioritizes non-synonymous mutations over synonymous ones. This ensures that functionally important mutations are always visible.

**Example**: If position 6965 has both a missense mutation (G>C, Gly18Arg) and a synonymous mutation (A>C, Gly18Gly), the visualization will show the **missense color (orange)** instead of being masked by the synonymous color (green).

#### **Dynamic Legend Generation**
The legend automatically adjusts to show only the mutation types present in your specific dataset, reducing visual clutter and focusing on relevant variants.

#### **Enhanced Mutation Tables**
Output TSV files now include:
- Mutation_Type: Human-readable classification (e.g., Missense, Nonsense/Stop Gained)
- Complete functional annotation from SnpEff
- Frequency and impact information

### Tested References with Perfect Visualization

 < /dev/null |  Virus | Accession | Status | Visualization Quality |
|-------|-----------|--------|----------------------|
| **Zika virus (Dakar)** | KU955591.1 | ⭐ Perfect | Gene labels optimally positioned |
| **Powassan virus** | HM440560.1 | ✅ Enhanced | Comprehensive mutation classification |
| **West Nile virus** | NC_075022.1 | ✅ Good | Standard visualization |

---

*Last updated: January 2025 with comprehensive mutation type classification and visual priority system*
