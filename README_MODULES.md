# Viral Genomics Pipeline - Complete Module Guide

## Initial Setup

```bash
# Clone the pipeline to your working directory
cd /scratch/sahlab/kathie/
git clone https://github.com/mihinduk/viral-genomics-pipeline.git

# Set pipeline path for easy reference
PIPELINE_DIR="/scratch/sahlab/kathie/viral-genomics-pipeline"

# Activate conda environment
source /ref/sahlab/software/anaconda3/bin/activate
conda activate viral_genomics
```

## Complete Workflow Examples

### A. Standard Sample (Most common)

For typical viral samples with moderate sequencing depth:

```bash
# Navigate to your analysis directory
cd /scratch/sahlab/kathie/my_analysis

# Module 1: Assembly & Analysis
/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics/run_pipeline_htcf_consolidated.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  NC_001477.1 \
  4

# Module 2: Create Output Directory & Generate Depth File
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics bash -c '
  mkdir -p sample_results &&
  echo -e "chrom\tposition\tdepth" > sample_results/sample_depth.txt &&
  samtools depth cleaned_seqs/mapping/*.lofreq.final.bam >> sample_results/sample_depth.txt'

# Module 3: Parse & Filter Mutations
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/parse_snpeff_tsv.py \
  "cleaned_seqs/variants/*.snpEFF.ann.tsv" \
  "sample_results/sample_filtered_mutations.tsv" \
  --quality 1000 \
  --depth 200 \
  --freq 0.01

# Module 4: Mutation Visualization
# NOTE: Use ONLY the filename with --output when using --outdir
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input sample_results/sample_filtered_mutations.tsv \
  --output sample_mutations.png \
  --accession NC_001477.1 \
  --cutoff 0.01 \
  --outdir sample_results

# Module 5: Depth Visualization
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_depth.py \
  --depth sample_results/sample_depth.txt \
  --output sample_results/sample_depth.png \
  --output-html sample_results/sample_depth.html \
  --accession NC_001477.1

# Module 6: Download Results (run from local machine)
rsync -avh htcf:/scratch/sahlab/kathie/my_analysis/sample_results/* ./local_results/

# Module 7: Diagnostic Report (Optional)
sbatch ${PIPELINE_DIR}/viral_pipeline/analysis/submit_viral_diagnostic.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  NC_001477.1 \
  sample \
  4

# Module 8: Generate Filtered Consensus (Optional)
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/utils/generate_filtered_consensus.py \
  --vcf sample_results/sample_filtered_mutations.tsv \
  --reference /path/to/reference.fasta \
  --accession NC_001477.1 \
  --quality 1000 \
  --depth 200 \
  --freq 0.05 \
  --output-prefix sample_results/sample_consensus

# Module 10: BLAST Gene Annotation (If needed for polyprotein viruses)
# Only run if you see generic gene names like "Gene_106_10377"
sbatch ${PIPELINE_DIR}/viral_pipeline/analysis/submit_blast_annotation.sh \
  NC_001477.1 \
  NC_001477.1

# Update SnpEff database
bash ${PIPELINE_DIR}/viral_pipeline/utils/update_snpeff_for_sample.sh NC_001477.1
```

### B. Large Sample (High depth sequencing)

For samples with high sequencing depth or complex populations:

```bash
# Module 1: Assembly & Analysis with large file support
sbatch /scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics/submit_viral_pipeline_large.sh \
  "./large_sample_R1.fastq.gz" \
  "./large_sample_R2.fastq.gz" \
  NC_009942.1 \
  8 \
  --large-files

# Module 2: Same as standard
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics bash -c '
  mkdir -p large_sample_results &&
  echo -e "chrom\tposition\tdepth" > large_sample_results/large_sample_depth.txt &&
  samtools depth cleaned_seqs/mapping/*.lofreq.final.bam >> large_sample_results/large_sample_depth.txt'

# Module 3: Parse with stricter filtering for high-depth data
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/parse_snpeff_tsv.py \
  "cleaned_seqs/variants/*.snpEFF.ann.tsv" \
  "large_sample_results/large_sample_filtered_mutations.tsv" \
  --quality 10000 \
  --depth 500 \
  --freq 0.05

# Module 4: Visualize with higher cutoff
# NOTE: Use ONLY the filename with --output when using --outdir
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input large_sample_results/large_sample_filtered_mutations.tsv \
  --output large_sample_mutations.png \
  --accession NC_009942.1 \
  --cutoff 0.10 \
  --outdir large_sample_results

# Module 5: Depth visualization (use default log scale for high depth)
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_depth.py \
  --depth large_sample_results/large_sample_depth.txt \
  --output large_sample_results/large_sample_depth.png \
  --output-html large_sample_results/large_sample_depth.html \
  --accession NC_009942.1

# Modules 6-8: Same as standard
```

### C. Extremely Large Sample (Ultra-deep sequencing)

For samples with extremely high coverage (>100,000x) or very large files:

```bash
# Module 1: Use extreme pipeline
sbatch /scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics/submit_viral_pipeline_extreme.sh \
  "./extreme_sample_R1.fastq.gz" \
  "./extreme_sample_R2.fastq.gz" \
  AY532665.1 \
  16 \
  --extremely-large-files

# Module 2: Same as standard
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics bash -c '
  mkdir -p extreme_sample_results &&
  echo -e "chrom\tposition\tdepth" > extreme_sample_results/extreme_sample_depth.txt &&
  samtools depth cleaned_seqs/mapping/*.lofreq.final.bam >> extreme_sample_results/extreme_sample_depth.txt'

# Module 3: Very strict filtering for ultra-deep data
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/parse_snpeff_tsv.py \
  "cleaned_seqs/variants/*.snpEFF.ann.tsv" \
  "extreme_sample_results/extreme_sample_filtered_mutations.tsv" \
  --quality 49314 \
  --depth 1000 \
  --freq 0.10

# Module 4: Visualize only high-frequency variants
# NOTE: Use ONLY the filename with --output when using --outdir
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input extreme_sample_results/extreme_sample_filtered_mutations.tsv \
  --output extreme_sample_mutations.png \
  --accession AY532665.1 \
  --cutoff 0.20 \
  --outdir extreme_sample_results

# Module 5: Use linear scale if needed to see actual depth values
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/visualization/visualize_depth.py \
  --depth extreme_sample_results/extreme_sample_depth.txt \
  --output extreme_sample_results/extreme_sample_depth.png \
  --output-html extreme_sample_results/extreme_sample_depth.html \
  --accession AY532665.1 \
  --linear-scale

# Module 8: Generate high-confidence consensus (very important for noisy data)
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_DIR}/viral_pipeline/utils/generate_filtered_consensus.py \
  --vcf extreme_sample_results/extreme_sample_filtered_mutations.tsv \
  --reference AY532665.1.fasta \
  --accession AY532665.1 \
  --quality 49314 \
  --depth 5000 \
  --freq 0.90 \
  --output-prefix extreme_sample_results/extreme_sample_consensus
```

## Parameter Guidelines

### Quality Cutoffs (--quality)
- **Standard**: 1000
- **Large**: 10000  
- **Extreme**: 49314

### Depth Cutoffs (--depth)
- **Standard**: 200
- **Large**: 500
- **Extreme**: 1000-5000

### Frequency Cutoffs
- **Parse (--freq)**: 0.01-0.10 (keep more for flexibility)
- **Visualize (--cutoff)**: 0.01-0.50 (adjust for clean plots)
- **Consensus (--freq)**: 0.50-0.95 (only major variants)

### Threads
- **Standard**: 4
- **Large**: 8
- **Extreme**: 16

## Common Viruses and Accessions

| Virus | Accession | Genome Length |
|-------|-----------|---------------|
| Dengue virus 1 | NC_001477.1 | 10,723 bp |
| West Nile virus (NY99) | NC_009942.1 | 11,029 bp |
| West Nile virus (NY99 alt) | AY532665.1 | 11,029 bp |
| Powassan virus | HM440560.1 | 10,839 bp |
| Venezuelan equine encephalitis | NC_075022.1 | 11,522 bp |

## Troubleshooting

### Too many mutations in visualization
- Increase `--freq` in parsing step
- Increase `--cutoff` in visualization step
- Check for reference mismatch or contamination

### Unknown virus error
- Use enhanced pipeline: it automatically adds new viruses to known_viruses.json
- For manual addition, see "Adding New Viral Genomes to SnpEff Database" section
- Ensure exact accession format (e.g., NC_001477.1, not NC_001477)

### Generic gene names (Gene_106_10377)
- This indicates the virus needs BLAST-based annotation
- Follow the "Polyprotein-only annotation" section above
- The enhanced pipeline handles this automatically

### Memory errors with large files
- Use appropriate pipeline variant (standard/large/extreme)
- Increase thread count
- Consider splitting into smaller chunks

### Depth visualization shows log scale
- Add `--linear-scale` flag for actual depth values
- Log scale is default and recommended for high-depth samples

### SnpEff database errors
- "Database not found": Run SnpEff database update process
- "Most Exons do not have sequences": GFF3 format issue, regenerate using BLAST annotation
- "No CDS checked": ID mismatch between GFF3 and FASTA files, use update_snpeff_database_simple.py

### Polyprotein-only annotation (Gene_106_10377 errors)
If your variants are annotated with generic names like "Gene_106_10377" instead of real protein names, this means your virus genome only has polyprotein annotation in GenBank. You need to **choose a well-annotated reference genome** from the same virus family to extract individual protein coordinates.

**⚠️ IMPORTANT: Choose the Right Reference**
You must select a reference genome that:
1. **Same virus family** (Flavivirus, Alphavirus, etc.)
2. **Has individual mat_peptide annotations** (not just polyprotein)
3. **Is well-curated** (preferably NC_ RefSeq accessions)

```bash
# Step 1: Run BLAST annotation against a well-annotated RefSeq
# Replace REFSEQ_ACCESSION.1 with appropriate reference from table below
sbatch /scratch/sahlab/kathie/viral-genomics-pipeline/viral_pipeline/analysis/submit_blast_annotation.sh \
    YOUR_ACCESSION.1 \
    REFSEQ_ACCESSION.1

# Step 2: Update SnpEff database with new annotations
bash /scratch/sahlab/kathie/viral-genomics-pipeline/viral_pipeline/utils/update_snpeff_for_sample.sh \
    YOUR_ACCESSION.1

# Step 3: Re-run variant calling to get proper annotations
# (Re-run your original pipeline from Module 1)
```

**Reference Selection Guide:**

| Your Virus | Virus Family | Use This Reference | Example |
|------------|--------------|-------------------|----------|
| Any Zika virus | Flavivirus | **NC_012532.1** | `KU955591.1` → use `NC_012532.1` |
| Any Dengue virus | Flavivirus | **NC_001477.1** (DENV1) | `MK732862.1` → use `NC_001477.1` |
| Any West Nile virus | Flavivirus | **NC_009942.1** | `AY532665.1` → use `NC_009942.1` |
| Any Venezuelan EEV | Alphavirus | **NC_001449.1** | `KY000454.1` → use `NC_001449.1` |
| Any Chikungunya | Alphavirus | **NC_004162.2** | `MH823616.1` → use `NC_004162.2` |

This will replace generic gene IDs with descriptive names like:
- `anchored_capsid_protein_ancC`
- `envelope_protein_E`
- `nonstructural_protein_NS1`
- `RNA-dependent_RNA_polymerase_NS5`

## Adding New Viral Genomes to SnpEff Database

### When Do You Need to Update SnpEff?

You need to add a new genome to the SnpEff database when:

1. **New virus accession not in known_viruses.json**: First time using a specific viral accession
2. **Generic gene names in results**: Seeing annotations like "Gene_106_10377" instead of protein names
3. **Polyprotein-only viruses**: Viruses where GenBank only has polyprotein annotation, not individual proteins
4. **Missing gene coordinates**: Variants show as "intergenic" when they should be in genes

### How to Add New Genomes

**Method 1: Automatic (Enhanced Pipeline)**
```bash
# The enhanced pipeline automatically detects and adds new viruses
sbatch /scratch/sahlab/kathie/viral-genomics-pipeline/viral_pipeline/analysis/submit_viral_pipeline_enhanced.sh \
    "reads_R1.fastq.gz" \
    "reads_R2.fastq.gz" \
    NEW_ACCESSION.1 \
    4

# If virus not in known_viruses.json, it will:
# 1. Download GenBank file
# 2. Extract gene coordinates
# 3. Add to known_viruses.json
# 4. Create SnpEff database
# 5. Continue with variant calling
```

**Method 2: Manual BLAST Annotation (for polyprotein viruses)**
```bash
# Step 1: Run BLAST-based annotation
# CRITICAL: Choose the correct reference accession for your virus family
sbatch /scratch/sahlab/kathie/viral-genomics-pipeline/viral_pipeline/analysis/submit_blast_annotation.sh \
    YOUR_ACCESSION.1 \
    REFERENCE_ACCESSION.1

# Example for Zika virus:
sbatch /scratch/sahlab/kathie/viral-genomics-pipeline/viral_pipeline/analysis/submit_blast_annotation.sh \
    KU955591.1 \
    NC_012532.1

# Step 2: Update SnpEff database
bash /scratch/sahlab/kathie/viral-genomics-pipeline/viral_pipeline/utils/update_snpeff_for_sample.sh \
    YOUR_ACCESSION.1

# Step 3: Run/re-run variant calling
# Your samples will now use proper gene annotations
```

### BLAST Reference Selection Guide

| Virus Family | Best RefSeq References | Notes |
|--------------|----------------------|-------|
| **Flavivirus** | NC_012532.1 (ZIKV)<br>NC_001477.1 (DENV1)<br>NC_009942.1 (WNV) | Well-annotated mat_peptides |
| **Alphavirus** | NC_001449.1 (VEEV)<br>NC_001547.1 (SINV) | Complete polyprotein processing |
| **Coronavirus** | NC_045512.2 (SARS-CoV-2) | Extensive protein annotation |
| **Picornavirus** | NC_001612.1 (Poliovirus) | Classic polyprotein model |

### Verification Steps

After updating SnpEff database, verify success:

```bash
# Check that variants now show protein names
grep -v '^#' cleaned_seqs/variants/*.snpEFF.ann.tsv | head -10 | cut -f9,10,11

# Should see:
# missense_variant    MODERATE    capsid_protein_C
# missense_variant    MODERATE    envelope_protein_E
# Instead of:
# missense_variant    MODERATE    Gene_106_10377
```

### Troubleshooting SnpEff Database Issues

**Problem**: "Database 'ACCESSION.1' not found"
**Solution**: Run the SnpEff database update process above

**Problem**: Still seeing "Gene_106_10377" after update
**Solution**: 
1. Check that SnpEff database was properly created
2. Verify GFF3 file has proper CDS features with product attributes
3. Re-run variant calling (Module 1) to use updated database

**Problem**: BLAST annotation finds no proteins
**Solution**: 
1. Try a different well-annotated reference from the same virus family
2. Lower minimum coverage threshold: `--min-coverage 0.6`
3. Check that your genome file exists and is properly formatted