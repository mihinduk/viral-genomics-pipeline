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
- Ensure accession is in `${PIPELINE_DIR}/viral_pipeline/visualization/known_viruses.json`
- Use exact accession format (e.g., NC_001477.1, not NC_001477)

### Memory errors with large files
- Use appropriate pipeline variant (standard/large/extreme)
- Increase thread count
- Consider splitting into smaller chunks

### Depth visualization shows log scale
- Add `--linear-scale` flag for actual depth values
- Log scale is default and recommended for high-depth samples