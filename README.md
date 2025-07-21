# Viral Genomics Pipeline

**End-to-end viral genomics pipeline: from raw sequencing reads to publication-ready mutation visualizations.**

The only integrated analysis-to-visualization toolkit for defensive viral surveillance.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)

## 🌟 Features

- **Complete workflow**: Raw reads → Assembly → Variant calling → Visualization
- **Family-based configuration**: Dynamic gene layouts for different virus families
- **Publication-quality plots**: Beautiful mutation visualizations with meaningful colors
- **Interactive reports**: HTML depth coverage with gene annotations
- **Dual filtering system**: Quality filtering + frequency filtering for clean results
- **Multi-virus support**: Flaviviruses (DENV, WNV, POWV), Alphaviruses (VEEV), and more

## 🚀 Quick Start

### Prerequisites

```bash
# Conda/Mamba environment with viral genomics tools
mamba create -n viral_genomics bwa samtools lofreq snpeff pandas matplotlib numpy
conda activate viral_genomics
```

### Installation

```bash
git clone https://github.com/mihinduk/viral-genomics-pipeline.git
cd viral-genomics-pipeline
chmod +x scripts/run_viral_pipeline.py
```

### Basic Usage

```bash
# Complete pipeline for DENV1
python scripts/run_viral_pipeline.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --accession NC_001477.1 \
  --threads 4 \
  --quality-cutoff 1000 \
  --freq-cutoff 0.01 \
  --output sample_results
```

## 📋 Step-by-Step Examples

### Example 1: Dengue Virus Type 1 (DENV1)

```bash
# 1. Assembly and Analysis
/scratch/sahlab/kathie/Diamond_test/shotgun_viral_genomics/run_pipeline_htcf_consolidated.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  NC_001477.1 \
  4

# 2. Create Output Directory
mkdir -p sample_results

# 3. Generate Depth File
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  samtools depth cleaned_seqs/mapping/*.lofreq.final.bam > sample_results/sample_depth.txt

# 4. Quality Filter VCF
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 viral_pipeline/visualization/parse_snpeff_vcf.py \
  -i cleaned_seqs/variants/*NC_001477*.snpEFF.ann.vcf \
  -d 200 \
  -q 1000 \
  -o sample_filtered_mutations.tsv \
  -O sample_results

# 5. Create Mutation Visualization
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input sample_results/sample_filtered_mutations.tsv \
  --output sample_mutations.png \
  --accession NC_001477.1 \
  --cutoff 0.01 \
  --outdir sample_results

# 6. Create Depth Visualization
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 viral_pipeline/visualization/visualize_depth.py \
  --depth sample_results/sample_depth.txt \
  --output sample_results/sample_depth.png \
  --output-html sample_results/sample_depth.html \
  --accession NC_001477.1

# 7. Generate Diagnostics Report
sbatch viral_pipeline/analysis/submit_viral_diagnostic.sh \
  "./sample_R1.fastq.gz" \
  "./sample_R2.fastq.gz" \
  NC_001477.1 \
  sample \
  4
```

### Example 2: Venezuelan Equine Encephalitis Virus (VEEV)

```bash
# Same workflow, different accession
python scripts/run_viral_pipeline.py \
  --r1 veev_R1.fastq.gz \
  --r2 veev_R2.fastq.gz \
  --accession NC_075022.1 \
  --quality-cutoff 49314 \
  --freq-cutoff 0.001 \
  --output veev_sample
```

### Example 3: Visualization Only (Existing Data)

```bash
# If you already have VCF and BAM files
python scripts/run_viral_pipeline.py \
  --visualization-only \
  --vcf sample.snpEFF.ann.vcf \
  --bam sample.lofreq.final.bam \
  --accession NC_001477.1 \
  --quality-cutoff 1000 \
  --freq-cutoff 0.01 \
  --output existing_sample
```

## 🎨 Visualization Features

### Mutation Plots
- **Gene-specific colors**: Different colors for each viral gene
- **Mutation type coloring**:
  - 🟢 Green: Synonymous variants (silent)
  - 🟠 Orange: Missense variants (amino acid changes)
  - 🔴 Red: Nonsense/Stop variants (severe)
- **Complete parity**: Every mutation line has corresponding table entry
- **Family-specific layouts**: Flavivirus vs Alphavirus gene arrangements

### Depth Coverage
- **Interactive HTML reports**: Zoomable, gene-annotated coverage plots
- **Coverage statistics**: Mean depth, coverage percentages
- **Gene-colored regions**: Visual gene boundaries with family colors
- **Quality thresholds**: Configurable minimum depth highlighting

## 🧬 Supported Viruses

| Virus | Accession | Family | Status |
|-------|-----------|---------|---------|
| Dengue virus type 1 | NC_001477.1 | Flavivirus | ✅ Tested |
| West Nile virus | NC_009942.1 | Flavivirus | ✅ Tested |
| Powassan virus | HM440560.1 | Flavivirus | ✅ Tested |
| Venezuelan equine encephalitis virus | NC_075022.1 | Alphavirus | ✅ Tested |

*New viruses are automatically configured based on GenBank annotations and family classification.*

## ⚙️ Configuration

### Quality Cutoffs

The pipeline uses dual filtering for optimal results:

1. **Quality Cutoff** (`--quality-cutoff`): Filters mutations by sequencing quality
   - DENV1: `1000` (stringent for clean samples)
   - VEEV: `49314` (very stringent for high-quality data)
   - Inspect your VCF QUAL column to choose appropriate values

2. **Frequency Cutoff** (`--freq-cutoff`): Filters mutations for visualization
   - `0.01` (1%): Show common variants
   - `0.001` (0.1%): Show rare variants
   - `0.05` (5%): Show only major variants

### Adding New Viruses

Add virus configurations to `viral_pipeline/config/known_viruses.json`:

```json
{
  "NC_XXXXXX.1": {
    "name": "Virus name",
    "family": "flavivirus",
    "genome_length": 10000,
    "gene_coords": {
      "gene1": [start, end],
      "gene2": [start, end]
    },
    "colors": {
      "gene1": "#color1",
      "gene2": "#color2"
    },
    "structural_genes": ["gene1"],
    "nonstructural_genes": ["gene2"]
  }
}
```

## 📊 Output Files

```
sample_results/
├── sample_mutations.png              # Publication-ready mutation plot
├── sample_mutations_mutations_table.tsv  # Complete mutation data
├── sample_depth.png                  # Depth coverage plot
├── sample_depth.html                 # Interactive depth visualization
├── sample_filtered_mutations.tsv     # Quality-filtered mutations
├── sample_depth.txt                  # Raw depth data
└── diagnostic_sample_report.html     # Comprehensive diagnostic report
```

## 🛠️ Pipeline Components

### Analysis Module (`viral_pipeline/analysis/`)
- **Assembly pipeline**: BWA alignment, LoFreq variant calling
- **Annotation pipeline**: SnpEff effect prediction
- **Quality control**: Coverage analysis, contamination detection

### Visualization Module (`viral_pipeline/visualization/`)
- **Mutation visualizer**: Dynamic gene layouts, publication-quality plots
- **Depth visualizer**: Coverage analysis with gene annotations
- **Configuration system**: Family-based virus configurations

### Utilities (`viral_pipeline/utils/`)
- **VCF parsers**: SnpEff output processing
- **Configuration managers**: Dynamic virus detection
- **Report generators**: HTML and table outputs

## 🔬 Scientific Applications

### Defensive Viral Surveillance
- **Outbreak monitoring**: Track viral genetic changes
- **Vaccine strain selection**: Identify antigenic variants
- **Drug resistance**: Monitor resistance mutations

### Research Applications
- **Comparative genomics**: Cross-family viral analysis
- **Evolution studies**: Mutation rate analysis
- **Diagnostic development**: Variant characterization

## 📚 Citation

If you use this pipeline in your research, please cite:

```
Viral Genomics Pipeline: End-to-end analysis and visualization toolkit
for defensive viral surveillance. 
GitHub: https://github.com/mihinduk/viral-genomics-pipeline
```

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## 📧 Support

- **Issues**: [GitHub Issues](https://github.com/mihinduk/viral-genomics-pipeline/issues)
- **Documentation**: [Wiki](https://github.com/mihinduk/viral-genomics-pipeline/wiki)
- **Discussions**: [GitHub Discussions](https://github.com/mihinduk/viral-genomics-pipeline/discussions)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Shotgun Viral Genomics Pipeline**: Foundation analysis workflow
- **Viral Mutation Visualizer**: Core visualization components
- **Claude Code**: Development assistance and optimization

---

**🦠 Making viral genomics accessible to everyone - from raw reads to publication-ready insights.**