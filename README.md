# Viral Genomics Pipeline v2.0

**A clean, modular, and extensible pipeline for viral genomics analysis**

## 🎯 Key Improvements in v2.0

- **Modular Architecture**: Clean separation of concerns
- **SnpEff Coordination**: Centralized gene naming and annotation management
- **Configuration-Driven**: YAML-based virus profiles with inheritance
- **No Hard-coding**: Dynamic virus family support
- **Clean Codebase**: Removed 93% redundancy from v1
- **Proper Versioning**: Semantic versioning and Git best practices

## 📁 Repository Structure

```
viral-genomics-pipeline-v2/
├── src/viral_pipeline/     # Core pipeline modules
│   ├── core/               # Core functionality
│   ├── analysis/           # Assembly and variant calling
│   ├── visualization/      # Plotting and reports
│   ├── utils/             # Utilities
│   └── config/            # Configuration management
├── config/                # Configuration files
│   ├── virus_profiles/    # Virus-specific configs
│   ├── family_templates/  # Family base templates
│   └── assembly_params/   # Assembly parameters
├── scripts/              # User-facing scripts
├── workflows/            # Snakemake/Nextflow workflows
├── tests/               # Test suite
├── docs/               # Documentation
└── examples/          # Example usage

```

## 🚀 Quick Start

### Installation

```bash
git clone https://github.com/mihinduk/viral-genomics-pipeline-v2.git
cd viral-genomics-pipeline-v2
conda env create -f environment.yml
conda activate viral_genomics
```

### Basic Usage

```bash
# Run complete pipeline
python scripts/run_viral_pipeline.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --accession NC_001477.1 \
  --threads 4 \
  --output results/
```

## 🔧 Pipeline Modules

1. **SnpEff Database Management**: Automatic database updates
2. **Reference Assembly**: BWA-based alignment
3. **Variant Calling**: LoFreq with quality filtering
4. **Mutation Parsing**: SnpEff annotation parsing
5. **Visualization**: Publication-ready plots
6. **Consensus Generation**: Filtered consensus sequences
7. **Quality Control**: Comprehensive QC reports

## 🧬 Supported Viruses

- **Flaviviruses**: DENV, ZIKA, WNV (all lineages), POWV
- **Alphaviruses**: VEEV, CHIKV
- **Extensible**: Easy addition of new virus families

## 📊 Configuration System

### Family Templates
- Base configuration for virus families
- Inherited by specific viruses
- Reduces redundancy

### Virus Profiles
- Specific virus configurations
- Override family defaults
- SnpEff gene mapping included

### Assembly Parameters
- Separate from virus configs
- Mix and match profiles
- Quality thresholds

## 🔄 Version History

- **v2.0.0** (2025-08): Complete restructure, modular design
- **v1.x.x**: Legacy version (deprecated)

## 📄 License

MIT License

## 🤝 Contributing

See [CONTRIBUTING.md](docs/CONTRIBUTING.md)

## 📧 Contact

- Issues: [GitHub Issues](https://github.com/mihinduk/viral-genomics-pipeline-v2/issues)
- Email: mihindu@wustl.edu
