# Virus Annotation Methods for SnpEff Database

## Overview

There are two methods for annotating viral genomes for SnpEff database updates:

1. **GenBank Parser Method** - For viruses with complete mat_peptide annotations
2. **BLAST Annotation Method** - For viruses without mat_peptide annotations

## Method 1: GenBank Parser (Recommended when available)

Use this method when the GenBank record contains `mat_peptide` features with gene coordinates.

### When to use:
- GenBank record has `mat_peptide` annotations (e.g., DENV1, DENV2, DENV3, DENV4)
- You need precise gene boundaries from curated annotations
- Faster than BLAST (seconds vs minutes)

### How to use:
```bash
sbatch submit_genbank_annotation.sh NC_001477.1
```

### What it does:
1. Downloads GenBank file (if not present)
2. Extracts mat_peptide features with coordinates
3. Creates clean gene IDs (e.g., "nonstructural protein NS2A" → "NS2A")
4. Generates GFF file for SnpEff
5. Creates JSON files for visualization

### Example mat_peptide in GenBank:
```
mat_peptide     3476..4129
                /gene="POLY"
                /product="nonstructural protein NS2A"
                /protein_id="NP_733808.1"
```

### Output files:
- `NC_001477.1_genes.gff` - Gene annotations for SnpEff
- `NC_001477.1_annotation_summary.json` - Complete annotation data
- `NC_001477.1_gene_coords.json` - Gene coordinates for visualization

## Method 2: BLAST Annotation

Use this method when the GenBank record lacks mat_peptide annotations.

### When to use:
- No mat_peptide features in GenBank record
- Only CDS/gene features available (e.g., ZIKV, WNV, POWV)
- Need to infer gene boundaries from protein sequences

### How to use:
```bash
sbatch submit_blast_annotation.sh KU955591.1
```

### What it does:
1. Downloads GenBank file and extracts protein sequence
2. BLASTs against reference viruses with known annotations
3. Transfers gene coordinates based on alignment
4. Creates GFF file for SnpEff

### Output files:
- `KU955591.1_blast_genes.gff` - Gene annotations for SnpEff
- `KU955591.1_blast_annotation_summary.json` - BLAST results and annotations

## Decision Tree

```
Does the GenBank record have mat_peptide features?
├─ YES → Use GenBank Parser (submit_genbank_annotation.sh)
│   └─ Examples: DENV1-4, YFV, JEV
└─ NO → Use BLAST Annotation (submit_blast_annotation.sh)
    └─ Examples: ZIKV, WNV, POWV, VEEV
```

## Checking for mat_peptide features

To check if a GenBank record has mat_peptide annotations:

```python
from Bio import SeqIO
record = SeqIO.read("NC_001477.1.gb", "genbank")
mat_peptides = [f for f in record.features if f.type == "mat_peptide"]
print(f"Found {len(mat_peptides)} mat_peptide features")
```

## Common Viruses and Recommended Method

| Virus | Accession | Method | Reason |
|-------|-----------|---------|---------|
| DENV1 | NC_001477.1 | GenBank Parser | Has mat_peptide annotations |
| DENV2 | NC_001474.2 | GenBank Parser | Has mat_peptide annotations |
| DENV3 | NC_001475.2 | GenBank Parser | Has mat_peptide annotations |
| DENV4 | NC_002640.1 | GenBank Parser | Has mat_peptide annotations |
| ZIKV | KU955591.1 | BLAST | No mat_peptide annotations |
| WNV | NC_009942.1 | BLAST | No mat_peptide annotations |
| POWV | HM440560.1 | BLAST | No mat_peptide annotations |
| YFV | NC_002031.1 | GenBank Parser | Has mat_peptide annotations |
| JEV | NC_001437.1 | GenBank Parser | Has mat_peptide annotations |

## Notes

1. **GenBank Parser is preferred** when available because:
   - Uses curated gene boundaries
   - Faster execution
   - More accurate gene names
   - No alignment uncertainties

2. **BLAST method is necessary** for:
   - Newly sequenced genomes
   - Genomes with only CDS annotations
   - Custom virus isolates

3. **Both methods produce** compatible outputs for:
   - SnpEff database building
   - Pipeline visualization
   - Variant annotation