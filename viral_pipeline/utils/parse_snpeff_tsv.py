#!/usr/bin/env python3
import pandas as pd
import re
import sys

if len(sys.argv) != 3:
    print("Usage: python parse_snpeff_tsv.py input.tsv output.tsv")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the SnpEff TSV file (skip lines starting with ## but keep #CHROM header)
lines = []
with open(input_file, 'r') as f:
    for line in f:
        if not line.startswith('##'):
            lines.append(line)

import io
df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

# Extract AF from INFO field
def extract_af(info):
    match = re.search(r'AF=([0-9.]+)', str(info))
    return float(match.group(1)) if match else 0.0

# Extract DP from INFO field  
def extract_dp(info):
    match = re.search(r'DP=([0-9]+)', str(info))
    return int(match.group(1)) if match else 0

# Parse the data
df['Allele_Frequency'] = df['INFO'].apply(extract_af)
df['Total_Depth'] = df['INFO'].apply(extract_dp)

# Create output format matching working VEEV file
output_df = pd.DataFrame({
    'CHROM': df['#CHROM'],
    'POS': df['POS'],
    'ID': df['ID'],
    'REF': df['REF'],
    'ALT': df['ALT'],
    'QUAL': df['QUAL'],
    'INFO': df['INFO'],
    'Total_Depth': df['Total_Depth'],
    'Allele_Frequency': df['Allele_Frequency'],
    'strand_bias': '0',  # placeholder
    'DP4': '0,0,0,0',  # placeholder
    'EFFECT': df['EFFECT'],
    'PUTATIVE_IMPACT': df['PUTATIVE_IMPACT'],
    'GENE_NAME': df['GENE_NAME'],
    'GENE_ID': df['GENE_ID'],
    'FEATURE_TYPE': df['FEATURE_TYPE'],
    'FEATURE_ID': df['FEATURE_ID'],
    'TRANSCRIPT_TYPE': df['TRANSCRIPT_TYPE'],
    'HGVSc': df['HGVSc'],
    'HGVSp': df['HGVSp'],
    'cDNA_POSITION_AND_LENGTH': df['cDNA_POSITION_AND_LENGTH'],
    'CDS_POSITION_AND_LENGTH': df['CDS_POSITION_AND_LENGTH'],
    'PROTEIN_POSITION_AND_LENGTH': df['PROTEIN_POSITION_AND_LENGTH'],
    'ERROR': df['ERROR']
})

# Save the parsed file
output_df.to_csv(output_file, sep='\t', index=False)
print(f"Parsed {len(output_df)} mutations to {output_file}")