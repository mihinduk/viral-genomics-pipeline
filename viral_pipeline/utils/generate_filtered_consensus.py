#!/usr/bin/env python3
"""
Generate consensus genome and proteins with only desired mutations
Filters mutations by quality, depth, and frequency before applying to reference
"""

import argparse
import sys; sys.path.append("viral_pipeline/visualization"); from viral_translator import viral_translate
import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

def read_virus_config(accession):
    """Read virus configuration from known_viruses.json"""
    # Try multiple possible locations for the config file
    possible_paths = [
        Path(__file__).parent / "known_viruses.json",  # Same directory as script
        Path(__file__).parent.parent / "visualization" / "known_viruses.json",  # Visualization directory
        Path(__file__).parent.parent.parent / "known_viruses.json",  # Project root
        Path("known_viruses.json"),  # Current working directory
    ]
    
    config_file = None
    for path in possible_paths:
        if path.exists():
            config_file = path
            break
    
    if not config_file:
        print(f"Warning: known_viruses.json not found in any expected location")
        return None
    
    with open(config_file, 'r') as f:
        viruses = json.load(f)
    
    return viruses.get(accession, None)

def filter_vcf(vcf_file, quality_cutoff, depth_cutoff, freq_cutoff):
    """Filter VCF/TSV file by quality, depth, and frequency"""
    print(f"Reading variants from {vcf_file}")
    
    # Check if it's a TSV or VCF
    with open(vcf_file, 'r') as f:
        first_line = f.readline()
        if first_line.startswith('#CHROM'):
            # It's a TSV
            df = pd.read_csv(vcf_file, sep='\t')
        else:
            # Assume it's our parsed TSV format
            df = pd.read_csv(vcf_file, sep='\t')
    
    print(f"Found {len(df)} total variants")
    
    # Apply filters
    filtered = df[
        (df['QUAL'] >= quality_cutoff) & 
        (df['Total_Depth'] >= depth_cutoff) &
        (df['Allele_Frequency'] >= freq_cutoff)
    ].copy()
    
    print(f"After filtering: {len(filtered)} variants pass criteria")
    return filtered

def apply_mutations_to_sequence(ref_seq, mutations_df):
    """Apply mutations to reference sequence"""
    # Convert to mutable list
    seq_list = list(str(ref_seq))
    
    # Sort mutations by position (important!)
    mutations_df = mutations_df.sort_values('POS')
    
    applied_count = 0
    for _, mut in mutations_df.iterrows():
        pos = int(mut['POS']) - 1  # Convert to 0-based
        ref_base = mut['REF']
        alt_base = mut['ALT']
        
        # Verify reference matches
        if pos < len(seq_list) and seq_list[pos] == ref_base:
            seq_list[pos] = alt_base
            applied_count += 1
        else:
            print(f"Warning: Reference mismatch at position {pos+1}: expected {ref_base}, found {seq_list[pos] if pos < len(seq_list) else 'out of range'}")
    
    print(f"Applied {applied_count} mutations to sequence")
    return ''.join(seq_list)

def translate_genes(sequence, virus_config):
    """Translate sequence to proteins based on gene coordinates"""
    proteins = {}
    
    if not virus_config or 'gene_coords' not in virus_config:
        print("Warning: No gene coordinates available, translating full sequence")
        # Translate as single polyprotein
        seq = Seq(sequence)
        protein = viral_translate(sequence, coordinates=None, stop_at_stop_codon=True)
        proteins['Polyprotein'] = str(protein)
        return proteins
    
    # Translate each gene
    for gene, (start, end) in virus_config['gene_coords'].items():
        # Extract gene sequence (coordinates are 1-based)
        gene_seq = sequence[start-1:end]
        
        # Use custom viral translator - NO start codon detection!
        protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
        proteins[gene] = protein
        
        print(f"Translated {gene}: {len(protein)} aa")
    
    return proteins

def main():
    parser = argparse.ArgumentParser(description='Generate consensus genome with filtered mutations')
    parser.add_argument('--vcf', required=True, help='Input VCF or filtered TSV file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--accession', required=True, help='Virus accession number')
    parser.add_argument('--quality', type=float, default=1000, help='Minimum quality score (default: 1000)')
    parser.add_argument('--depth', type=int, default=200, help='Minimum depth (default: 200)')
    parser.add_argument('--freq', type=float, default=0.05, help='Minimum allele frequency (default: 0.05)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--mutations-only', action='store_true', help='Only include positions with mutations in output')
    
    args = parser.parse_args()
    
    print("Consensus Genome Generator")
    print("=" * 50)
    print(f"Reference: {args.reference}")
    print(f"VCF/TSV: {args.vcf}")
    print(f"Quality cutoff: {args.quality}")
    print(f"Depth cutoff: {args.depth}")
    print(f"Frequency cutoff: {args.freq}")
    
    # Read reference sequence
    print(f"\nReading reference genome...")
    ref_record = SeqIO.read(args.reference, "fasta")
    ref_seq = str(ref_record.seq).upper()
    print(f"Reference length: {len(ref_seq)} bp")
    
    # Filter mutations
    filtered_mutations = filter_vcf(args.vcf, args.quality, args.depth, args.freq)
    
    if len(filtered_mutations) == 0:
        print("No mutations pass filtering criteria!")
        sys.exit(0)
    
    # Apply mutations
    print(f"\nApplying mutations to reference...")
    consensus_seq = apply_mutations_to_sequence(ref_seq, filtered_mutations)
    
    # Save consensus genome
    consensus_file = f"{args.output_prefix}_consensus.fasta"
    consensus_record = SeqRecord(
        Seq(consensus_seq),
        id=f"{args.accession}_filtered_consensus",
        description=f"Consensus with {len(filtered_mutations)} mutations (Q>={args.quality}, D>={args.depth}, F>={args.freq})"
    )
    SeqIO.write(consensus_record, consensus_file, "fasta")
    print(f"\nWrote consensus genome to: {consensus_file}")
    
    # Get virus configuration
    virus_config = read_virus_config(args.accession)
    
    # Translate to proteins
    print(f"\nTranslating to proteins...")
    proteins = translate_genes(consensus_seq, virus_config)
    
    # Save proteins
    protein_file = f"{args.output_prefix}_proteins.fasta"
    protein_records = []
    for gene, protein_seq in proteins.items():
        record = SeqRecord(
            Seq(protein_seq),
            id=f"{args.accession}_{gene}",
            description=f"{gene} protein from filtered consensus"
        )
        protein_records.append(record)
    
    SeqIO.write(protein_records, protein_file, "fasta")
    print(f"Wrote {len(proteins)} proteins to: {protein_file}")
    
    # Save mutation summary
    summary_file = f"{args.output_prefix}_mutations_applied.tsv"
    summary_df = filtered_mutations[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'Total_Depth', 'Allele_Frequency', 'EFFECT', 'HGVSp']].copy()
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Wrote mutation summary to: {summary_file}")
    
    print("\n✅ Consensus generation complete!")

if __name__ == "__main__":
    main()