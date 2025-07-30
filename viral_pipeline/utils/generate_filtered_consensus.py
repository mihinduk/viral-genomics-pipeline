#\!/usr/bin/env python3
"""
Generate consensus genome and proteins with multi-allelic site handling
Creates separate consensus sequences for each allele at multi-allelic positions
"""

import argparse
import sys; import os; sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "visualization")); from viral_translator import viral_translate
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

def read_virus_config(accession):
    """Read virus configuration from known_viruses.json"""
    possible_paths = [
        Path(__file__).parent / "known_viruses.json",
        Path(__file__).parent.parent / "visualization" / "known_viruses.json",
        Path(__file__).parent.parent.parent / "known_viruses.json",
        Path("known_viruses.json"),
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
    
    with open(vcf_file, 'r') as f:
        first_line = f.readline()
        if first_line.startswith('#CHROM'):
            df = pd.read_csv(vcf_file, sep='\t')
        else:
            df = pd.read_csv(vcf_file, sep='\t')
    
    print(f"Found {len(df)} total variants")
    
    filtered = df[
        (df['QUAL'] >= quality_cutoff) & 
        (df['Total_Depth'] >= depth_cutoff) &
        (df['Allele_Frequency'] >= freq_cutoff)
    ].copy()
    
    print(f"After filtering: {len(filtered)} variants pass criteria")
    return filtered

def identify_multiallelic_sites(mutations_df):
    """Identify positions with multiple alleles"""
    position_groups = mutations_df.groupby('POS')
    multiallelic = {}
    
    for pos, group in position_groups:
        if len(group) > 1:
            alleles = []
            for _, mut in group.iterrows():
                alleles.append({
                    'alt': mut['ALT'],
                    'freq': mut['Allele_Frequency'],
                    'mutation': mut
                })
            multiallelic[pos] = alleles
    
    return multiallelic

def apply_mutations_single_sequence(ref_seq, mutations_df, allele_id=""):
    """Apply mutations to single reference sequence"""
    seq_list = list(str(ref_seq))
    mutations_df = mutations_df.sort_values('POS')
    
    applied_count = 0
    for _, mut in mutations_df.iterrows():
        pos = int(mut['POS']) - 1
        ref_base = mut['REF']
        alt_base = mut['ALT']
        
        if pos < len(seq_list) and seq_list[pos] == ref_base:
            seq_list[pos] = alt_base
            applied_count += 1
        else:
            print(f"Warning: Reference mismatch at position {pos+1}: expected {ref_base}, found {seq_list[pos] if pos < len(seq_list) else 'out of range'}")
    
    allele_suffix = f" [{allele_id}]" if allele_id else ""
    print(f"Applied {applied_count} mutations to sequence{allele_suffix}")
    return ''.join(seq_list)

def generate_multiallelic_consensus(ref_seq, mutations_df):
    """Generate separate consensus sequences for each allele at multi-allelic sites"""
    multiallelic_sites = identify_multiallelic_sites(mutations_df)
    
    if not multiallelic_sites:
        # No multi-allelic sites, return single consensus
        consensus_seq = apply_mutations_single_sequence(ref_seq, mutations_df)
        return [{
            'sequence': consensus_seq,
            'mutations': mutations_df,
            'allele_id': '',
            'total_mutations': len(mutations_df)
        }]
    
    print(f"Found multi-allelic sites at positions: {list(multiallelic_sites.keys())}")
    
    # Generate separate consensus for each allele at multi-allelic sites
    consensus_results = []
    
    for pos, alleles in multiallelic_sites.items():
        for allele_info in alleles:
            # Create mutations dataframe for this specific allele
            other_mutations = mutations_df[mutations_df['POS'] != pos].copy()
            this_allele = pd.DataFrame([allele_info['mutation']])
            allele_mutations = pd.concat([other_mutations, this_allele], ignore_index=True)
            
            # Generate consensus for this allele combination
            mut_id = f"{pos}>{allele_info['alt']}"
            consensus_seq = apply_mutations_single_sequence(ref_seq, allele_mutations, mut_id)
            
            consensus_results.append({
                'sequence': consensus_seq,
                'mutations': allele_mutations,
                'allele_id': mut_id,
                'total_mutations': len(allele_mutations)
            })
    
    return consensus_results

def translate_genes_with_allele_info(sequence, virus_config, allele_id=""):
    """Translate sequence to proteins based on gene coordinates"""
    proteins = {}
    
    if not virus_config or 'gene_coords' not in virus_config:
        print("Warning: No gene coordinates available, translating full sequence")
        seq = Seq(sequence)
        protein = viral_translate(sequence, coordinates=None, stop_at_stop_codon=True)
        proteins['Polyprotein'] = str(protein)
        return proteins
    
    # Translate each gene
    for gene, (start, end) in virus_config['gene_coords'].items():
        gene_seq = sequence[start-1:end]
        protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
        proteins[gene] = {
            'sequence': protein,
            'allele_id': allele_id
        }
        
        print(f"Translated {gene}: {len(protein)} aa [{allele_id}]" if allele_id else f"Translated {gene}: {len(protein)} aa")
    
    return proteins

def main():
    parser = argparse.ArgumentParser(description='Generate consensus genome with multi-allelic handling')
    parser.add_argument('--vcf', required=True, help='Input VCF or filtered TSV file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--accession', required=True, help='Virus accession number')
    parser.add_argument('--quality', type=float, default=1000, help='Minimum quality score (default: 1000)')
    parser.add_argument('--depth', type=int, default=200, help='Minimum depth (default: 200)')
    parser.add_argument('--freq', type=float, default=0.01, help='Minimum allele frequency (default: 0.01)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--mutations-only', action='store_true', help='Only include positions with mutations in output')
    
    args = parser.parse_args()
    
    print("Multi-Allelic Consensus Generator")
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
        print("No mutations pass filtering criteria\!")
        sys.exit(0)
    
    # Generate consensus sequences (one per allele at multi-allelic sites)
    print(f"\nGenerating consensus sequences...")
    consensus_results = generate_multiallelic_consensus(ref_seq, filtered_mutations)
    
    # Get virus configuration
    virus_config = read_virus_config(args.accession)
    
    # Process each consensus sequence
    all_consensus_records = []
    all_protein_records = []
    
    for i, result in enumerate(consensus_results):
        sequence = result['sequence']
        allele_id = result['allele_id']
        total_muts = result['total_mutations']
        
        # Create consensus record
        suffix = f" [{allele_id}]" if allele_id else ""
        consensus_record = SeqRecord(
            Seq(sequence),
            id=f"{args.accession}_filtered_consensus",
            description=f"Consensus with {total_muts} mutations{suffix} (Q>={args.quality}, D>={args.depth}, F>={args.freq})"
        )
        all_consensus_records.append(consensus_record)
        
        # Translate to proteins
        print(f"\nTranslating to proteins{suffix}...")
        proteins = translate_genes_with_allele_info(sequence, virus_config, allele_id)
        
        # Create protein records
        for gene, protein_info in proteins.items():
            if isinstance(protein_info, dict):
                protein_seq = protein_info['sequence']
                gene_allele_id = protein_info['allele_id']
                gene_suffix = f" [{gene_allele_id}]" if gene_allele_id else ""
            else:
                protein_seq = protein_info
                gene_suffix = suffix
            
            record = SeqRecord(
                Seq(protein_seq),
                id=f"{args.accession}_{gene}",
                description=f"{gene} protein from filtered consensus{gene_suffix}"
            )
            all_protein_records.append(record)
    
    # Save all consensus genomes
    consensus_file = f"{args.output_prefix}_consensus.fasta"
    SeqIO.write(all_consensus_records, consensus_file, "fasta")
    print(f"\nWrote {len(all_consensus_records)} consensus sequences to: {consensus_file}")
    
    # Save all proteins
    protein_file = f"{args.output_prefix}_proteins.fasta"
    SeqIO.write(all_protein_records, protein_file, "fasta")
    print(f"Wrote {len(all_protein_records)} protein sequences to: {protein_file}")
    
    print(f"\n✅ Multi-allelic consensus generation complete\!")
    print(f"Generated {len(consensus_results)} consensus sequences from multi-allelic sites")

if __name__ == "__main__":
    main()
