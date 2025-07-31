#!/usr/bin/env python3
"""
Generate consensus genome and proteins with intelligent multi-allelic site handling
Creates separate sequences only when amino acid changes differ
Generates comprehensive summary report
"""

import argparse
import sys; import os; sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "visualization")); from viral_translator import viral_translate
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from collections import defaultdict

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
                    'mutation': mut,
                    'effect': mut.get('EFFECT', ''),
                    'hgvsp': mut.get('HGVSp', '')
                })
            multiallelic[pos] = alleles
    
    return multiallelic

def find_gene_for_position(pos, gene_coords):
    """Find which gene contains a given position"""
    for gene, (start, end) in gene_coords.items():
        if start <= pos <= end:
            return gene, start, end
    return None, None, None

def parse_aa_change(hgvsp):
    """Parse amino acid change from HGVSp annotation"""
    if not hgvsp or pd.isna(hgvsp):
        return None, None, None
    
    # Handle various formats: p.Leu287Ile, p.L287I, etc
    import re
    match = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', hgvsp)
    if match:
        ref_aa = match.group(1)
        position = int(match.group(2))
        alt_aa = match.group(3)
        
        # Convert 3-letter to 1-letter if needed
        aa_3to1 = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
            'Ter': '*', 'TER': '*'
        }
        
        if len(ref_aa) > 1:
            ref_aa = aa_3to1.get(ref_aa, ref_aa)
        if len(alt_aa) > 1:
            alt_aa = aa_3to1.get(alt_aa, alt_aa)
            
        return ref_aa, position, alt_aa
    
    return None, None, None

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

def generate_multiallelic_consensus(ref_seq, mutations_df, virus_config):
    """Generate consensus sequences with intelligent multi-allelic handling"""
    multiallelic_sites = identify_multiallelic_sites(mutations_df)
    
    if not multiallelic_sites:
        # No multi-allelic sites, return single consensus
        consensus_seq = apply_mutations_single_sequence(ref_seq, mutations_df)
        return [{
            'sequence': consensus_seq,
            'mutations': mutations_df,
            'allele_id': '',
            'total_mutations': len(mutations_df),
            'multiallelic_info': {}
        }]
    
    print(f"Found multi-allelic sites at positions: {list(multiallelic_sites.keys())}")
    
    # Analyze which genes are affected
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    multiallelic_analysis = {}
    
    for pos, alleles in multiallelic_sites.items():
        gene, start, end = find_gene_for_position(pos, gene_coords)
        aa_changes = []
        
        for allele in alleles:
            ref_aa, aa_pos, alt_aa = parse_aa_change(allele['hgvsp'])
            aa_change = f"{ref_aa}{aa_pos}{alt_aa}" if ref_aa and alt_aa else "Unknown"
            aa_changes.append({
                'nucleotide': f"{pos}>{allele['alt']}",
                'aa_change': aa_change,
                'effect': allele['effect'],
                'freq': allele['freq']
            })
        
        multiallelic_analysis[pos] = {
            'gene': gene,
            'aa_changes': aa_changes,
            'unique_aa_changes': len(set(ch['aa_change'] for ch in aa_changes))
        }
    
    # Generate consensus sequences
    consensus_results = []
    
    for pos, alleles in multiallelic_sites.items():
        for allele_info in alleles:
            # Create mutations dataframe for this specific allele
            other_mutations = mutations_df[mutations_df['POS'] != pos].copy()
            this_allele = pd.DataFrame([allele_info['mutation']])
            allele_mutations = pd.concat([other_mutations, this_allele], ignore_index=True)
            
            # Generate consensus for this allele combination
            # Get reference base
            ref_base = allele_info["mutation"]["REF"]
            mut_id = f"{pos}{ref_base}>{allele_info['alt']}"
            consensus_seq = apply_mutations_single_sequence(ref_seq, allele_mutations, mut_id)
            
            consensus_results.append({
                'sequence': consensus_seq,
                'mutations': allele_mutations,
                'allele_id': mut_id,
                'total_mutations': len(allele_mutations),
                'multiallelic_info': multiallelic_analysis
            })
    
    return consensus_results

def translate_and_deduplicate_proteins(consensus_results, virus_config, summary_data):
    """Translate proteins and deduplicate identical sequences"""
    all_proteins = {}
    protein_mutation_summary = defaultdict(lambda: defaultdict(list))
    
    for result in consensus_results:
        sequence = result['sequence']
        allele_id = result['allele_id']
        mutations = result['mutations']
        
        # Get gene coordinates
        gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
        
        if not gene_coords:
            # Translate as single polyprotein
            protein = viral_translate(sequence, coordinates=None, stop_at_stop_codon=True)
            key = ('Polyprotein', protein)
            if key not in all_proteins:
                all_proteins[key] = {
                    'sequence': protein,
                    'allele_ids': [allele_id] if allele_id else [],
                    'mutations': []
                }
            else:
                if allele_id:
                    all_proteins[key]['allele_ids'].append(allele_id)
        else:
            # Translate each gene
            for gene, (start, end) in gene_coords.items():
                gene_seq = sequence[start-1:end]
                protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
                
                # Find mutations in this gene
                gene_mutations = mutations[(mutations['POS'] >= start) & (mutations['POS'] <= end)]
                
                # Parse AA changes for this gene
                aa_changes = []
                for _, mut in gene_mutations.iterrows():
                    ref_aa, aa_pos, alt_aa = parse_aa_change(mut.get('HGVSp', ''))
                    if ref_aa and alt_aa:
                        aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
                        aa_changes.append(aa_change)
                        protein_mutation_summary[gene][aa_change].append({
                            'nucleotide': f"{mut['POS']}>{mut['ALT']}",
                            'effect': mut.get('EFFECT', 'unknown'),
                            'freq': mut.get('Allele_Frequency', 0)
                        })
                
                # Create unique key for this protein
                key = (gene, protein, tuple(sorted(aa_changes)))
                
                if key not in all_proteins:
                    all_proteins[key] = {
                        'sequence': protein,
                        'allele_ids': [allele_id] if allele_id else [],
                        'mutations': aa_changes,
                        'gene': gene
                    }
                else:
                    if allele_id:
                        all_proteins[key]['allele_ids'].append(allele_id)
    
    # Update summary data
    summary_data['protein_mutation_summary'] = dict(protein_mutation_summary)
    
    return all_proteins

def generate_summary_report(summary_data, output_prefix):
    """Generate comprehensive summary report"""
    report_file = f"{output_prefix}_summary_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MULTI-ALLELIC CONSENSUS GENERATION SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        # Basic statistics
        f.write(f"Total mutations passing filters: {summary_data['total_mutations']}\n")
        f.write(f"Multi-allelic sites found: {len(summary_data['multiallelic_sites'])}\n")
        if summary_data['multiallelic_sites']:
            f.write(f"Multi-allelic positions: {summary_data['multiallelic_sites']}\n")
        f.write(f"Consensus sequences generated: {summary_data['consensus_count']}\n")
        f.write(f"Unique protein sequences: {summary_data['unique_proteins']}\n\n")
        
        # Multi-allelic site details
        if summary_data['multiallelic_analysis']:
            f.write("-" * 80 + "\n")
            f.write("MULTI-ALLELIC SITE ANALYSIS\n")
            f.write("-" * 80 + "\n\n")
            
            for pos, info in summary_data['multiallelic_analysis'].items():
                f.write(f"Position {pos} (Gene: {info['gene'] or 'Unknown'}):\n")
                for change in info['aa_changes']:
                    f.write(f"  - {change['nucleotide']} → {change['aa_change']} ")
                    f.write(f"({change['effect']}, freq: {change['freq']:.2%})\n")
                f.write(f"  Unique amino acid changes: {info['unique_aa_changes']}\n\n")
        
        # Protein mutation summary
        if summary_data['protein_mutation_summary']:
            f.write("-" * 80 + "\n")
            f.write("MUTATIONS PER PROTEIN\n")
            f.write("-" * 80 + "\n\n")
            
            for gene, mutations in summary_data['protein_mutation_summary'].items():
                f.write(f"{gene}:\n")
                if mutations:
                    for aa_change, occurrences in mutations.items():
                        f.write(f"  - {aa_change}:")
                        for occ in occurrences:
                            f.write(f" {occ['nucleotide']} ({occ['effect']}, {occ['freq']:.2%})")
                        f.write("\n")
                else:
                    f.write("  No mutations\n")
                f.write("\n")
        
        # Protein deduplication summary
        if summary_data.get('protein_dedup_info'):
            f.write("-" * 80 + "\n")
            f.write("PROTEIN DEDUPLICATION SUMMARY\n")
            f.write("-" * 80 + "\n\n")
            
            for info in summary_data['protein_dedup_info']:
                f.write(f"{info['gene']}: ")
                if len(info['allele_ids']) > 1:
                    f.write(f"Identical protein from {len(info['allele_ids'])} alleles: ")
                    f.write(f"Unique protein from {info['allele_ids'][0]}\n")
                if info['mutations']:
                    f.write(f"  Mutations: {', '.join(info['mutations'])}\n")
                f.write("\n")
    
    print(f"\n📄 Summary report saved to: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate consensus genome with intelligent multi-allelic handling')
    parser.add_argument('--vcf', required=True, help='Input VCF or filtered TSV file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--accession', required=True, help='Virus accession number')
    parser.add_argument('--quality', type=float, default=1000, help='Minimum quality score (default: 1000)')
    parser.add_argument('--depth', type=int, default=200, help='Minimum depth (default: 200)')
    parser.add_argument('--freq', type=float, default=0.01, help='Minimum allele frequency (default: 0.01)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument("--sample-name", help="Sample name for sequence IDs (default: extracted from output-prefix)")
    parser.add_argument('--mutations-only', action='store_true', help='Only include positions with mutations in output')
    
    args = parser.parse_args()

    # Determine sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        # Extract from output prefix
        base_name = os.path.basename(args.output_prefix)
        parts = base_name.split("_")
        if len(parts) >= 2:
            sample_name = "_".join(parts[:2])
        else:
            sample_name = parts[0]
    
    print(f"Using sample name: {sample_name}")
    
    print("Multi-Allelic Consensus Generator v2.0")
    print("=" * 50)
    print(f"Reference: {args.reference}")
    print(f"VCF/TSV: {args.vcf}")
    print(f"Quality cutoff: {args.quality}")
    print(f"Depth cutoff: {args.depth}")
    print(f"Frequency cutoff: {args.freq}")
    
    # Initialize summary data
    summary_data = {
        'total_mutations': 0,
        'multiallelic_sites': [],
        'multiallelic_analysis': {},
        'consensus_count': 0,
        'unique_proteins': 0,
        'protein_mutation_summary': {},
        'protein_dedup_info': []
    }
    
    # Read reference sequence
    print(f"\nReading reference genome...")
    ref_record = SeqIO.read(args.reference, "fasta")
    ref_seq = str(ref_record.seq).upper()
    print(f"Reference length: {len(ref_seq)} bp")
    
    # Filter mutations
    filtered_mutations = filter_vcf(args.vcf, args.quality, args.depth, args.freq)
    summary_data['total_mutations'] = len(filtered_mutations)
    
    if len(filtered_mutations) == 0:
        print("No mutations pass filtering criteria!")
        sys.exit(0)
    
    # Get virus configuration
    virus_config = read_virus_config(args.accession)
    
    # Generate consensus sequences
    print(f"\nGenerating consensus sequences...")
    consensus_results = generate_multiallelic_consensus(ref_seq, filtered_mutations, virus_config)
    summary_data['consensus_count'] = len(consensus_results)
    
    # Extract multi-allelic info for summary
    if consensus_results[0]['multiallelic_info']:
        summary_data['multiallelic_sites'] = list(consensus_results[0]['multiallelic_info'].keys())
        summary_data['multiallelic_analysis'] = consensus_results[0]['multiallelic_info']
    
    # Save consensus genomes
    consensus_records = []
    for result in consensus_results:
        suffix = f" [{result['allele_id']}]" if result['allele_id'] else ""
        consensus_record = SeqRecord(
            Seq(result['sequence']),
            id=f"{sample_name}_filtered_consensus",
            description=f"Consensus with {result['total_mutations']} mutations{suffix} (Q>={args.quality}, D>={args.depth}, F>={args.freq})"
        )
        consensus_records.append(consensus_record)
    
    consensus_file = f"{args.output_prefix}_consensus.fasta"
    SeqIO.write(consensus_records, consensus_file, "fasta")
    print(f"\nWrote {len(consensus_records)} consensus sequences to: {consensus_file}")
    
    # Translate and deduplicate proteins
    print(f"\nTranslating and deduplicating proteins...")
    unique_proteins = translate_and_deduplicate_proteins(consensus_results, virus_config, summary_data)
    summary_data['unique_proteins'] = len(unique_proteins)
    
    # Save proteins
    protein_records = []
    for key, protein_data in unique_proteins.items():
        # Handle both polyprotein (2 elements) and gene-specific (3 elements) keys
        if len(key) == 2:
            # Polyprotein case: (gene_name, protein_seq)
            gene_info, protein_seq = key
            mutations = tuple()  # No mutations for polyprotein
        elif len(key) == 3:
            # Gene-specific case: (gene_name, protein_seq, mutations_tuple)
            gene_info, protein_seq, mutations = key
        else:
            print(f"Warning: Unexpected key format: {key}")
            continue
        gene = protein_data.get('gene', gene_info)
        
        # Create description
        if protein_data['mutations']:
            mutation_str = ", ".join(protein_data['mutations'])
            desc = f"{gene} protein [{mutation_str}]"
        else:
            desc = f"{gene} protein"
        
        # Add allele info if multiple alleles produce same protein
        if len(protein_data['allele_ids']) > 1:
            desc += f" (quasispecies: {', '.join(protein_data['allele_ids'])})"
        elif protein_data['allele_ids']:
            desc += f" (from {protein_data['allele_ids'][0]})"
        
        record = SeqRecord(
            Seq(protein_data['sequence']),
            id=f"{sample_name}_{gene}",
            description=desc
        )
        protein_records.append(record)
        
        # Add to dedup info
        summary_data['protein_dedup_info'].append({
            'gene': gene,
            'allele_ids': protein_data['allele_ids'],
            'mutations': protein_data['mutations']
        })
    
    protein_file = f"{args.output_prefix}_proteins.fasta"
    SeqIO.write(protein_records, protein_file, "fasta")
    print(f"Wrote {len(protein_records)} unique protein sequences to: {protein_file}")
    
    # Generate summary report
    generate_summary_report(summary_data, args.output_prefix)
    
    print(f"\n✅ Multi-allelic consensus generation complete!")
    print(f"Generated {len(consensus_results)} consensus sequences")
    print(f"Generated {len(unique_proteins)} unique proteins")

if __name__ == "__main__":
    main()
