#!/usr/bin/env python3
"""
Module 9: Add Virus to Configuration
Automatically add a virus to known_viruses.json from its GenBank file
Can be run standalone or integrated into the pipeline
"""

import argparse
import json
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import re

def determine_virus_family(features, description):
    """Determine virus family based on gene names and description"""
    gene_names = []
    for feature in features:
        if feature.type in ['CDS', 'gene']:
            qualifiers = feature.qualifiers
            if 'gene' in qualifiers:
                gene_names.extend(qualifiers['gene'])
            elif 'product' in qualifiers:
                gene_names.extend(qualifiers['product'])
    
    gene_names_str = ' '.join(gene_names).lower()
    desc_lower = description.lower()
    
    # Check for flavivirus genes (NS1-NS5 pattern)
    if any(gene in gene_names_str for gene in ['ns1', 'ns2a', 'ns2b', 'ns3', 'ns4a', 'ns4b', 'ns5']):
        return 'flavivirus'
    
    # Check for alphavirus genes (nsP1-nsP4 pattern)
    if any(gene in gene_names_str for gene in ['nsp1', 'nsp2', 'nsp3', 'nsp4']):
        return 'alphavirus'
    
    # Check description
    if 'flavivirus' in desc_lower or 'dengue' in desc_lower or 'zika' in desc_lower or 'west nile' in desc_lower or 'yellow fever' in desc_lower:
        return 'flavivirus'
    if 'alphavirus' in desc_lower or 'chikungunya' in desc_lower or 'venezuelan' in desc_lower or 'semliki' in desc_lower:
        return 'alphavirus'
    
    # Default to flavivirus if uncertain
    print(f"Warning: Could not definitively determine virus family, defaulting to flavivirus")
    return 'flavivirus'

def get_gene_colors(family, genes):
    """Assign colors to genes based on virus family"""
    if family == 'flavivirus':
        structural_colors = {
            'C': '#4575b4',
            'prM': '#74add1', 
            'M': '#74add1',
            'Env': '#abd9e9',
            'E': '#abd9e9'
        }
        nonstructural_colors = {
            'NS1': '#fdae61',
            'NS2a': '#f46d43',
            'NS2A': '#f46d43',
            'NS2b': '#d73027',
            'NS2B': '#d73027',
            'NS3': '#a50026',
            'NS4a': '#762a83',
            'NS4A': '#762a83',
            'NS4b': '#9970ab',
            'NS4B': '#9970ab',
            'NS5': '#c2a5cf'
        }
    else:  # alphavirus
        structural_colors = {
            'C': '#4575b4',
            'E3': '#5ba5c9',
            'E2': '#7fb3d3',
            'E1': '#abd9e9'
        }
        nonstructural_colors = {
            'nsP1': '#fdae61',
            'nsP2': '#f46d43',
            'nsP3': '#d73027',
            'nsP4': '#a50026'
        }
    
    colors = {}
    structural_genes = []
    nonstructural_genes = []
    
    for gene in genes:
        if gene in structural_colors:
            colors[gene] = structural_colors[gene]
            structural_genes.append(gene)
        elif gene in nonstructural_colors:
            colors[gene] = nonstructural_colors[gene]
            nonstructural_genes.append(gene)
        else:
            # Default color for unknown genes
            colors[gene] = '#808080'
            if gene.startswith('NS') or gene.startswith('ns'):
                nonstructural_genes.append(gene)
            else:
                structural_genes.append(gene)
    
    return colors, structural_genes, nonstructural_genes

def extract_virus_info(genbank_file):
    """Extract virus information from GenBank file"""
    record = SeqIO.read(genbank_file, "genbank")
    
    virus_info = {
        'name': record.description,
        'genome_length': len(record.seq)
    }
    
    # Extract gene coordinates
    gene_coords = {}
    all_features = []
    
    for feature in record.features:
        if feature.type == 'CDS':
            # Get gene name
            gene_name = None
            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0]
            elif 'product' in feature.qualifiers:
                # Try to extract gene name from product
                product = feature.qualifiers['product'][0]
                # Handle various naming patterns
                if 'polyprotein' in product.lower():
                    continue  # Skip polyprotein entries
                elif 'protein' in product.lower():
                    # Extract gene names from product descriptions
                    parts = product.split()
                    for part in parts:
                        if part.upper() in ['C', 'E', 'M', 'E1', 'E2', 'E3'] or \
                           part.upper().startswith('NS') or part.lower().startswith('nsp'):
                            gene_name = part
                            break
                    
                    # Handle special cases
                    if not gene_name:
                        if 'capsid' in product.lower():
                            gene_name = 'C'
                        elif 'envelope' in product.lower() and 'E3' not in product:
                            gene_name = 'Env' if 'E2' not in product else 'E2'
                        elif 'membrane' in product.lower():
                            gene_name = 'prM' if 'precursor' in product.lower() else 'M'
                        elif re.search(r'NS\d[AB]?', product):
                            match = re.search(r'(NS\d[AB]?)', product)
                            gene_name = match.group(1)
                        elif re.search(r'nsP\d', product):
                            match = re.search(r'(nsP\d)', product)
                            gene_name = match.group(1)
                
            if gene_name:
                # Normalize gene names
                gene_name = gene_name.replace(' ', '')
                if gene_name.upper() in ['PRE-M', 'PREM', 'PRM']:
                    gene_name = 'prM'
                elif gene_name == 'E':
                    gene_name = 'Env'
                elif gene_name.upper().startswith('NS') and len(gene_name) > 2:
                    # Normalize NS genes (NS2a -> NS2a, NS2A -> NS2a)
                    gene_name = 'NS' + gene_name[2:]
                    if len(gene_name) == 4 and gene_name[3].isupper():
                        gene_name = gene_name[:3] + gene_name[3].lower()
                
                start = int(feature.location.start) + 1  # Convert to 1-based
                end = int(feature.location.end)
                
                # Only add if we don't already have this gene (avoid duplicates)
                if gene_name not in gene_coords:
                    gene_coords[gene_name] = [start, end]
                    all_features.append(feature)
    
    if not gene_coords:
        print("Warning: No gene coordinates found in GenBank file")
        print("The file may use polyprotein annotation. Manual configuration may be required.")
    
    # Determine virus family
    virus_info['family'] = determine_virus_family(all_features, record.description)
    
    # Get colors and gene classifications
    colors, structural, nonstructural = get_gene_colors(virus_info['family'], gene_coords.keys())
    
    virus_info['gene_coords'] = gene_coords
    virus_info['colors'] = colors
    virus_info['structural_genes'] = structural
    virus_info['nonstructural_genes'] = nonstructural
    
    return record.id, virus_info

def update_known_viruses(accession, virus_info, config_file):
    """Update known_viruses.json with new virus"""
    # Read existing config
    if config_file.exists():
        with open(config_file, 'r') as f:
            known_viruses = json.load(f)
    else:
        known_viruses = {}
    
    # Add or update virus
    known_viruses[accession] = virus_info
    
    # Write back with nice formatting
    with open(config_file, 'w') as f:
        json.dump(known_viruses, f, indent=2, sort_keys=True)
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Module 9: Add virus to configuration from GenBank file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add virus using GenBank file in current directory
  python3 add_virus_to_config.py KU955591.1.gb
  
  # Add virus from specific path
  python3 add_virus_to_config.py cleaned_seqs/variants/KU955591.1.gb
  
  # Specify custom config location
  python3 add_virus_to_config.py KU955591.1.gb --config /path/to/known_viruses.json
  
  # Force overwrite existing entry
  python3 add_virus_to_config.py KU955591.1.gb --force
        """
    )
    parser.add_argument('genbank_file', help='Path to GenBank file')
    parser.add_argument('--config', help='Path to known_viruses.json (will search standard locations if not specified)')
    parser.add_argument('--force', action='store_true', help='Overwrite existing virus entry')
    parser.add_argument('--pipeline-dir', help='Path to viral-genomics-pipeline directory', 
                       default='/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline')
    
    args = parser.parse_args()
    
    # Verify GenBank file exists
    gb_path = Path(args.genbank_file)
    if not gb_path.exists():
        print(f"❌ Error: GenBank file not found: {args.genbank_file}")
        sys.exit(1)
    
    # Find config file
    if args.config:
        config_file = Path(args.config)
    else:
        # Try to find it in standard locations
        possible_paths = [
            Path(args.pipeline_dir) / "viral_pipeline" / "visualization" / "known_viruses.json",
            Path(__file__).parent / "known_viruses.json",
            Path("known_viruses.json")
        ]
        config_file = None
        for path in possible_paths:
            if path.exists():
                config_file = path
                break
        
        if not config_file:
            # Use default location
            config_file = Path(args.pipeline_dir) / "viral_pipeline" / "visualization" / "known_viruses.json"
            config_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Module 9: Add Virus to Configuration")
    print("=" * 50)
    print(f"GenBank file: {gb_path}")
    print(f"Config file: {config_file}")
    
    # Extract virus info
    try:
        accession, virus_info = extract_virus_info(gb_path)
        print(f"\n📊 Extracted information for {accession}:")
        print(f"  Name: {virus_info['name']}")
        print(f"  Family: {virus_info['family']}")
        print(f"  Genome length: {virus_info['genome_length']:,} bp")
        print(f"  Genes found: {len(virus_info['gene_coords'])}")
        
        if virus_info['gene_coords']:
            print("\n  Gene coordinates:")
            for gene, coords in sorted(virus_info['gene_coords'].items()):
                print(f"    {gene}: {coords[0]:,}-{coords[1]:,}")
        
        # Check if already exists
        if config_file.exists():
            with open(config_file, 'r') as f:
                existing = json.load(f)
                if accession in existing and not args.force:
                    print(f"\n⚠️  {accession} already exists in config. Use --force to overwrite.")
                    sys.exit(1)
        
        # Update config
        if update_known_viruses(accession, virus_info, config_file):
            print(f"\n✅ Successfully added {accession} to {config_file}")
            
            # Show next steps
            print("\n📝 Next steps:")
            print(f"  1. Review the gene annotations in {config_file}")
            print(f"  2. Run Module 8 (consensus generation) to use these annotations")
            print(f"  3. If genes are missing, manually edit the config file")
        else:
            print(f"\n❌ Failed to update config file")
            sys.exit(1)
            
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()