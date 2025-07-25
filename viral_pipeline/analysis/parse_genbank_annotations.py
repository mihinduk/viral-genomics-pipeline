#!/usr/bin/env python3
"""
Parse GenBank mat_peptide annotations for SnpEff database update
Extracts gene coordinates and names directly from GenBank records
No BLAST required - uses existing annotations in the record
"""

import argparse
import json
import os
import sys
from Bio import SeqIO
from pathlib import Path

def parse_mat_peptides(genbank_file):
    """Extract mat_peptide features from GenBank file"""
    
    gene_annotations = {}
    
    with open(genbank_file, 'r') as f:
        record = SeqIO.read(f, 'genbank')
        
    print(f"Parsing {record.id} - {record.description}")
    print(f"Sequence length: {len(record.seq)} bp")
    
    # Extract mat_peptide features
    mat_peptides = [f for f in record.features if f.type == 'mat_peptide']
    
    if not mat_peptides:
        print("Warning: No mat_peptide features found in GenBank record")
        print("Available feature types:", set(f.type for f in record.features))
        return None
    
    print(f"\nFound {len(mat_peptides)} mat_peptide features:")
    
    for i, feature in enumerate(mat_peptides, 1):
        # Get coordinates
        start = int(feature.location.start) + 1  # Convert to 1-based
        end = int(feature.location.end)
        
        # Get gene name from product field
        product = feature.qualifiers.get('product', ['Unknown'])[0]
        
        # Create a clean gene ID from the product name
        # e.g., "nonstructural protein NS2A" -> "NS2A"
        gene_id = create_gene_id(product)
        
        # Store gene info
        gene_annotations[gene_id] = {
            'start': start,
            'end': end,
            'product': product,
            'type': classify_gene_type(product),
            'strand': '+' if feature.location.strand == 1 else '-'
        }
        
        # Include protein_id if available
        if 'protein_id' in feature.qualifiers:
            gene_annotations[gene_id]['protein_id'] = feature.qualifiers['protein_id'][0]
        
        print(f"  {i}. {gene_id}: {start}-{end} ({end-start+1} bp) - {product}")
    
    return {
        'accession': record.id,
        'organism': record.annotations.get('organism', 'Unknown'),
        'length': len(record.seq),
        'genes': gene_annotations
    }

def create_gene_id(product_name):
    """Create a clean gene ID from product description"""
    
    # Common patterns to extract gene names
    product_lower = product_name.lower()
    
    # Handle structural proteins
    if 'capsid' in product_lower or 'core protein c' in product_lower:
        return 'C'
    elif 'envelope' in product_lower:
        return 'E'
    elif 'membrane glycoprotein m' in product_lower:
        return 'M'
    elif 'premembrane' in product_lower or 'prem protein' in product_lower:
        return 'prM'
    
    # Handle non-structural proteins
    elif 'ns1' in product_lower or 'nonstructural protein 1' in product_lower:
        return 'NS1'
    elif 'ns2a' in product_lower:
        return 'NS2A'
    elif 'ns2b' in product_lower:
        return 'NS2B'
    elif 'ns3' in product_lower:
        return 'NS3'
    elif 'ns4a' in product_lower:
        return 'NS4A'
    elif 'ns4b' in product_lower:
        return 'NS4B'
    elif 'ns5' in product_lower:
        return 'NS5'
    
    # Handle special cases
    elif 'protein pr' in product_lower:
        return 'pr'
    elif 'protein 2k' in product_lower:
        return '2K'
    elif 'anchored capsid' in product_lower:
        return 'ancC'
    
    # If no pattern matches, try to extract from the product name
    parts = product_name.split()
    for part in reversed(parts):  # Check from end first
        if part.upper() in ['C', 'E', 'M', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', 'NS4B', 'NS5']:
            return part.upper()
    
    # Default: use full product name
    return product_name.replace(' ', '_')

def classify_gene_type(product_name):
    """Classify gene as structural or nonstructural"""
    product_lower = product_name.lower()
    
    structural_keywords = ['capsid', 'envelope', 'membrane', 'core', 'anchored']
    nonstructural_keywords = ['nonstructural', 'ns1', 'ns2', 'ns3', 'ns4', 'ns5']
    
    for keyword in structural_keywords:
        if keyword in product_lower:
            return 'structural'
    
    for keyword in nonstructural_keywords:
        if keyword in product_lower:
            return 'nonstructural'
    
    return 'unknown'

def create_snpeff_files(annotation_data, output_dir):
    """Create files needed for SnpEff database update"""
    
    accession = annotation_data['accession']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Create genes.gff file for SnpEff
    gff_file = os.path.join(output_dir, f"{accession}_genes.gff")
    with open(gff_file, 'w') as f:
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {accession} 1 {annotation_data['length']}\n")
        
        # Write gene and CDS features
        for gene_id, gene_info in annotation_data['genes'].items():
            # Gene feature
            f.write(f"{accession}\tGenBank\tgene\t{gene_info['start']}\t{gene_info['end']}\t.\t{gene_info['strand']}\t.\tID=gene_{gene_id};Name={gene_id}\n")
            
            # CDS feature
            f.write(f"{accession}\tGenBank\tCDS\t{gene_info['start']}\t{gene_info['end']}\t.\t{gene_info['strand']}\t0\tID=cds_{gene_id};Parent=gene_{gene_id};Name={gene_info['product']}\n")
    
    print(f"\nCreated GFF file: {gff_file}")
    
    # 2. Create JSON summary for pipeline
    json_file = os.path.join(output_dir, f"{accession}_annotation_summary.json")
    with open(json_file, 'w') as f:
        json.dump(annotation_data, f, indent=2)
    
    print(f"Created JSON summary: {json_file}")
    
    # 3. Create visualization-ready gene coordinates
    vis_file = os.path.join(output_dir, f"{accession}_gene_coords.json")
    gene_coords = {gene_id: [info['start'], info['end']] 
                   for gene_id, info in annotation_data['genes'].items()}
    
    vis_data = {
        'accession': accession,
        'organism': annotation_data['organism'],
        'gene_coords': gene_coords,
        'structural_genes': [g for g, info in annotation_data['genes'].items() 
                            if info['type'] == 'structural'],
        'nonstructural_genes': [g for g, info in annotation_data['genes'].items() 
                               if info['type'] == 'nonstructural']
    }
    
    with open(vis_file, 'w') as f:
        json.dump(vis_data, f, indent=2)
    
    print(f"Created visualization file: {vis_file}")
    
    return gff_file

def main():
    parser = argparse.ArgumentParser(
        description='Parse mat_peptide annotations from GenBank file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
  python3 parse_genbank_annotations.py NC_001477.1.gb --output-dir NC_001477.1_annotations

This will create:
  - NC_001477.1_genes.gff: GFF file for SnpEff database update
  - NC_001477.1_annotation_summary.json: Complete annotation data
  - NC_001477.1_gene_coords.json: Gene coordinates for visualization
        '''
    )
    
    parser.add_argument('genbank_file', help='GenBank file to parse')
    parser.add_argument('--output-dir', default=None, 
                       help='Output directory (default: ACCESSION_annotations)')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.genbank_file):
        print(f"Error: GenBank file not found: {args.genbank_file}")
        sys.exit(1)
    
    # Parse GenBank file
    print(f"Parsing GenBank file: {args.genbank_file}")
    annotation_data = parse_mat_peptides(args.genbank_file)
    
    if not annotation_data:
        print("Error: No mat_peptide annotations found")
        sys.exit(1)
    
    # Set output directory
    if args.output_dir is None:
        args.output_dir = f"{annotation_data['accession']}_annotations"
    
    # Create output files
    gff_file = create_snpeff_files(annotation_data, args.output_dir)
    
    print(f"\nSuccess! Gene annotations extracted for {annotation_data['accession']}")
    print(f"Total genes annotated: {len(annotation_data['genes'])}")
    print(f"\nNext step: Update SnpEff database with:")
    print(f"  java -jar snpEff.jar build -gff3 -v {annotation_data['accession']}")

if __name__ == '__main__':
    main()