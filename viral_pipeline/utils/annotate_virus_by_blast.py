#!/usr/bin/env python3
"""
Module 10: Annotate Virus by BLAST
Uses BLASTX against RefSeq proteins to identify gene coordinates for viruses
with only polyprotein annotations in GenBank
"""

import argparse
import json
import sys
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO, Entrez
# Removed NCBIXMLT import - using tabular output instead
from collections import defaultdict
import re

# Set email for NCBI
Entrez.email = "viral_pipeline@example.com"

def get_refseq_proteins(refseq_accession, output_dir):
    """Download RefSeq proteins for a virus"""
    print(f"Downloading RefSeq proteins for {refseq_accession}...")
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download GenBank record
    handle = Entrez.efetch(db="nucleotide", id=refseq_accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    # Extract protein sequences from both CDS and mat_peptide features
    proteins = {}
    protein_coords = {}
    
    # Debug: show all feature types
    feature_types = set(f.type for f in record.features)
    print(f"  Feature types in GenBank: {', '.join(sorted(feature_types))}")
    
    # Count features
    cds_count = sum(1 for f in record.features if f.type == "CDS")
    mat_peptide_count = sum(1 for f in record.features if f.type == "mat_peptide")
    print(f"  CDS features: {cds_count}, mat_peptide features: {mat_peptide_count}")
    
    for feature in record.features:
        if feature.type in ["CDS", "mat_peptide"]:
            # Get protein info
            gene_name = None
            protein_seq = None
            
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
                if feature.type == "mat_peptide":
                    # For mat_peptide, always show what we're processing
                    print(f"  Processing mat_peptide: {product}")
                
                # Skip polyprotein - it's not informative
                if 'polyprotein' in product.lower():
                    print(f"    Skipping polyprotein: {product}")
                    continue
                    
                # Extract gene name from product (prioritize this over /gene for mat_peptide)
                gene_name = extract_gene_name(product)
            elif 'gene' in feature.qualifiers and feature.type != "mat_peptide":
                # Only use /gene for non-mat_peptide features
                gene_name = feature.qualifiers['gene'][0]
            else:
                # Debug missing qualifiers
                print(f"  Feature {feature.type} has no gene/product qualifiers")
                # Check all qualifiers for mat_peptide
                if feature.type == "mat_peptide":
                    print(f"    Available qualifiers: {list(feature.qualifiers.keys())}")
                    for qual_key, qual_values in feature.qualifiers.items():
                        print(f"      {qual_key}: {qual_values}")
            
            # Get protein sequence
            if 'translation' in feature.qualifiers:
                protein_seq = feature.qualifiers['translation'][0]
            elif 'peptide' in feature.qualifiers:
                protein_seq = feature.qualifiers['peptide'][0]
            elif feature.type == "mat_peptide":
                # For mat_peptide, extract the nucleotide sequence and translate
                try:
                    # Find the parent CDS feature to get the protein sequence
                    parent_location = None
                    for f in record.features:
                        if f.type == "CDS" and 'translation' in f.qualifiers:
                            parent_location = f.location
                            parent_translation = f.qualifiers['translation'][0]
                            break
                    
                    if parent_location and parent_translation:
                        # Calculate the position within the CDS
                        mat_start = int(feature.location.start)
                        mat_end = int(feature.location.end)
                        cds_start = int(parent_location.start)
                        
                        # Calculate protein positions (3 nucleotides = 1 amino acid)
                        prot_start = (mat_start - cds_start) // 3
                        prot_end = (mat_end - cds_start) // 3
                        
                        # Extract the protein subsequence
                        protein_seq = parent_translation[prot_start:prot_end]
                        # gene_name should already be set from product, but let's make sure
                        if not gene_name:
                            gene_name = f"mat_peptide_{prot_start}_{prot_end}"
                        print(f"    Extracted mat_peptide {gene_name}: {prot_start}-{prot_end} ({len(protein_seq)} aa)")
                    else:
                        # Fallback: translate the nucleotide sequence
                        seq = feature.extract(record.seq)
                        protein_seq = str(seq.translate())
                except Exception as e:
                    print(f"  Error extracting mat_peptide: {e}")
                    continue
            
            if gene_name and protein_seq:
                # Normalize gene names
                gene_name = normalize_gene_name(gene_name)
                # For mat_peptides from polyprotein, make sure each has unique name
                if gene_name in proteins and feature.type == "mat_peptide":
                    print(f"    Warning: Duplicate gene name {gene_name}, keeping first occurrence")
                else:
                    proteins[gene_name] = protein_seq
                    protein_coords[gene_name] = (int(feature.location.start) + 1, int(feature.location.end))
    
    # Write protein sequences
    protein_file = output_dir / f"{refseq_accession}_proteins.fasta"
    with open(protein_file, 'w') as f:
        for gene, seq in proteins.items():
            f.write(f">{gene}\n{seq}\n")
    
    print(f"  Downloaded {len(proteins)} proteins: {', '.join(proteins.keys())}")
    return protein_file, proteins, protein_coords, record.description

def extract_gene_name(product_desc):
    """Extract gene name from product description - use full name with underscores"""
    # Clean up the product description
    cleaned = product_desc.strip()
    
    # Remove common prefixes/suffixes that aren't part of the actual name
    cleaned = re.sub(r'^(putative|hypothetical)\s+', '', cleaned, flags=re.IGNORECASE)
    
    # Replace spaces with underscores for valid identifiers
    cleaned = cleaned.replace(' ', '_')
    
    # Remove special characters that might cause issues
    cleaned = re.sub(r'[^\w_-]', '', cleaned)
    
    return cleaned

def normalize_gene_name(gene_name):
    """Normalize gene names to standard format"""
    gene_name = gene_name.strip()
    
    # Envelope
    if gene_name.upper() == 'E':
        return 'Env'
    
    # NS genes
    if gene_name.upper().startswith('NS') and len(gene_name) > 2:
        # Normalize NS2a, NS2A -> NS2a
        base = gene_name[:3]
        if len(gene_name) == 4 and gene_name[3].isalpha():
            return base + gene_name[3].lower()
        return gene_name
    
    return gene_name

def run_tblastn(protein_file, genome_file, output_file):
    """Run tBLASTn to find protein locations in genome"""
    print(f"\nRunning tBLASTn...")
    
    # Create BLAST database
    db_name = genome_file + ".blastdb"
    
    # Setup commands to run with mamba
    # First check if BLAST is available, if not, install it
    makeblastdb_cmd = f"""
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics which makeblastdb > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Installing BLAST+ in viral_genomics environment..."
    /home/mihindu/miniforge3/bin/mamba install -n viral_genomics -y -c bioconda blast
fi
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics makeblastdb -in {genome_file} -dbtype nucl -out {db_name}
"""
    
    # Create database
    try:
        result = subprocess.run(makeblastdb_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error creating BLAST database: {result.stderr}")
            return False
    except Exception as e:
        print(f"Error creating BLAST database: {e}")
        return False
    
    # Run tBLASTn
    tblastn_cmd = f"""
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics tblastn -query {protein_file} -db {db_name} -out {output_file} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -max_target_seqs 1 -evalue 1e-5
"""
    
    try:
        result = subprocess.run(tblastn_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error running tBLASTn: {result.stderr}")
            return False
        print("  tBLASTn completed successfully")
        return True
    except Exception as e:
        print(f"Error running tBLASTn: {e}")
        return False

def parse_blast_results(blast_file, min_coverage=0.8):
    """Parse BLAST results to get gene coordinates"""
    gene_coords = {}
    
    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 13:
                gene = parts[0]
                pident = float(parts[2])
                length = int(parts[3])
                qlen = int(parts[12])
                sstart = int(parts[8])
                send = int(parts[9])
                
                # Calculate coverage
                coverage = length / qlen
                
                # Only accept high-quality matches
                if pident >= 80 and coverage >= min_coverage:
                    # Handle reverse strand
                    if sstart > send:
                        start, end = send, sstart
                    else:
                        start, end = sstart, send
                    
                    # Store best match for each gene
                    if gene not in gene_coords or pident > gene_coords[gene][2]:
                        gene_coords[gene] = (start, end, pident, coverage)
    
    # Convert to simple coordinate dict
    final_coords = {}
    for gene, (start, end, pident, coverage) in gene_coords.items():
        final_coords[gene] = [start, end]
        print(f"  {gene}: {start}-{end} (identity: {pident:.1f}%, coverage: {coverage*100:.1f}%)")
    
    return final_coords

def determine_virus_family(gene_names):
    """Determine virus family based on gene names"""
    gene_str = ' '.join(gene_names).lower()
    
    if any(g in gene_str for g in ['ns1', 'ns2a', 'ns2b', 'ns3', 'ns4a', 'ns4b', 'ns5']):
        return 'flavivirus'
    elif any(g in gene_str for g in ['nsp1', 'nsp2', 'nsp3', 'nsp4']):
        return 'alphavirus'
    else:
        return 'flavivirus'  # default

def get_gene_colors(family, genes):
    """Get standard colors for genes"""
    if family == 'flavivirus':
        color_map = {
            'C': '#4575b4',
            'prM': '#74add1',
            'M': '#74add1',
            'Env': '#abd9e9',
            'NS1': '#fdae61',
            'NS2a': '#f46d43',
            'NS2b': '#d73027',
            'NS3': '#a50026',
            'NS4a': '#762a83',
            'NS4b': '#9970ab',
            'NS5': '#c2a5cf'
        }
    else:  # alphavirus
        color_map = {
            'C': '#4575b4',
            'E3': '#5ba5c9',
            'E2': '#7fb3d3',
            'E1': '#abd9e9',
            'nsP1': '#fdae61',
            'nsP2': '#f46d43',
            'nsP3': '#d73027',
            'nsP4': '#a50026'
        }
    
    colors = {}
    structural_genes = []
    nonstructural_genes = []
    
    for gene in genes:
        colors[gene] = color_map.get(gene, '#808080')
        
        if gene in ['C', 'prM', 'M', 'Env', 'E1', 'E2', 'E3']:
            structural_genes.append(gene)
        else:
            nonstructural_genes.append(gene)
    
    return colors, structural_genes, nonstructural_genes

def create_snpeff_gff(accession, gene_coords, output_file):
    """Create GFF3 file for SnpEff with proper gene names"""
    with open(output_file, 'w') as f:
        # Write header
        f.write("##gff-version 3\n")
        
        # Sort genes by start position
        sorted_genes = sorted(gene_coords.items(), key=lambda x: x[1][0])
        
        # Write gene and CDS features
        for i, (gene_name, (start, end)) in enumerate(sorted_genes, 1):
            gene_id = f"{gene_name}_gene"
            transcript_id = f"{gene_name}_mRNA"
            
            # Gene feature
            f.write(f"{accession}\tBLAST\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id};Name={gene_name}\n")
            
            # mRNA feature
            f.write(f"{accession}\tBLAST\tmRNA\t{start}\t{end}\t.\t+\t.\tID={transcript_id};Parent={gene_id};Name={gene_name}\n")
            
            # CDS feature
            f.write(f"{accession}\tBLAST\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene_name}_cds;Parent={transcript_id};Name={gene_name}\n")
    
    print(f"\nCreated GFF3 file for SnpEff: {output_file}")
    return output_file

def update_virus_config(accession, virus_name, gene_coords, family, config_file):
    """Update known_viruses.json with BLAST-derived coordinates"""
    
    # Get colors and gene classifications
    colors, structural, nonstructural = get_gene_colors(family, gene_coords.keys())
    
    # Create virus entry
    virus_info = {
        "name": virus_name,
        "family": family,
        "genome_length": 0,  # Will be updated if we have the genome
        "gene_coords": gene_coords,
        "colors": colors,
        "structural_genes": structural,
        "nonstructural_genes": nonstructural
    }
    
    # Read existing config
    if config_file.exists():
        with open(config_file, 'r') as f:
            known_viruses = json.load(f)
    else:
        known_viruses = {}
    
    # Update entry
    known_viruses[accession] = virus_info
    
    # Write back
    with open(config_file, 'w') as f:
        json.dump(known_viruses, f, indent=2, sort_keys=True)
    
    print(f"\n✅ Updated {config_file} with BLAST-derived gene coordinates")

def main():
    parser = argparse.ArgumentParser(
        description='Annotate virus genes using BLAST against RefSeq proteins',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annotate ZIKV using RefSeq NC_012532.1 as reference
  python3 annotate_virus_by_blast.py KU955591.1.fasta --refseq NC_012532.1
  
  # Specify custom paths
  python3 annotate_virus_by_blast.py genome.fasta --refseq NC_012532.1 --name "My virus"
  
  # Update SnpEff database
  python3 annotate_virus_by_blast.py genome.fasta --refseq NC_012532.1 --update-snpeff
        """
    )
    
    parser.add_argument('genome', help='Path to genome FASTA file to annotate')
    parser.add_argument('--refseq', required=True, help='RefSeq accession with proper gene annotation')
    parser.add_argument('--name', help='Virus name (default: from RefSeq)')
    parser.add_argument('--config', help='Path to known_viruses.json to update')
    parser.add_argument('--output-dir', default='blast_annotation', help='Output directory')
    parser.add_argument('--update-snpeff', action='store_true', help='Create GFF3 for SnpEff update')
    parser.add_argument('--min-coverage', type=float, default=0.8, help='Minimum protein coverage (default: 0.8)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Virus Gene Annotation by BLAST")
    print("=" * 50)
    
    # Get genome accession from filename if possible
    genome_path = Path(args.genome)
    # Keep full accession including version (e.g., KU955591.1)
    accession = genome_path.stem
    
    # Step 1: Download RefSeq proteins
    try:
        protein_file, proteins, ref_coords, ref_name = get_refseq_proteins(
            args.refseq, output_dir
        )
    except Exception as e:
        print(f"❌ Error downloading RefSeq proteins: {e}")
        return 1
    
    # Step 2: Run tBLASTn
    blast_output = output_dir / f"{accession}_vs_{args.refseq}.blast"
    if not run_tblastn(str(protein_file), args.genome, str(blast_output)):
        return 1
    
    # Step 3: Parse BLAST results
    print(f"\nParsing BLAST results (minimum coverage: {args.min_coverage*100}%)...")
    gene_coords = parse_blast_results(blast_output, args.min_coverage)
    
    if not gene_coords:
        print("❌ No significant matches found!")
        return 1
    
    print(f"\n✅ Found coordinates for {len(gene_coords)} genes")
    
    # Step 4: Determine virus family
    family = determine_virus_family(gene_coords.keys())
    print(f"Virus family: {family}")
    
    # Step 5: Update configuration
    if args.config:
        config_file = Path(args.config)
    else:
        config_file = Path("/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/visualization/known_viruses.json")
    
    virus_name = args.name or f"{accession} (annotated by BLAST)"
    update_virus_config(accession, virus_name, gene_coords, family, config_file)
    
    # Step 6: Create GFF for SnpEff if requested
    if args.update_snpeff:
        gff_file = output_dir / f"{accession}.gff3"
        create_snpeff_gff(accession, gene_coords, gff_file)
        
        print("\nTo update SnpEff database:")
        print(f"  1. Copy {gff_file} to SnpEff data directory")
        print(f"  2. Rebuild the database with the new annotation")
    
    # Save results summary
    summary_file = output_dir / f"{accession}_annotation_summary.json"
    summary = {
        "accession": accession,
        "name": virus_name,
        "family": family,
        "refseq_used": args.refseq,
        "genes_found": len(gene_coords),
        "gene_coordinates": gene_coords
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n📊 Summary saved to: {summary_file}")
    print("\n🎉 Annotation complete!")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())