#!/usr/bin/env python3
"""
Update SnpEff database with BLAST-derived gene annotations
Replaces generic polyprotein annotation with individual protein annotations
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path
import tempfile

def find_snpeff_directory():
    """Find SnpEff installation directory"""
    # Common SnpEff locations
    possible_paths = [
        "/home/mihindu/software/snpEff",
        "/scratch/sahlab/kathie/software/snpEff",
        "/opt/snpEff",
        "/usr/local/snpEff"
    ]
    
    for path in possible_paths:
        snpeff_path = Path(path)
        if snpeff_path.exists() and (snpeff_path / "snpEff.jar").exists():
            return snpeff_path
    
    return None

def backup_existing_database(snpeff_dir, accession):
    """Backup existing database if it exists"""
    db_path = snpeff_dir / "data" / accession
    if db_path.exists():
        backup_path = db_path.with_suffix('.backup')
        if backup_path.exists():
            shutil.rmtree(backup_path)
        shutil.copytree(db_path, backup_path)
        print(f"  Backed up existing database to: {backup_path}")
        return backup_path
    return None

def update_snpeff_config(snpeff_dir, accession, genome_file):
    """Add or update entry in snpEff.config"""
    config_file = snpeff_dir / "snpEff.config"
    
    # Read existing config
    if config_file.exists():
        with open(config_file, 'r') as f:
            lines = f.readlines()
    else:
        lines = []
    
    # Check if entry already exists
    entry_exists = False
    for i, line in enumerate(lines):
        if line.strip().startswith(f"{accession}.genome"):
            lines[i] = f"{accession}.genome : {accession}\n"
            entry_exists = True
            break
    
    # Add entry if it doesn't exist
    if not entry_exists:
        lines.append(f"\n# Custom virus database\n")
        lines.append(f"{accession}.genome : {accession}\n")
    
    # Write updated config
    with open(config_file, 'w') as f:
        f.writelines(lines)
    
    print(f"  Updated SnpEff config: {config_file}")

def create_snpeff_database(snpeff_dir, accession, genome_file, gff_file):
    """Create SnpEff database with new annotation"""
    
    # Create database directory
    db_dir = snpeff_dir / "data" / accession
    db_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy files to database directory
    sequences_file = db_dir / "sequences.fa"
    genes_file = db_dir / "genes.gff"
    
    shutil.copy2(genome_file, sequences_file)
    shutil.copy2(gff_file, genes_file)
    
    print(f"  Copied genome to: {sequences_file}")
    print(f"  Copied GFF to: {genes_file}")
    
    # Build database using mamba
    build_cmd = f"""
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
    java -jar {snpeff_dir}/snpEff.jar build -gff3 -v {accession}
"""
    
    print(f"  Building SnpEff database...")
    try:
        result = subprocess.run(build_cmd, shell=True, capture_output=True, text=True, cwd=str(snpeff_dir))
        if result.returncode != 0:
            print(f"Error building database: {result.stderr}")
            return False
        print("  Database built successfully!")
        return True
    except Exception as e:
        print(f"Error running SnpEff build: {e}")
        return False

def test_database(snpeff_dir, accession):
    """Test the database by running SnpEff info"""
    test_cmd = f"""
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
    java -jar {snpeff_dir}/snpEff.jar databases | grep {accession}
"""
    
    try:
        result = subprocess.run(test_cmd, shell=True, capture_output=True, text=True, cwd=str(snpeff_dir))
        if accession in result.stdout:
            print(f"  ✅ Database {accession} is available")
            return True
        else:
            print(f"  ⚠️  Database {accession} may not be properly installed")
            return False
    except Exception as e:
        print(f"  Error testing database: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Update SnpEff database with BLAST-derived annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Update ZIKV database with new annotations
  python3 update_snpeff_database.py \\
    --accession KU955591.1 \\
    --genome cleaned_seqs/KU955591.1.fasta \\
    --gff blast_annotation_KU955591.1/KU955591.1.gff3
  
  # Specify custom SnpEff directory
  python3 update_snpeff_database.py \\
    --accession KU955591.1 \\
    --genome cleaned_seqs/KU955591.1.fasta \\
    --gff blast_annotation_KU955591.1/KU955591.1.gff3 \\
    --snpeff-dir /path/to/snpEff
        """
    )
    
    parser.add_argument('--accession', required=True, help='Virus accession (e.g., KU955591.1)')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--gff', required=True, help='Path to GFF3 annotation file')
    parser.add_argument('--snpeff-dir', help='Path to SnpEff directory (auto-detected if not provided)')
    parser.add_argument('--backup', action='store_true', help='Backup existing database before updating')
    parser.add_argument('--test', action='store_true', help='Test database after creation')
    
    args = parser.parse_args()
    
    print("SnpEff Database Updater")
    print("=" * 50)
    
    # Validate input files
    genome_path = Path(args.genome)
    gff_path = Path(args.gff)
    
    if not genome_path.exists():
        print(f"❌ Error: Genome file not found: {args.genome}")
        return 1
    
    if not gff_path.exists():
        print(f"❌ Error: GFF file not found: {args.gff}")
        return 1
    
    # Find SnpEff directory
    if args.snpeff_dir:
        snpeff_dir = Path(args.snpeff_dir)
    else:
        snpeff_dir = find_snpeff_directory()
    
    if not snpeff_dir or not snpeff_dir.exists():
        print(f"❌ Error: SnpEff directory not found")
        print("Available options:")
        print("  1. Specify path with --snpeff-dir")
        print("  2. Install SnpEff in a standard location")
        return 1
    
    snpeff_jar = snpeff_dir / "snpEff.jar"
    if not snpeff_jar.exists():
        print(f"❌ Error: snpEff.jar not found in {snpeff_dir}")
        return 1
    
    print(f"Using SnpEff directory: {snpeff_dir}")
    print(f"Accession: {args.accession}")
    print(f"Genome: {genome_path}")
    print(f"GFF: {gff_path}")
    
    # Backup existing database
    if args.backup:
        print(f"\nStep 1: Backing up existing database...")
        backup_path = backup_existing_database(snpeff_dir, args.accession)
    
    # Update SnpEff config
    print(f"\nStep 2: Updating SnpEff configuration...")
    update_snpeff_config(snpeff_dir, args.accession, genome_path)
    
    # Create database
    print(f"\nStep 3: Creating SnpEff database...")
    if not create_snpeff_database(snpeff_dir, args.accession, genome_path, gff_path):
        print("❌ Failed to create database")
        return 1
    
    # Test database
    if args.test:
        print(f"\nStep 4: Testing database...")
        test_database(snpeff_dir, args.accession)
    
    print(f"\n✅ SnpEff database update complete!")
    print(f"\nNext steps:")
    print(f"  1. Re-run variant calling on your samples")
    print(f"  2. Check that variants are now annotated with individual protein names")
    print(f"  3. Run Module 8 for proper consensus generation")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())