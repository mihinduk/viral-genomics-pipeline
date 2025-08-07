#!/usr/bin/env python3
"""
Viral Genomics Pipeline
End-to-end viral genomics: from raw reads to publication-ready visualizations

The only integrated analysis-to-visualization toolkit for defensive viral surveillance.
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
# Configuration - can be overridden with environment variables
SHOTGUN_PIPELINE_PATH = os.environ.get("SHOTGUN_PIPELINE_PATH", "{SHOTGUN_PIPELINE_PATH}")
MAMBA_PATH = os.environ.get("MAMBA_PATH", "{MAMBA_PATH}")


def run_command(cmd, description=""):
    """Run shell command with error handling"""
    print(f"\n🔄 {description}")
    print(f"   Command: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ Error: {result.stderr}")
        sys.exit(1)
    print(f"✅ {description} completed")
    return result.stdout

def main():
    parser = argparse.ArgumentParser(
        description="End-to-end viral genomics pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete pipeline for DENV1
  python run_viral_pipeline.py \\
    --r1 sample_R1.fastq.gz \\
    --r2 sample_R2.fastq.gz \\
    --accession NC_001477.1 \\
    --threads 4 \\
    --quality-cutoff 1000 \\
    --freq-cutoff 0.01 \\
    --output sample_results

  # Visualization only (if you have existing data)
  python run_viral_pipeline.py \\
    --visualization-only \\
    --vcf sample.vcf \\
    --bam sample.bam \\
    --accession NC_001477.1 \\
    --quality-cutoff 1000 \\
    --freq-cutoff 0.01 \\
    --output sample_results
        """
    )
    
    # Input files
    parser.add_argument('--r1', help='Forward reads FASTQ file')
    parser.add_argument('--r2', help='Reverse reads FASTQ file')
    parser.add_argument('--vcf', help='Input VCF file (for visualization-only mode)')
    parser.add_argument('--bam', help='Input BAM file (for depth visualization)')
    
    # Parameters
    parser.add_argument('--accession', required=True, help='GenBank accession (e.g., NC_001477.1)')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    parser.add_argument('--quality-cutoff', type=int, default=1000, help='Quality score cutoff (default: 1000)')
    parser.add_argument('--freq-cutoff', type=float, default=0.01, help='Allele frequency cutoff for visualization (default: 0.01)')
    parser.add_argument('--output', required=True, help='Output directory prefix')
    
    # Modes
    parser.add_argument('--visualization-only', action='store_true', help='Skip analysis, only create visualizations')
    parser.add_argument('--analysis-only', action='store_true', help='Skip visualization, only run analysis')
    
    args = parser.parse_args()
    
    # Get pipeline directory
    pipeline_dir = Path(__file__).parent
    
    print("🦠 Viral Genomics Pipeline")
    print("=" * 50)
    print(f"📁 Accession: {args.accession}")
    print(f"📊 Quality cutoff: {args.quality_cutoff}")
    print(f"📈 Frequency cutoff: {args.freq_cutoff}")
    print(f"📂 Output: {args.output}")
    
    # Create output directory
    output_dir = f"{args.output}_results"
    os.makedirs(output_dir, exist_ok=True)
    
    if not args.visualization_only:
        # Run analysis pipeline
        if not args.r1 or not args.r2:
            print("❌ Error: --r1 and --r2 required for analysis")
            sys.exit(1)
            
        analysis_cmd = f"""{SHOTGUN_PIPELINE_PATH} \\
            "{args.r1}" \\
            "{args.r2}" \\
            {args.accession} \\
            {args.threads}"""
        
        run_command(analysis_cmd, "Running viral genomics analysis")
        
        # Set paths for visualization
        vcf_file = f"cleaned_seqs/variants/*{args.accession.replace('.', '_')}*.snpEFF.ann.vcf"
        bam_file = f"cleaned_seqs/mapping/*.lofreq.final.bam"
    else:
        if not args.vcf:
            print("❌ Error: --vcf required for visualization-only mode")
            sys.exit(1)
        vcf_file = args.vcf
        bam_file = args.bam
    
    if not args.analysis_only:
        # Generate depth file if BAM provided
        if bam_file:
            depth_cmd = f"""{MAMBA_PATH} run -n viral_genomics \\
                echo -e "chrom\tposition\tdepth" > {output_dir}/{args.output}_depth.txt && samtools depth {bam_file} >> {output_dir}/{args.output}_depth.txt"""
            run_command(depth_cmd, "Generating depth file")
        
        # Parse VCF with quality filtering
        parse_cmd = f"""{MAMBA_PATH} run -n viral_genomics \\
            python3 {pipeline_dir}/viral_pipeline/visualization/parse_snpeff_vcf.py \\
            -i {vcf_file} \\
            -d 200 \\
            -q {args.quality_cutoff} \\
            -o {args.output}_filtered_mutations.tsv \\
            -O {output_dir}"""
        run_command(parse_cmd, "Parsing and filtering VCF")
        
        # Create mutation visualization
        mutation_cmd = f"""{MAMBA_PATH} run -n viral_genomics \\
            python3 {pipeline_dir}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \\
            --input {output_dir}/{args.output}_filtered_mutations.tsv \\
            --output {args.output}_mutations.png \\
            --accession {args.accession} \\
            --cutoff {args.freq_cutoff} \\
            --outdir {output_dir}"""
        run_command(mutation_cmd, "Creating mutation visualization")
        
        # Create depth visualization if depth file exists
        depth_file = f"{output_dir}/{args.output}_depth.txt"
        if os.path.exists(depth_file):
            depth_cmd = f"""{MAMBA_PATH} run -n viral_genomics \\
                python3 {pipeline_dir}/viral_pipeline/visualization/visualize_depth.py \\
                --depth {depth_file} \\
                --output {output_dir}/{args.output}_depth.png \\
                --output-html {output_dir}/{args.output}_depth.html \\
                --accession {args.accession}"""
            run_command(depth_cmd, "Creating depth visualization")
    
    print("\n🎉 Pipeline completed successfully!")
    print(f"📁 Results in: {output_dir}")
    print("   - Mutation visualization with perfect gene layout")
    print("   - Depth coverage visualization")
    print("   - Interactive HTML reports")
    print("   - Complete mutation tables")

if __name__ == "__main__":
    main()