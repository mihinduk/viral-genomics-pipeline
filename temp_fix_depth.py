import pandas as pd

def read_depth_file(depth_file):
    """Read depth file created by samtools depth"""
    print(f"📖 Reading depth file: {depth_file}")
    
    # Read first line to check if it's a header
    with open(depth_file, 'r') as f:
        first_line = f.readline().strip()
    
    # Check if first line starts with 'chrom' (header)
    if first_line.startswith('chrom'):
        # Skip header line
        df = pd.read_csv(depth_file, sep="\t", header=0, 
                         names=["reference", "position", "depth"],
                         dtype={"position": int, "depth": int})
    else:
        # No header, read directly
        df = pd.read_csv(depth_file, sep="\t", header=None, 
                         names=["reference", "position", "depth"],
                         dtype={"position": int, "depth": int})
    
    print(f"   Found {len(df)} positions")
    return df

# Test with both files
print("Testing SMS_5 (headerless):")
df5 = read_depth_file("/scratch/sahlab/kathie/NovaSeq_N1027_DENV1/SMS_5_DENV1_results/SMS_5_DENV1_depth.txt")
print(f"SMS_5 columns: {list(df5.columns)}")
print(f"SMS_5 first 3 rows:\n{df5.head(3)}")

print("\nTesting SMS_6 (with header):")
df6 = read_depth_file("/scratch/sahlab/kathie/NovaSeq_N1027_DENV1/SMS_6_DENV1_results/SMS_6_DENV1_depth.txt")
print(f"SMS_6 columns: {list(df6.columns)}")
print(f"SMS_6 first 3 rows:\n{df6.head(3)}")
