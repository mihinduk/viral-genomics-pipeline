import re

# Read the original file
with open('viral_pipeline/visualization/visualize_depth.py', 'r') as f:
    content = f.read()

# Define the new function
new_function = '''def read_depth_file(depth_file):
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
    return df'''

# Replace the old function with the new one
pattern = r'def read_depth_file\(depth_file\):.*?return df'
content = re.sub(pattern, new_function, content, flags=re.DOTALL)

# Write back to file
with open('viral_pipeline/visualization/visualize_depth.py', 'w') as f:
    f.write(content)

print("Function replaced successfully\!")
