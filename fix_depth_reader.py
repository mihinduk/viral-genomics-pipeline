# Fix the depth file reader to handle headers properly

with open("viral_pipeline/visualization/visualize_depth.py", "r") as f:
    lines = f.readlines()

# Find and replace the line that reads the depth file
for i, line in enumerate(lines):
    if "pd.read_csv(depth_file, sep=" in line and "header=None" in line:
        # Replace with header-aware version
        lines[i] = """    # Read depth file - handle both with and without headers
    first_line = open(depth_file).readline().strip()
    if first_line.startswith(#):
        df = pd.read_csv(depth_file, sep=\t, header=0, comment=#)
        df.columns = [chrom, position, depth]  # Standardize column names
    else:
        df = pd.read_csv(depth_file, sep=\t, header=None, names=[chrom, position, depth])
"""
        break

# Write the fixed version
with open("viral_pipeline/visualization/visualize_depth.py", "w") as f:
    f.writelines(lines)

print("✅ Fixed depth file reader to handle headers properly")
