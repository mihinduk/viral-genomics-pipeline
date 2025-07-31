import re

# Read the file
with open('viral_pipeline/utils/generate_filtered_consensus.py', 'r') as f:
    content = f.read()

# Replace the problematic line with logic that handles both 2 and 3 element keys
old_line = r'    for \(gene_info, protein_seq, mutations\), protein_data in unique_proteins\.items\(\):'

new_section = '''    for key, protein_data in unique_proteins.items():
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
            continue'''

# Replace the line
content = re.sub(old_line, new_section, content, flags=re.MULTILINE)

# Write back
with open('viral_pipeline/utils/generate_filtered_consensus.py', 'w') as f:
    f.write(content)

print("Fixed ZIKV consensus generation unpacking\!")
