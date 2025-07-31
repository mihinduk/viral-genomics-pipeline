import re

# Read the file
with open('viral_pipeline/visualization/visualize_depth.py', 'r') as f:
    content = f.read()

# Find and replace the sorting logic section
old_section = r'''    # Create a mapping of labels to their genomic start positions
    label_positions = {}
    for gene, \(start, end\) in gene_coords\.items\(\):
        display_name = display_names\.get\(gene, gene\)
        label_positions\[display_name\] = start
    
    # Sort handles and labels by genomic position
    sorted_pairs = sorted\(zip\(handles, labels\), key=lambda x: label_positions\.get\(x\[1\], float\('inf'\)\)\)'''

new_section = '''    # Create a mapping of labels to their genomic positions (start, end)
    label_positions = {}
    for gene, (start, end) in gene_coords.items():
        display_name = display_names.get(gene, gene)
        label_positions[display_name] = (start, end)
    
    # Sort handles and labels by genomic position (start, then end)
    sorted_pairs = sorted(zip(handles, labels), key=lambda x: label_positions.get(x[1], (float("inf"), float("inf"))))'''

# Replace the section
content = re.sub(old_section, new_section, content, flags=re.MULTILINE)

# Write back
with open('viral_pipeline/visualization/visualize_depth.py', 'w') as f:
    f.write(content)

print("Legend sorting logic updated\!")
