# Quick targeted fix - just comment out the POWV hardcoded section
with open('viral_pipeline/visualization/visualize_mutations_python_outdir.py', 'r') as f:
    content = f.read()

# Replace the POWV hardcoded section with a simple pass
old_section = '''        if accession == 'HM440560.1' or accession == 'HM440560':
            # Hardcoded Powassan virus gene coordinates
            gene_coords = {
                'C': [110, 357],
                'prM': [358, 858],
                'Env': [859, 2361],
                'NS1': [2362, 3417],
                'NS2a': [3418, 4110],
                'NS2b': [4111, 4503],
                'NS3': [4504, 6360],
                'NS4a': [6361, 6807],
                'NS4b': [6808, 7563],
                'NS5': [7564, 10287]
            }
            # Use standard flavivirus colors
            gene_colors = {
                'C': '#4575b4',
                'prM': '#74add1',
                'Env': '#abd9e9',
                'NS1': '#fee090',
                'NS2a': '#fdae61',
                'NS2b': '#f46d43',
                'NS3': '#d73027',
                'NS4a': '#a50026',
                'NS4b': '#762a83',
                'NS5': '#5aae61'
            }
            structural_genes = ['C', 'prM', 'Env']
            nonstructural_genes = ['NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5']
            print('Using hardcoded Powassan virus gene coordinates')'''

new_section = '''        # POWV now uses database configuration (removed hardcoded coordinates)
        pass'''

content = content.replace(old_section, new_section)

with open('viral_pipeline/visualization/visualize_mutations_python_outdir.py', 'w') as f:
    f.write(content)

print(Fixed POWV - now uses database)
