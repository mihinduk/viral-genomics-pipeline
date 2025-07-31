import sys
sys.path.append('viral_pipeline/visualization')
from visualize_depth import get_gene_info

# Check what colors are assigned
gene_info = get_gene_info('NC_001477.1')
gene_coords, gene_colors, structural_genes, nonstructural_genes, display_names = gene_info

print("Gene coordinates and colors:")
key_genes = ['ancC', 'C', 'prM', 'pr', 'M']
for gene in key_genes:
    if gene in gene_coords:
        coords = gene_coords[gene]
        color = gene_colors.get(gene, 'NOT_FOUND')
        display = display_names.get(gene, gene)
        print(f"  {gene}: coords={coords}, color={color}, display='{display}'")
