import sys
sys.path.append('viral_pipeline/visualization')
from visualize_depth import get_gene_info

# Check DENV1 color assignments
gene_info = get_gene_info('NC_001477.1')
gene_coords, gene_colors, structural_genes, nonstructural_genes, display_names = gene_info

print("DENV1 Gene Colors:")
for gene in ['pr', 'M', 'prM']:
    if gene in gene_colors:
        print(f"  {gene}: {gene_colors[gene]}")
    else:
        print(f"  {gene}: NOT FOUND")
        
print("\nAll genes in gene_colors:")
for gene, color in gene_colors.items():
    print(f"  {gene}: {color}")
