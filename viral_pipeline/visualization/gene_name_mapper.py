#!/usr/bin/env python3
"""
Gene Name Mapper for Viral Genomics Pipeline
Maps various gene/protein names from SnpEff output to standardized names for visualization
"""

import json
from pathlib import Path

class GeneNameMapper:
    def __init__(self, standards_file=None):
        """Initialize with virus gene standards"""
        if standards_file is None:
            # Try to find the standards file
            possible_paths = [
                Path(__file__).parent / "virus_gene_standards.json",
                Path(__file__).parent.parent / "visualization" / "virus_gene_standards.json",
                Path("/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline/viral_pipeline/visualization/virus_gene_standards.json"),
                Path("virus_gene_standards.json")
            ]
            
            for path in possible_paths:
                if path.exists():
                    standards_file = path
                    break
        
        if standards_file and Path(standards_file).exists():
            with open(standards_file, 'r') as f:
                self.standards = json.load(f)
        else:
            # Fallback to hardcoded minimal mapping
            self.standards = self._get_default_standards()
        
        # Build reverse mapping for quick lookup
        self._build_alias_mapping()
    
    def _build_alias_mapping(self):
        """Build reverse mapping from aliases to standard names"""
        self.alias_to_standard = {}
        
        for family, family_data in self.standards.items():
            self.alias_to_standard[family] = {}
            
            for gene_key, gene_data in family_data.get('standard_genes', {}).items():
                # Map the standard name to itself
                self.alias_to_standard[family][gene_key] = gene_key
                
                # Map all aliases to the standard name
                for alias in gene_data.get('aliases', []):
                    self.alias_to_standard[family][alias] = gene_key
    
    def get_standard_name(self, gene_name, virus_family='flavivirus'):
        """Get standardized gene name for visualization"""
        if virus_family not in self.alias_to_standard:
            return gene_name
        
        family_mapping = self.alias_to_standard[virus_family]
        return family_mapping.get(gene_name, gene_name)
    
    def get_display_name(self, gene_name, virus_family='flavivirus'):
        """Get display name for a gene"""
        standard_name = self.get_standard_name(gene_name, virus_family)
        
        if virus_family in self.standards and 'standard_genes' in self.standards[virus_family]:
            gene_info = self.standards[virus_family]['standard_genes'].get(standard_name, {})
            return gene_info.get('display_name', standard_name)
        
        return standard_name
    
    def get_gene_color(self, gene_name, virus_family='flavivirus'):
        """Get color for a gene"""
        standard_name = self.get_standard_name(gene_name, virus_family)
        
        if virus_family in self.standards and 'standard_genes' in self.standards[virus_family]:
            gene_info = self.standards[virus_family]['standard_genes'].get(standard_name, {})
            return gene_info.get('color', '#808080')
        
        return '#808080'
    
    def get_gene_type(self, gene_name, virus_family='flavivirus'):
        """Get gene type (structural/nonstructural)"""
        standard_name = self.get_standard_name(gene_name, virus_family)
        
        if virus_family in self.standards and 'standard_genes' in self.standards[virus_family]:
            gene_info = self.standards[virus_family]['standard_genes'].get(standard_name, {})
            return gene_info.get('type', 'unknown')
        
        return 'unknown'
    
    def get_all_mappings(self, virus_family='flavivirus'):
        """Get all gene mappings for a virus family"""
        mappings = {}
        
        if virus_family in self.standards and 'standard_genes' in self.standards[virus_family]:
            for gene_key, gene_data in self.standards[virus_family]['standard_genes'].items():
                mappings[gene_key] = {
                    'display_name': gene_data.get('display_name', gene_key),
                    'color': gene_data.get('color', '#808080'),
                    'type': gene_data.get('type', 'unknown'),
                    'aliases': gene_data.get('aliases', [])
                }
        
        return mappings
    
    def _get_default_standards(self):
        """Minimal fallback standards if file not found"""
        return {
            'flavivirus': {
                'standard_genes': {
                    'C': {'display_name': 'C', 'color': '#4575b4', 'type': 'structural', 'aliases': []},
                    'E': {'display_name': 'E', 'color': '#abd9e9', 'type': 'structural', 'aliases': []},
                    'NS1': {'display_name': 'NS1', 'color': '#fdae61', 'type': 'nonstructural', 'aliases': []},
                    'NS3': {'display_name': 'NS3', 'color': '#a50026', 'type': 'nonstructural', 'aliases': []},
                    'NS5': {'display_name': 'NS5', 'color': '#c2a5cf', 'type': 'nonstructural', 'aliases': []}
                }
            }
        }


# Example usage
if __name__ == "__main__":
    mapper = GeneNameMapper()
    
    # Test mappings
    test_names = [
        "anchored_capsid_protein_ancC",
        "RNA-dependent_RNA_polymerase_NS5",
        "envelope_protein_E",
        "nonstructural_protein_NS3"
    ]
    
    print("Gene Name Mapping Examples:")
    print("-" * 50)
    for name in test_names:
        standard = mapper.get_standard_name(name)
        display = mapper.get_display_name(name)
        color = mapper.get_gene_color(name)
        print(f"{name:40} → {display:6} (color: {color})")