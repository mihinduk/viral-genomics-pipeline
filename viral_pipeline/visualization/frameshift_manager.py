#\!/usr/bin/env python3
"""
Frameshift Gene Manager for Viral Genomics Pipeline
Handles frameshift genes like NS1' that are produced by ribosomal frameshifting
"""

import json
from pathlib import Path

class FrameshiftGeneManager:
    def __init__(self, config_file=None):
        """Initialize frameshift gene manager"""
        if config_file is None:
            config_file = Path(__file__).parent / "frameshift_genes.json"
        
        self.config_file = config_file
        self.config = self._load_config()
    
    def _load_config(self):
        """Load frameshift gene configuration"""
        if self.config_file.exists():
            with open(self.config_file, 'r') as f:
                return json.load(f)
        else:
            return {"frameshift_genes": {}, "gene_types": {}}
    
    def is_frameshift_gene(self, gene_name, accession):
        """Check if a gene is a frameshift gene"""
        accession_config = self.config.get("frameshift_genes", {}).get(accession, {})
        return gene_name in accession_config
    
    def get_frameshift_info(self, gene_name, accession):
        """Get frameshift gene information"""
        accession_config = self.config.get("frameshift_genes", {}).get(accession, {})
        return accession_config.get(gene_name, {})
    
    def calculate_frameshift_coordinates(self, gene_name, accession, parent_coords):
        """Calculate frameshift gene coordinates based on parent gene"""
        info = self.get_frameshift_info(gene_name, accession)
        if not info:
            return None
        
        start_rule = info.get("start", "same_as_parent")
        end_extension = info.get("end_extension", 0)
        
        if start_rule == "same_as_parent":
            start = parent_coords[0]
        else:
            start = int(start_rule)
        
        end = parent_coords[1] + end_extension
        
        return [start, end]
    
    def get_visualization_properties(self, gene_name, accession):
        """Get visualization properties for frameshift genes"""
        info = self.get_frameshift_info(gene_name, accession)
        if not info:
            return {}
        
        gene_type_config = self.config.get("gene_types", {}).get("frameshift", {})
        
        return {
            "pattern": info.get("visualization", "hatched"),
            "transparency": gene_type_config.get("transparency", 0.7),
            "color_inherit": gene_type_config.get("color_inherit", True),
            "legend_suffix": gene_type_config.get("legend_suffix", " (frameshift)")
        }
    
    def get_parent_gene(self, gene_name, accession):
        """Get the parent gene for a frameshift gene"""
        info = self.get_frameshift_info(gene_name, accession)
        return info.get("parent_gene")
    
    def validate_frameshift_genes(self, virus_data):
        """Validate that frameshift genes are properly configured"""
        errors = []
        
        for accession, virus_info in virus_data.items():
            if accession not in self.config.get("frameshift_genes", {}):
                continue
            
            frameshift_config = self.config["frameshift_genes"][accession]
            gene_coords = virus_info.get("gene_coords", {})
            
            for frameshift_gene, config in frameshift_config.items():
                parent_gene = config.get("parent_gene")
                
                # Check if parent gene exists
                if parent_gene not in gene_coords:
                    errors.append(f"{accession}: Parent gene '{parent_gene}' not found for frameshift gene '{frameshift_gene}'")
                
                # Check if frameshift gene exists in coordinates
                if frameshift_gene not in gene_coords:
                    errors.append(f"{accession}: Frameshift gene '{frameshift_gene}' not found in gene_coords")
                else:
                    # Validate coordinates are calculated correctly
                    parent_coords = gene_coords[parent_gene]
                    expected_coords = self.calculate_frameshift_coordinates(frameshift_gene, accession, parent_coords)
                    actual_coords = gene_coords[frameshift_gene]
                    
                    if expected_coords != actual_coords:
                        errors.append(f"{accession}: {frameshift_gene} coordinates {actual_coords} don't match expected {expected_coords}")
        
        return errors

# Global instance for easy import
frameshift_manager = FrameshiftGeneManager()
