#!/usr/bin/env python3
"""
Simple virus configuration manager without BioPython dependency.
Uses the known_viruses.json file for virus information.
"""
import json
import os

class SimpleVirusConfigManager:
    def __init__(self, config_dir=None):
        if config_dir is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            config_dir = script_dir
        
        self.config_file = os.path.join(config_dir, 'known_viruses.json')
        self.virus_data = self.load_virus_data()
    
    def load_virus_data(self):
        """Load virus data from JSON file"""
        try:
            with open(self.config_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Could not find {self.config_file}")
            return {}
        except json.JSONDecodeError:
            print(f"Warning: Invalid JSON in {self.config_file}")
            return {}
    
    def get_virus_info(self, accession):
        """Get virus information for a given accession"""
        # Clean accession (remove version number if present)
        clean_accession = accession.split('.')[0] + '.1' if '.' not in accession else accession
        
        if clean_accession in self.virus_data:
            return self.virus_data[clean_accession]
        
        # Fallback for unknown viruses
        return {
            'name': f'Unknown virus ({accession})',
            'family': 'unknown',
            'genome_length': 11000,
            'gene_coords': {'Polyprotein': [1, 10000]},
            'colors': {'Polyprotein': '#808080'},
            'structural_genes': [],
            'nonstructural_genes': ['Polyprotein']
        }

# Create global instance
simple_virus_manager = SimpleVirusConfigManager()