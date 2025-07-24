#\!/usr/bin/env python3
"""
Family Template Manager
======================
Manages family-based templates for virus visualization configuration.
Provides lookup and application of templates for new viruses.

Usage:
  python3 template_manager.py --list-templates
  python3 template_manager.py --find-template ZIKV
  python3 template_manager.py --apply-template KU955591.1 Flaviviridae_mosquito_borne
"""

import json
import argparse
from pathlib import Path

class TemplateManager:
    def __init__(self):
        self.templates_file = Path(__file__).parent.parent.parent / 'viral_pipeline' / 'visualization' / 'virus_configs' / 'family_templates.json'
        self.known_viruses_file = Path(__file__).parent.parent.parent / 'viral_pipeline' / 'visualization' / 'virus_configs' / 'known_viruses.json'
        
        with open(self.templates_file, 'r') as f:
            self.templates_data = json.load(f)
    
    def list_templates(self):
        """List all available templates"""
        templates = self.templates_data['templates']
        
        print("\n🧬 Available Family Templates:\n")
        for template_id, template in templates.items():
            print(f"📋 {template_id}")
            print(f"   Description: {template['description']}")
            print(f"   Gene Count: {template['gene_count']}")
            print(f"   Reference: {template['reference_accession']} ({template['reference_name']})")
            print(f"   Applies to: {', '.join(template['applied_to'])}")
            if 'key_differences' in template:
                print(f"   Key differences: {template['key_differences']}")
            print()
    
    def find_template_for_virus(self, virus_name):
        """Find appropriate template for a virus"""
        templates = self.templates_data['templates']
        
        matches = []
        for template_id, template in templates.items():
            if virus_name.upper() in [v.upper() for v in template['applied_to']]:
                matches.append((template_id, template))
        
        if matches:
            print(f"\n🎯 Templates for {virus_name}:\n")
            for template_id, template in matches:
                print(f"✅ {template_id}")
                print(f"   Reference: {template['reference_accession']}")
                print(f"   Gene Count: {template['gene_count']}")
                print()
        else:
            print(f"\n❌ No specific template found for {virus_name}")
            print("\n🔍 Suggested approach:")
            print("1. Identify virus family (Flaviviridae, Togaviridae, etc.)")
            print("2. Check if mosquito-borne or tick-borne (for flaviviruses)")
            print("3. Use closest template as starting point")
            print("4. Run BLAST annotation to get exact coordinates")
    
    def get_template(self, template_id):
        """Get specific template by ID"""
        templates = self.templates_data['templates']
        
        if template_id in templates:
            return templates[template_id]
        else:
            available = list(templates.keys())
            raise ValueError(f"Template '{template_id}' not found. Available: {available}")
    
    def create_virus_from_template(self, accession, template_id, custom_coords=None):
        """Create a new virus entry based on template"""
        template = self.get_template(template_id)
        
        # Load existing known viruses
        with open(self.known_viruses_file, 'r') as f:
            known_viruses = json.load(f)
        
        # Get reference virus configuration  
        ref_accession = template['reference_accession']
        if ref_accession not in known_viruses:
            raise ValueError(f"Reference virus {ref_accession} not found in known_viruses.json")
        
        ref_config = known_viruses[ref_accession]
        
        # Create new virus configuration based on template
        new_config = {
            "name": f"New virus (template: {template_id})",
            "family": ref_config['family'],
            "gene_coords": custom_coords if custom_coords else ref_config['gene_coords'].copy(),
            "colors": ref_config['colors'].copy(),
            "structural_genes": template['typical_genes'][:len([g for g in template['typical_genes'] if 'structural' in self._classify_gene(g)])],
            "nonstructural_genes": [g for g in template['typical_genes'] if 'nonstructural' in g or 'nsP' in g],
            "template_used": template_id,
            "needs_curation": True
        }
        
        return new_config
    
    def _classify_gene(self, gene_name):
        """Classify gene as structural or nonstructural"""
        structural_keywords = ['capsid', 'membrane', 'envelope', 'protein_pr', 'protein_6K']
        nonstructural_keywords = ['nonstructural', 'nsP', 'polymerase']
        
        for keyword in structural_keywords:
            if keyword in gene_name:
                return 'structural'
        
        for keyword in nonstructural_keywords:
            if keyword in gene_name:
                return 'nonstructural'
        
        return 'unknown'

def main():
    parser = argparse.ArgumentParser(description='Manage family templates for virus visualization')
    parser.add_argument('--list-templates', action='store_true', help='List all available templates')
    parser.add_argument('--find-template', help='Find template for specific virus (e.g., ZIKV, POWV)')
    parser.add_argument('--get-template', help='Get details for specific template ID')
    
    args = parser.parse_args()
    
    manager = TemplateManager()
    
    if args.list_templates:
        manager.list_templates()
    elif args.find_template:
        manager.find_template_for_virus(args.find_template)
    elif args.get_template:
        template = manager.get_template(args.get_template)
        print(json.dumps(template, indent=2))
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
