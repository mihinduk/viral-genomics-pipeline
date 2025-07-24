#\!/usr/bin/env python3
"""
POWV Database Entry Update Script
=================================
Development tool for updating Powassan virus entry in known_viruses.json
Based on literature analysis showing POWV has 10 proteins (not 13 like ZIKV).

Usage: python3 update_powv_entry.py
Note: This is a one-time curation script for developers only.
"""

import json
import os
from pathlib import Path

def update_powv_entry():
    """Update POWV entry with correct 10-protein structure"""
    
    # Find the known_viruses.json file
    script_dir = Path(__file__).parent.parent.parent
    config_file = script_dir / 'viral_pipeline' / 'visualization' / 'virus_configs' / 'known_viruses.json'
    
    if not config_file.exists():
        print(f"Error: Configuration file not found: {config_file}")
        return False
    
    # Load current configuration
    with open(config_file, 'r') as f:
        known_viruses = json.load(f)
    
    # Correct POWV entry based on literature
    powv_entry = {
        "name": "Powassan virus",
        "family": "flavivirus", 
        "gene_coords": {
            "capsid_protein_C": [110, 424],
            "membrane_glycoprotein_precursor_prM": [425, 919], 
            "membrane_glycoprotein_M": [670, 919],
            "envelope_protein_E": [920, 2404],
            "nonstructural_protein_NS1": [2405, 3460],
            "nonstructural_protein_NS2A": [3461, 4114],
            "nonstructural_protein_NS2B": [4115, 4231],
            "nonstructural_protein_NS3": [4232, 6085],
            "nonstructural_protein_NS4A": [6086, 6466],
            "nonstructural_protein_NS4B": [6467, 7210],
            "RNA-dependent_RNA_polymerase_NS5": [7211, 10111]
        },
        "colors": {
            "capsid_protein_C": "#4575b4",
            "membrane_glycoprotein_precursor_prM": "#74add1", 
            "membrane_glycoprotein_M": "#abd9e9",
            "envelope_protein_E": "#e0f3f8",
            "nonstructural_protein_NS1": "#fee090",
            "nonstructural_protein_NS2A": "#fdae61",
            "nonstructural_protein_NS2B": "#f46d43",
            "nonstructural_protein_NS3": "#d73027",
            "nonstructural_protein_NS4A": "#a50026",
            "nonstructural_protein_NS4B": "#762a83",
            "RNA-dependent_RNA_polymerase_NS5": "#5aae61"
        },
        "structural_genes": ["capsid_protein_C", "membrane_glycoprotein_precursor_prM", "membrane_glycoprotein_M", "envelope_protein_E"],
        "nonstructural_genes": ["nonstructural_protein_NS1", "nonstructural_protein_NS2A", "nonstructural_protein_NS2B", "nonstructural_protein_NS3", "nonstructural_protein_NS4A", "nonstructural_protein_NS4B", "RNA-dependent_RNA_polymerase_NS5"]
    }
    
    known_viruses["HM440560.1"] = powv_entry
    
    # Save updated file
    with open(config_file, 'w') as f:
        json.dump(known_viruses, f, indent=2)
    
    print(f"✅ Updated POWV entry in {config_file}")
    print("📊 POWV now configured with correct 10-protein tick-borne structure")
    return True

if __name__ == "__main__":
    update_powv_entry()
