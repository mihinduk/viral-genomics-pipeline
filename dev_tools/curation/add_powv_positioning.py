#\!/usr/bin/env python3
"""
Add POWV-specific gene label positioning
=======================================
Development tool for adding Powassan virus label positioning to visualization scripts.
Handles prM/M overlap for tick-borne flaviviruses.

Usage: python3 add_powv_positioning.py
"""

import re
from pathlib import Path

def add_powv_positioning():
    """Add POWV-specific positioning to visualization scripts"""
    
    script_dir = Path(__file__).parent.parent.parent
    
    # Update mutation visualization script
    mut_script = script_dir / 'viral_pipeline' / 'visualization' / 'visualize_mutations_python_outdir.py'
    depth_script = script_dir / 'viral_pipeline' / 'visualization' / 'visualize_depth.py'
    
    for script_path in [mut_script, depth_script]:
        if not script_path.exists():
            print(f"Warning: {script_path} not found")
            continue
            
        with open(script_path, 'r') as f:
            content = f.read()
        
        # Add POWV-specific check after the offset_genes definition
        if 'HM440560' not in content:
            # Find the location after offset_genes definition
            pattern = r'(offset_genes = \{[^}]+\})'
            replacement = r'\1\n    \n    # Add POWV-specific offsets for tick-borne flaviviruses\n    if accession and "HM440560" in accession:\n        offset_genes.update({\n            "membrane_glycoprotein_precursor_prM": {"offset_x": 0, "offset_y": -0.02, "fontsize": 8},\n            "membrane_glycoprotein_M": {"offset_x": 0, "offset_y": 0.02, "fontsize": 8}\n        })'
            
            content = re.sub(pattern, replacement, content, flags=re.DOTALL)
            
            with open(script_path, 'w') as f:
                f.write(content)
            
            print(f"✅ Added POWV positioning to {script_path.name}")
        else:
            print(f"ℹ️  POWV positioning already exists in {script_path.name}")

if __name__ == "__main__":
    add_powv_positioning()
