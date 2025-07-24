#\!/usr/bin/env python3
"""
Remove hardcoded POWV coordinates to use database
===============================================
Removes old hardcoded Powassan virus coordinates from visualization scripts
so they use the proper database configuration instead.
"""

def remove_hardcoded_powv():
    """Remove hardcoded POWV section from visualization scripts"""
    
    scripts = [
        'viral_pipeline/visualization/visualize_mutations_python_outdir.py',
        'viral_pipeline/visualization/visualize_depth.py'
    ]
    
    for script_path in scripts:
        with open(script_path, 'r') as f:
            lines = f.readlines()
        
        # Find and remove the hardcoded POWV section
        new_lines = []
        skip_section = False
        indent_level = 0
        
        for line in lines:
            # Start of POWV hardcoded section
            if "accession == 'HM440560.1' or accession == 'HM440560'" in line:
                skip_section = True
                indent_level = len(line) - len(line.lstrip())
                continue
                
            # End of POWV section (when we reach unindented code or different condition)
            if skip_section:
                line_indent = len(line) - len(line.lstrip())
                if line.strip() and line_indent <= indent_level and not line.strip().startswith('#'):
                    if 'elif' in line or 'else' in line or line_indent < indent_level:
                        skip_section = False
                        new_lines.append(line)
                    elif "print('Using hardcoded" in line:
                        continue  # Skip this line and end section
                    else:
                        continue
                elif not line.strip():  # Keep blank lines
                    new_lines.append(line)
                # Skip all other lines in the section
                continue
            else:
                new_lines.append(line)
        
        # Write back the cleaned file
        with open(script_path, 'w') as f:
            f.writelines(new_lines)
        
        print(f"✅ Removed hardcoded POWV from {script_path}")

if __name__ == "__main__":
    remove_hardcoded_powv()
