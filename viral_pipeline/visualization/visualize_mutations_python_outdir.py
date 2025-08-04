#!/usr/bin/env python3
"""
visualize_mutations_python.py
Python version of the viral mutation visualization tool.
Creates identical output to the R version but uses only standard Python packages.
Designed for M2 Mac compatibility and ease of use for bench scientists.
"""

import argparse
import requests
import json
import os
from datetime import datetime

def download_virus_configs():
    """Download the latest virus configuration file from GitHub for QC tracking"""
    config_url = "https://raw.githubusercontent.com/mihinduk/viral-genomics-pipeline/main/virus_configs/known_viruses.json"
    
    try:
        print("📥 Downloading latest virus configuration...")
        response = requests.get(config_url, timeout=30)
        response.raise_for_status()
        
        # Create virus_configs directory if it doesnt exist
        os.makedirs("virus_configs", exist_ok=True)
        
        # Save with timestamp for QC tracking
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        config_file = f"virus_configs/known_viruses_{timestamp}.json"
        
        with open(config_file, "w") as f:
            f.write(response.text)
        
        # Also save as the main config file
        with open("virus_configs/known_viruses.json", "w") as f:
            f.write(response.text)
            
        print(f"✅ Virus configuration downloaded and saved to {config_file}")
        return json.loads(response.text)
        
    except Exception as e:
        print(f"⚠️  Failed to download virus config: {e}")
        print("   Falling back to local configuration if available...")
        return None


# Import gene name mapper for standardized display names
try:
    from gene_name_mapper import GeneNameMapper
    gene_mapper = GeneNameMapper()
except ImportError:
    print("Warning: Gene name mapper not available, using full names")
    gene_mapper = None

# Load virus configurations with GitHub download
KNOWN_VIRUSES = {}
try:
    # First try to download the latest config
    downloaded_config = download_virus_configs()
    if downloaded_config:
        KNOWN_VIRUSES = downloaded_config
    else:
        # Fall back to local file
        try:
            with open('virus_configs/known_viruses.json', 'r') as f:
                KNOWN_VIRUSES = json.load(f)
            print('📁 Using local virus configuration')
        except:
            # Final fallback to root directory
            try:
                with open('known_viruses.json', 'r') as f:
                    KNOWN_VIRUSES = json.load(f)
                print('📁 Using root directory virus configuration')
            except:
                print('⚠️  No virus configuration found')
                KNOWN_VIRUSES = {}
except Exception as e:
    print(f'⚠️  Error loading virus configuration: {e}')
    KNOWN_VIRUSES = {}

# Define configuration functions using the downloaded config
def get_virus_name(accession):
    info = KNOWN_VIRUSES.get(accession, {})
    return info.get('name', f'Unknown virus ({accession})')

def get_gene_info(accession):
    info = KNOWN_VIRUSES.get(accession, {})
    gene_coords = info.get('gene_coords', {})
    colors = info.get('colors', {})
    
    # Apply gene name mapping if available
    if gene_mapper:
        family = info.get('family', 'flavivirus')
        display_names = {}
        mapped_colors = {}
        for gene in gene_coords:
            display_name = gene_mapper.get_display_name(gene, family)
            display_names[gene] = display_name
            # Get standardized color
            mapped_colors[gene] = gene_mapper.get_gene_color(gene, family)
    else:
        display_names = {gene: gene for gene in gene_coords}
        mapped_colors = colors
        
    return (gene_coords, 
            mapped_colors,
            info.get('structural_genes', []),
            info.get('nonstructural_genes', []),
            display_names)

print('Using downloaded virus configuration')
import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import warnings

# Suppress matplotlib warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)

def classify_mutation_type(effect_string, default_color="#1f77b4"):
    """Classify mutation type based on SnpEff EFFECT annotation
    Returns (color, category) tuple for comprehensive mutation classification
    """
    if not effect_string or pd.isna(effect_string):
        return default_color, "Other"
    
    effect = str(effect_string).lower()
    
    # High impact mutations (most severe)
    if any(term in effect for term in ['stop_gained', 'nonsense']):
        return "#d62728", "Nonsense/Stop Gained"  # Dark red
    elif 'stop_lost' in effect:
        return "#8b0000", "Stop Lost"  # Dark red
    elif any(term in effect for term in ['start_lost', 'initiator_codon']):
        return "#800080", "Start Lost"  # Purple
    elif 'frameshift' in effect:
        return "#ff1493", "Frameshift"  # Deep pink
    elif any(term in effect for term in ['splice_acceptor', 'splice_donor', 'splice_region']):
        return "#4169e1", "Splice Site"  # Royal blue
    elif any(term in effect for term in ['chromosome_large_deletion', 'structural_variant']):
        return "#ff4500", "Other High Impact"  # Orange red
    
    # Moderate impact mutations
    elif 'missense_variant' in effect:
        return "#ff7f0e", "Missense"  # Orange
    elif any(term in effect for term in ['inframe_deletion', 'inframe_insertion', 'disruptive_inframe']):
        return "#ffa500", "In-frame Indel"  # Orange
    
    # Low impact mutations
    elif 'synonymous_variant' in effect:
        return "#2ca02c", "Synonymous"  # Green
    
    # Modifier mutations (regulatory)
    elif any(term in effect for term in ['5_prime_utr', '3_prime_utr', 'upstream', 'downstream']):
        return "#32cd32", "Regulatory"  # Lime green
    elif 'intron' in effect:
        return "#90ee90", "Intronic"  # Light green
    elif 'intergenic' in effect:
        return "#98fb98", "Intergenic"  # Pale green
    
    # Default for unclassified
    return "#708090", "Other"  # Slate gray

def check_dependencies():
    """Check if all required Python packages are available"""
    required_packages = {
        'pandas': 'pandas',
        'numpy': 'numpy', 
        'matplotlib': 'matplotlib'
    }
    
    missing = []
    versions = {}
    
    for package_name, import_name in required_packages.items():
        try:
            module = __import__(import_name)
            if hasattr(module, '__version__'):
                versions[package_name] = module.__version__
            else:
                versions[package_name] = 'unknown'
        except ImportError:
            missing.append(package_name)
    
    if missing:
        print("❌ Missing required packages:")
        for pkg in missing:
            print(f"   - {pkg}")
        print("\nTo install missing packages:")
        print("   pip install " + " ".join(missing))
        print("\nOr if using conda:")
        print("   conda install " + " ".join(missing))
        return False
    
    print("✅ All required packages found:")
    for pkg, version in versions.items():
        print(f"   - {pkg}: {version}")
    
    return True

# Colors matching the R version exactly
# Genome length for West Nile virus
GENOME_LENGTH = 11029

def map_position_to_gene(position, accession):
    """Map a genomic position to the corresponding gene"""
    gene_info = get_gene_info(accession)
    if len(gene_info) == 5:
        gene_coords, _, _, _, _ = gene_info
    else:
        gene_coords, _, _, _ = gene_info
    for gene, (start, end) in gene_coords.items():
        if start <= position <= end:
            return gene
    return None

def parse_amino_acid_change(hgvs_p):
    """Parse amino acid change from HGVSp notation"""
    if pd.isna(hgvs_p) or hgvs_p == '':
        return None
    
    # Handle various HGVSp formats
    if 'p.' in str(hgvs_p):
        aa_change = str(hgvs_p).split('p.')[1]
        # Convert * to TER for stop codons
        aa_change = aa_change.replace('*', 'TER')
        return aa_change
    
    return None

def filter_mutations(df, cutoff, effect_filter=None):
    """Filter mutations by allele frequency and optionally by effect"""
    # Convert allele frequency to numeric if it's not already
    if 'Allele_Frequency' in df.columns:
        df['Allele_Frequency'] = pd.to_numeric(df['Allele_Frequency'], errors='coerce')
        df = df[df['Allele_Frequency'] >= cutoff]
    
    # Filter by effect type if specified
    if effect_filter:
        if 'EFFECT' in df.columns:
            df = df[df['EFFECT'].isin(effect_filter)]
    
    return df

def create_genome_diagram(ax, mutations_df, title, gene_filter="all", highlight_freq=0.5, accession=None):
    """Create the linear genome diagram with gene blocks"""
    # Get virus-specific gene information
    gene_info = get_gene_info(accession)
    if len(gene_info) == 5:
        gene_coords, gene_colors, _, _, display_names = gene_info
    else:
        gene_coords, gene_colors, _, _ = gene_info
        display_names = {gene: gene for gene in gene_coords}

    """Create the linear genome diagram with gene blocks"""
    # Get dynamic gene info for this virus
    gene_info = get_gene_info(accession)
    if len(gene_info) == 5:
        gene_coords, gene_colors, structural_genes, nonstructural_genes, display_names = gene_info
    else:
        # Fallback for old format
        gene_coords, gene_colors, structural_genes, nonstructural_genes = gene_info
        display_names = {gene: gene for gene in gene_coords}    
    # If we only have Polyprotein, use hardcoded coordinates for known viruses
    if len(gene_coords) == 1 and 'Polyprotein' in gene_coords:
        # POWV now uses database configuration (removed hardcoded coordinates)
        pass

        # All viruses now use database configuration
    
    
    # Set up the plot
    ax.set_xlim(0, GENOME_LENGTH)
    # Tighter spacing - bring position line closer to genome for better space utilization
    ax.set_ylim(-0.35, 1.0)
    
    # Draw gene blocks
    gene_height = 0.3
    gene_y = 0.35
    
    # Draw 5' UTR
    # Find first gene start
    first_gene_start = min(coord[0] for coord in gene_coords.values()) if gene_coords else 1
    if first_gene_start > 1:
        utr5_rect = Rectangle((0, gene_y), first_gene_start-1, gene_height, 
                         facecolor='lightgray', edgecolor='black', linewidth=0.5)
        ax.add_patch(utr5_rect)
    
    # Draw genes and handle overlapping labels
    # Define genes that need label offsetting for flaviviruses
    offset_genes = {
        # Also add standardized gene names
        "C": {"offset_x": 200, "offset_y": 0.1, "fontsize": 9},
        "prM": {"offset_x": 20, "offset_y": 0.1, "fontsize": 9},
        "2K": {"offset_x": 0, "offset_y": 0.1, "fontsize": 9},
        'anchored_capsid_protein_ancC': {'offset_x': -15, 'offset_y': -0.02, 'fontsize': 7},  # ancC above and left
        'capsid_protein_C': {'offset_x': 200, 'offset_y': 0.1, 'fontsize': 9},  # C moved further right to center of gray bar between ancC and pr
        'membrane_glycoprotein_M': {'offset_x': 0, 'offset_y': 0, 'fontsize': 9},  # M same size as E, horizontally aligned
        'protein_pr': {'offset_x': -20, 'offset_y': -0.02, 'fontsize': 6},  # pr above and left
        'membrane_glycoprotein_precursor_prM': {'offset_x': 20, 'offset_y': 0.1, 'fontsize': 6},  # prM below and right
        'protein_2K': {'offset_x': 0, 'offset_y': 0.1, 'fontsize': 7},  # 2K above
        'nonstructural_protein_NS4A': {'offset_x': 0, 'offset_y': 0.02, 'fontsize': 7}  # NS4A below
    }
    # Get gene properties for hatching
    gene_properties = KNOWN_VIRUSES.get(accession, {}).get("gene_properties", {})
    for gene, (start, end) in gene_coords.items():
        width = end - start + 1
    
        # Check if gene should have hatching
        props = gene_properties.get(gene, {})
        should_hatch = props.get("hatching", False)
        hatch_pattern = props.get("hatch_pattern", "...") if should_hatch else None
        
        rect = Rectangle((start-1, gene_y), width, gene_height,
                        facecolor=gene_colors.get(gene, "#808080"), edgecolor="black", linewidth=0.5,
                        hatch=hatch_pattern)
        ax.add_patch(rect)
        
        # Add gene label with offset handling
        gene_center = start + width/2 - 1
        # Use display name if available
        display_name = display_names.get(gene, gene)
        
        # Check if this gene needs offset
        if gene in offset_genes:
            x_offset = offset_genes[gene].get('offset_x', 0)
            y_offset = offset_genes[gene].get('offset_y', 0)
            font_size = offset_genes[gene].get('fontsize', 8)
            label_x = gene_center + x_offset
            label_y = gene_y + gene_height/2 + y_offset
        else:
            label_x = gene_center
            label_y = gene_y + gene_height/2
            font_size = 9
            
        ax.text(label_x, label_y, display_name, 
               ha='center', va='center', fontsize=font_size, fontweight='bold')
    
    # Draw 3' UTR
    # Find last gene end
    last_gene_end = max(coord[1] for coord in gene_coords.values()) if gene_coords else GENOME_LENGTH
    utr3_start = last_gene_end
    if utr3_start < GENOME_LENGTH:
        utr3_rect = Rectangle((utr3_start, gene_y), GENOME_LENGTH - utr3_start, gene_height,
                             facecolor='lightgray', edgecolor='black', linewidth=0.5)
        ax.add_patch(utr3_rect)
    
    # Add structural/non-structural labels
    # Find structural gene boundaries (use display names for classification)
    struct_coords = []
    nonstruct_coords = []
    
    for gene, (start, end) in gene_coords.items():
        display_name = display_names.get(gene, gene)
        # Check if the display name indicates structural protein
        if any(struct_gene in gene for struct_gene in structural_genes) or display_name in ['C', 'ancC', 'pr', 'prM', 'M', 'E', 'Env']:
            struct_coords.append((start, end))
        else:
            nonstruct_coords.append((start, end))
    
    if struct_coords:
        struct_start = min(coord[0] for coord in struct_coords) - 1
        struct_end = max(coord[1] for coord in struct_coords) - 1
    else:
        struct_start = 0
        struct_end = 0
    struct_center = (struct_start + struct_end) / 2
    
    # Find non-structural gene boundaries
    if nonstruct_coords:
        nonstruct_start = min(coord[0] for coord in nonstruct_coords) - 1
        nonstruct_end = max(coord[1] for coord in nonstruct_coords) - 1
    else:
        nonstruct_start = 0
        nonstruct_end = 0  
    nonstruct_center = (nonstruct_start + nonstruct_end) / 2
    
    # Structural proteins label and underline
    ax.text(struct_center, 0.8, 'Structural proteins', ha='center', va='center', 
           fontsize=12, fontweight='bold', color='#4575b4')
    ax.plot([struct_start, struct_end], [0.75, 0.75], color='#4575b4', linewidth=3)
    
    # Non-structural proteins label and underline  
    ax.text(nonstruct_center, 0.8, 'Non-structural proteins', ha='center', va='center',
           fontsize=12, fontweight='bold', color='#d73027')
    ax.plot([nonstruct_start, nonstruct_end], [0.75, 0.75], color='#d73027', linewidth=3)
    
    # Draw mutation lines
    if not mutations_df.empty:
        # Apply gene filtering to mutations for display
        display_mutations = filter_genes_for_display(mutations_df.copy(), gene_filter)
        
        # Group mutations by position to handle multiple mutations at same position
        position_mutations = {}
        for _, mutation in display_mutations.iterrows():
            pos = mutation['POS']
            if pos not in position_mutations:
                position_mutations[pos] = []
            position_mutations[pos].append(mutation)
        
        # Draw lines, prioritizing non-synonymous mutations at same position
        for pos, mutations in position_mutations.items():
            gene = map_position_to_gene(pos, accession)
            if gene:
                # Determine the highest priority mutation type at this position
                priority_scores = {
                    'stop_gained': 6, 'nonsense': 6, 'stop_lost': 5, 'start_lost': 4, 
                    'frameshift': 3, 'missense_variant': 2, 'synonymous_variant': 1
                }
                
                best_mutation = None
                best_priority = 0
                best_freq = 0
                
                for mutation in mutations:
                    mutation_type = mutation.get("EFFECT", "").lower()
                    # Calculate priority for this mutation
                    priority = 0
                    for effect, score in priority_scores.items():
                        if effect in mutation_type:
                            priority = max(priority, score)
                    
                    freq = mutation.get('Allele_Frequency', 0)
                    
                    # Choose mutation with highest priority, or highest frequency if same priority
                    if priority > best_priority or (priority == best_priority and freq > best_freq):
                        best_mutation = mutation
                        best_priority = priority
                        best_freq = freq
                
                if best_mutation is not None:
                    # Determine color based on best mutation type
                    mutation_type = best_mutation.get("EFFECT", "")
                    if "synonymous_variant" in mutation_type:
                        color = "#2ca02c"  # Green for synonymous
                    elif "missense_variant" in mutation_type:
                        color = "#ff7f0e"  # Orange for missense
                    elif "stop_gained" in mutation_type or "nonsense" in mutation_type:
                        color = "#d62728"  # Red for stop/nonsense
                    else:
                        color = gene_colors.get(gene, "#1f77b4")  # Blue default
                    
                    # Highlight high-frequency mutations with thicker lines
                    freq = best_mutation.get('Allele_Frequency', 0)
                    if freq >= highlight_freq:
                        linewidth = 3.0
                        alpha = 1.0
                    else:
                        linewidth = 1.5
                        alpha = 0.8
                    ax.plot([pos, pos], [0.25, 0.35], color=color, linewidth=linewidth, alpha=alpha)
    
    # Create dynamic legend based on mutation types found in data
    mutation_types_found = set()
    if not mutations_df.empty:
        display_mutations = filter_genes_for_display(mutations_df.copy(), gene_filter)
        for pos, mutations in position_mutations.items():
            for mutation in mutations:
                mutation_type = mutation.get("EFFECT", "")
                _, category = classify_mutation_type(mutation_type)
                mutation_types_found.add(category)
    
    # Build legend elements based on what's actually in the data
    legend_elements = []
    legend_mapping = {
        'Synonymous': '#2ca02c',
        'Missense': '#ff7f0e', 
        'Nonsense/Stop Gained': '#d62728',
        'Stop Lost': '#8b0000',
        'Start Lost': '#800080',
        'Frameshift': '#ff1493',
        'Splice Site': '#4169e1',
        'Other High Impact': '#ff4500',
        'In-frame Indel': '#ffa500',
        'Regulatory': '#32cd32',
        'Intronic': '#90ee90',
        'Intergenic': '#98fb98',
        'Other': '#708090'
    }
    
    for category in sorted(mutation_types_found):
        if category in legend_mapping:
            color = legend_mapping[category]
            legend_elements.append(plt.Line2D([0], [0], color=color, linewidth=2, label=category))
    
    if legend_elements:
        # Use 2 columns if more than 6 items, otherwise 1 column
        ncol = 2 if len(legend_elements) > 6 else 1
        ax.legend(handles=legend_elements, loc='upper right', frameon=True, 
                  fancybox=True, shadow=True, fontsize=9, ncol=ncol)
    
    # Format x-axis
    ax.set_xlabel('Genome Position', fontsize=12, fontweight='bold')
    ax.set_xticks(np.arange(0, GENOME_LENGTH+1, 1000))
    ax.set_xticklabels([f'{x:,}' for x in np.arange(0, GENOME_LENGTH+1, 1000)])
    
    # Remove y-axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add title
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

def filter_genes_for_display(mutations_df, gene_filter):
    """Filter mutations based on gene selection"""
    if gene_filter == 'all':
        return mutations_df
    elif gene_filter == 'structural':
        # Filter for structural genes
        mutations_df['Gene'] = mutations_df['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
        return mutations_df[mutations_df['Gene'].isin(structural_genes)]
    elif gene_filter == 'non-structural':
        # Filter for non-structural genes  
        mutations_df['Gene'] = mutations_df['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
        return mutations_df[mutations_df['Gene'].isin(nonstructural_genes)]
    else:
        # Comma-separated gene list
        selected_genes = [g.strip() for g in gene_filter.split(',')]
        mutations_df['Gene'] = mutations_df['POS'].apply(lambda pos: map_position_to_gene(pos, accession))
        return mutations_df[mutations_df['Gene'].isin(selected_genes)]

def create_mutation_tables(fig, mutations_df, start_row=0.4, gene_filter="all", accession=None):
    """Create tables showing ALL mutations for each gene with complete parity"""
    
    # Show ALL mutations (synonymous + non-synonymous) for complete parity with lines
    # Include all mutation types to match the vertical lines displayed
    all_mutations = mutations_df.copy()
    
    # Apply gene filtering
    all_mutations = filter_genes_for_display(all_mutations, gene_filter)
    
    if all_mutations.empty:
        # Add text saying no mutations found
        fig.text(0.5, 0.25, 'No mutations found above the specified cutoff',
                ha='center', va='center', fontsize=14, style='italic')
        return
    
    # Add gene column based on position
    # Add gene column - use GENE_NAME if available, otherwise map by position
    if "GENE_NAME" in all_mutations.columns:
        # Use existing gene names from SnpEff output
        all_mutations["Gene"] = all_mutations["GENE_NAME"]
        
        # Apply gene name mapping for known aliases
        gene_mapping = KNOWN_VIRUSES.get(accession, {}).get("snpeff_gene_mapping", {})
        alias_mapping = {
            "Env": "E",  # Common alias for envelope protein
            "env": "E",
            "envelope": "E"
        }
        # Combine both mappings
        all_mapping = {**gene_mapping, **alias_mapping}
        all_mutations["Gene"] = all_mutations["Gene"].map(lambda x: all_mapping.get(x, x))
    else:
        # Fallback to position-based mapping
        all_mutations["Gene"] = all_mutations["POS"].apply(lambda pos: map_position_to_gene(pos, accession))
    all_mutations = all_mutations.dropna(subset=['Gene'])
    
    # Parse amino acid changes
    all_mutations['AA_Change'] = all_mutations['HGVSp'].apply(parse_amino_acid_change)
    
    # Group by gene
    # Get dynamic gene list for this virus
    gene_info = get_gene_info(accession)
    if len(gene_info) == 5:
        gene_coords, gene_colors, _, _, display_names = gene_info
    else:
        gene_coords, gene_colors, _, _ = gene_info
        display_names = {gene: gene for gene in gene_coords}
    genes_with_mutations = []
    for gene in gene_coords.keys():
        gene_mutations = all_mutations[all_mutations['Gene'] == gene]
        if not gene_mutations.empty:
            genes_with_mutations.append((gene, gene_mutations))
    
    if not genes_with_mutations:
        fig.text(0.5, 0.25, 'No mutations found above the specified cutoff',
                ha='center', va='center', fontsize=14, style='italic')
        return
    
    # Smart table positioning with collision avoidance and row alignment
    table_width = 0.18
    table_height = 0.15
    table_spacing = 0.02
    
    # Get gene positions for spatial awareness
    gene_positions = {}
    for gene, mutations in genes_with_mutations:
        if gene in gene_coords:
            start, end = gene_coords[gene]
            center = (start + end) / 2
            # Get genome length from KNOWN_VIRUSES config
            genome_length = KNOWN_VIRUSES.get(accession, {}).get("genome_length", 10000)
            gene_positions[gene] = center / genome_length  # Normalize to 0-1
    
    # Sort genes by their genomic position for logical layout
    sorted_genes = sorted(genes_with_mutations, 
                         key=lambda x: gene_positions.get(x[0], 0))
    
    # Calculate grid layout with collision avoidance
    max_cols = 5
    first_row_positions = []  # Track occupied positions in first row
    second_row_positions = []  # Available positions for second row
    
    # Position first row - evenly distributed
    first_row_count = min(len(sorted_genes), max_cols)
    available_width = 1.0 - 2 * table_spacing
    total_table_width = first_row_count * table_width
    remaining_space = available_width - total_table_width
    spacing = remaining_space / (first_row_count + 1)
    
    # Calculate first row positions
    for i in range(first_row_count):
        x = table_spacing + spacing + i * (table_width + spacing)
        first_row_positions.append(x)
    
    # For second row, identify collision zones and find safe positions
    if len(sorted_genes) > max_cols:
        # Calculate potential collision zones (where first row table labels extend down)
        collision_zones = []
        for x in first_row_positions:
            # Each table label extends slightly below the table
            label_left = x - 0.01  # Small buffer
            label_right = x + table_width + 0.01
            collision_zones.append((label_left, label_right))
        
        # Generate safe positions for second row
        # Try to place tables in gaps between collision zones
        safe_positions = []
        
        # Start from left edge
        current_x = table_spacing + spacing
        
        while current_x + table_width <= 1.0 - table_spacing and len(safe_positions) < len(sorted_genes) - first_row_count:
            # Check if this position would cause collision
            table_left = current_x
            table_right = current_x + table_width
            
            collision = False
            for zone_left, zone_right in collision_zones:
                # Check if table or its label would overlap with collision zone
                if not (table_right < zone_left or table_left > zone_right):
                    collision = True
                    break
            
            if not collision:
                safe_positions.append(current_x)
                current_x += table_width + spacing
            else:
                # Move past the collision zone
                for zone_left, zone_right in collision_zones:
                    if current_x < zone_right:
                        current_x = zone_right + table_spacing
                        break
        
        second_row_positions = safe_positions
    
    # Apply positions to genes
    for i, (gene, gene_mutations) in enumerate(sorted_genes):
        if i < first_row_count:
            # First row - aligned at same Y
            x = first_row_positions[i]
            y = start_row - 0.05  # Top row alignment
        else:
            # Second row - use safe positions
            second_row_index = i - first_row_count
            if second_row_index < len(second_row_positions):
                x = second_row_positions[second_row_index]
            else:
                # Fallback to regular grid if we run out of safe positions
                col = second_row_index % max_cols
                x = table_spacing + spacing + col * (table_width + spacing)
            y = start_row - 0.25  # Second row
        
        # Create table data
        table_data = []
        for _, mut in gene_mutations.iterrows():
            nt_change = f"{mut['POS']}{mut['REF']}>{mut['ALT']}"
            aa_change = mut['AA_Change'] if mut['AA_Change'] else 'Unknown'
            table_data.append([nt_change, aa_change])
        
        # Create table
        if table_data:
            # Create subplot for this table
            table_ax = fig.add_axes([x, y-0.15, table_width, 0.15])
            table_ax.axis('off')
            
            # Create table
            table = table_ax.table(cellText=table_data,
                                 colLabels=['Position', 'AA Change'],
                                 cellLoc='center',
                                 loc='center',
                                 colWidths=[0.5, 0.5])
            
            # Style the table
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
            
            # Color the header
            for i in range(2):
                table[(0, i)].set_facecolor(gene_colors.get(gene, "#808080"))
                table[(0, i)].set_text_props(weight='bold', color='white')
            
            # Add gene name above table
            # Use display name for table header
            display_name = display_names.get(gene, gene)
            fig.text(x + table_width/2, y + 0.02, display_name, 
                    ha='center', va='center', fontsize=12, 
                    fontweight='bold', color=gene_colors.get(gene, "#808080"))

def save_mutations_table(mutations_df, output_path, cutoff, accession=None):
    """Save all mutations to a TSV file"""
    if mutations_df.empty:
        return
    
    # Prepare output data
    output_data = []
    for _, mut in mutations_df.iterrows():
        gene = map_position_to_gene(mut['POS'], accession)
        aa_change = parse_amino_acid_change(mut.get('HGVSp', ''))
        
        freq_percent = f"{mut['Allele_Frequency']*100:.1f}%" if 'Allele_Frequency' in mut else "Unknown"
        
        output_data.append({
            'Position': mut['POS'],
            'Gene': gene if gene else 'Unknown',
            'Nucleotide_Change': f"{mut['REF']}>{mut['ALT']}",
            'Amino_Acid_Change': aa_change if aa_change else 'N/A',
            'Effect': mut.get('EFFECT', 'Unknown'),
            'Impact': mut.get('PUTATIVE_IMPACT', 'Unknown'),
            'Frequency_Percent': freq_percent
        })
    
    # Create DataFrame and save
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_path, sep='\t', index=False)
    print(f"✅ Mutations table saved to: {output_path}")


def generate_html_report(mutations_df, image_path, html_path, accession, cutoff):
    """Generate HTML report with embedded image and mutation data"""
    import os
    
    # Get virus info for the report
    virus_name = get_virus_name(accession)
    
    # Count mutations
    total_mutations = len(mutations_df)
    
    # Create HTML content
    html_content = f"""<\!DOCTYPE html>
<html>
<head>
    <title>Viral Mutation Report - {virus_name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
        .image-container {{ text-align: center; margin: 20px 0; }}
        .image-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        .stats {{ display: flex; justify-content: space-around; margin: 20px 0; }}
        .stat-box {{ background-color: #e9ecef; padding: 15px; border-radius: 5px; text-align: center; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Viral Mutation Analysis Report</h1>
        <h2>{virus_name} ({accession})</h2>
        <p>Generated on: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="stats">
        <div class="stat-box">
            <h3>Total Mutations</h3>
            <p style="font-size: 24px; font-weight: bold;">{total_mutations}</p>
        </div>
        <div class="stat-box">
            <h3>Cutoff Threshold</h3>
            <p style="font-size: 24px; font-weight: bold;">{cutoff*100:.1f}%</p>
        </div>
    </div>
    
    <div class="image-container">
        <h3>Genome Visualization</h3>
        <img src="{os.path.basename(image_path)}" alt="Viral mutation visualization">
    </div>
    
    <div style="margin-top: 40px; font-size: 12px; color: #666;">
        <p>Generated with Viral Mutation Visualizer  < /dev/null |  Family-based configuration system</p>
    </div>
</body>
</html>"""
    
    # Write HTML file
    with open(html_path, 'w') as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(
        description='Visualize viral genome mutations (Python version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --input mutations.tsv --output plot.png --cutoff 0.05
  %(prog)s --input data.tsv --output figure.pdf --cutoff 0.9 --mutation-genes structural
  %(prog)s --input data.tsv --output plot.png --cutoff 0.05 --mutation-genes "NS3,NS5" --width 20 --height 12
  %(prog)s --input data.tsv --output plot.pdf --cutoff 0.9 --mutation-genes non-structural --highlight-freq 0.8
        """
    )
    
    parser.add_argument('--input', required=True,
                       help='Input TSV file with SnpEff-annotated LoFreq output')
    parser.add_argument('--output', required=True, 
                       help='Output file path (PNG or PDF)')
    parser.add_argument('--cutoff', type=float, default=0.05,
                       help='Allele frequency cutoff (default: 0.05)')
    parser.add_argument('--mutation-genes', default='all',
                       help='Genes to display mutations for: "all", "structural", "non-structural", or comma-separated gene names (default: all)')
    parser.add_argument('--width', type=float, default=16,
                       help='Figure width in inches (default: 16)')
    parser.add_argument('--height', type=float, default=10,
                       help='Figure height in inches (default: 10)')
    parser.add_argument('--highlight-freq', type=float, default=0.5,
                       help='Highlight mutations with frequency above this threshold (default: 0.5)')
    parser.add_argument('--accession', 
                       help='Manually specify GenBank accession (overrides accession in CHROM column)')
    parser.add_argument("--outdir", help="Output directory (default: current directory)")
    
    args = parser.parse_args()
    
    print("Viral Mutation Visualizer (Python Version)")
    print("=" * 50)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Check input file
    if not Path(args.input).exists():
        print(f"❌ Input file not found: {args.input}")
        sys.exit(1)
    
    print(f"📁 Input file: {args.input}")
    print(f"🎯 Output file: {args.output}")
    print(f"📊 Allele frequency cutoff: {args.cutoff}")
    
    try:
        # Read input data
        print("\n📖 Reading input data...")
        df = pd.read_csv(args.input, sep='\t')
        print(f"   Found {len(df)} total mutations")
        
        # Required columns
        required_cols = ['POS', 'REF', 'ALT']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"❌ Missing required columns: {missing_cols}")
            sys.exit(1)
        
        # Filter mutations by cutoff
        filtered_df = filter_mutations(df, args.cutoff)
        print(f"   {len(filtered_df)} mutations above cutoff {args.cutoff}")
        
        # Create the visualization
        print("\n🎨 Creating visualization...")
        fig = plt.figure(figsize=(args.width, args.height))
        
        # Main genome diagram (top 60% of figure)
        genome_ax = fig.add_axes([0.1, 0.5, 0.8, 0.4])
        
        # Get accession from command line or CHROM column
        if args.accession:
            accession = args.accession
        else:
            accession = df['CHROM'].iloc[0] if 'CHROM' in df.columns and not df.empty else "Unknown"
        
        # Create title based on gene filter
        gene_desc = ""
        if args.mutation_genes == 'structural':
            gene_desc = " (structural proteins)"
        elif args.mutation_genes == 'non-structural':
            gene_desc = " (non-structural proteins)"
        elif args.mutation_genes != 'all':
            gene_desc = f" ({args.mutation_genes})"
        
        virus_name = get_virus_name(accession)
        title = f"Mutations in {accession} - {virus_name}, complete genome{gene_desc} (cutoff: {args.cutoff*100:.1f}%)"
        
        create_genome_diagram(genome_ax, filtered_df, title, args.mutation_genes, args.highlight_freq, accession)
        
        # Mutation tables (bottom 40% of figure)
        create_mutation_tables(fig, filtered_df, start_row=0.4, gene_filter=args.mutation_genes, accession=accession)
        
        # Save the figure
        print(f"\n💾 Saving figure...")
        fig.savefig(args.output, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        
        # Handle output directory  
        if args.outdir:
            # Create output directory if specified
            outdir = Path(args.outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            
            # Use output filename but place in specified directory
            output_filename = Path(args.output).name
            final_output_path = outdir / output_filename
            
            # Table path in the same directory
            output_base = Path(args.output).stem
            table_path = outdir / f"{output_base}_mutations_table.tsv"
            
            print(f"\\n📁 Using output directory: {outdir}")
            print(f"   Full path: {outdir.absolute()}")
        else:
            # Use original output path
            final_output_path = Path(args.output)
            output_base = final_output_path.stem
            output_dir = final_output_path.parent
            table_path = output_dir / f"{output_base}_mutations_table.tsv"
        
        # Save the figure
        print(f"\\n💾 Saving figure...")
        fig.savefig(final_output_path, dpi=300, bbox_inches="tight", 
                   facecolor="white", edgecolor="none")
        print(f"✅ Figure saved to: {final_output_path}")
        
        # Save mutations table
        save_mutations_table(filtered_df, table_path, args.cutoff, accession)
        
        plt.close()
        
        print("\\n🎉 Visualization complete\!")
        print(f"   Main figure: {final_output_path}")
        print(f"   Data table: {table_path}")
        
    except Exception as e:
        print(f"\\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
