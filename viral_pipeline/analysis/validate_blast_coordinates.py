#\!/usr/bin/env python3
"""
Validate BLAST-mapped coordinates against our JSON configuration
Helps ensure coordinate conservation between sequenced and reference strains
"""

import argparse
import json
import sys
from pathlib import Path

def compare_coordinates(blast_coords, json_coords, tolerance=3):
    """Compare BLAST vs JSON coordinates with small tolerance for edge effects"""
    validation_results = {
        'matches': [],
        'mismatches': [],
        'missing_in_blast': [],
        'missing_in_json': []
    }
    
    all_genes = set(blast_coords.keys()) | set(json_coords.keys())
    
    for gene in sorted(all_genes):
        if gene not in blast_coords:
            validation_results['missing_in_blast'].append(gene)
        elif gene not in json_coords:
            validation_results['missing_in_json'].append(gene)
        else:
            blast_start, blast_end = blast_coords[gene]
            json_start, json_end = json_coords[gene]
            
            # Allow small tolerance for coordinate differences
            if (abs(blast_start - json_start) <= tolerance and 
                abs(blast_end - json_end) <= tolerance):
                validation_results['matches'].append({
                    'gene': gene,
                    'blast': [blast_start, blast_end],
                    'json': [json_start, json_end],
                    'exact_match': blast_start == json_start and blast_end == json_end
                })
            else:
                validation_results['mismatches'].append({
                    'gene': gene,
                    'blast': [blast_start, blast_end],
                    'json': [json_start, json_end],
                    'start_diff': blast_start - json_start,
                    'end_diff': blast_end - json_end
                })
    
    return validation_results

def generate_validation_report(results, query_acc, ref_acc):
    """Generate human-readable validation report"""
    print(f"\n🧬 Coordinate Validation Report")
    print(f"{'='*60}")
    print(f"Query: {query_acc} → Reference: {ref_acc}")
    print(f"{'='*60}\n")
    
    # Summary statistics
    total_genes = (len(results['matches']) + len(results['mismatches']) + 
                  len(results['missing_in_blast']) + len(results['missing_in_json']))
    
    print(f"Summary:")
    print(f"  Total genes checked: {total_genes}")
    print(f"  ✅ Matches: {len(results['matches'])}")
    print(f"  ❌ Mismatches: {len(results['mismatches'])}")
    print(f"  ⚠️  Missing in BLAST: {len(results['missing_in_blast'])}")
    print(f"  ⚠️  Missing in JSON: {len(results['missing_in_json'])}")
    
    # Validation passed?
    is_valid = len(results['mismatches']) == 0 and len(results['missing_in_blast']) == 0
    
    if is_valid:
        print(f"\n✅ VALIDATION PASSED: Coordinates are conserved\!")
    else:
        print(f"\n❌ VALIDATION FAILED: Manual review required")
    
    return is_valid

# Placeholder for actual implementation
print("Coordinate validation script created\!")
print("This script will validate BLAST results against JSON coordinates")
print("TODO: Implement parse_blast_annotation_output() based on your BLAST output format")
