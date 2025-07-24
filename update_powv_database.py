import json

# Update reference genomes database with correct POWV information
with open('virus_configs/reference_genomes_database.json', 'r') as f:
    db = json.load(f)

# Update POWV entry in Flaviviridae_tickborne family
powv_entry = {
    "accession": "HM440560.1",
    "virus_name": "Powassan virus",
    "strain": "LB",
    "genome_length": 10839,
    "description": "Tick-borne flavivirus with 10-protein polyprotein (no ancC, pr, 2K)",
    "year_isolated": 1958,
    "host": "Ixodes ticks, mammals",
    "curator_notes": "Tick-adapted flavivirus with simplified protein repertoire compared to mosquito-borne relatives",
    "visualization_settings": {
        "figure_width": 18,
        "figure_height": 12,
        "gene_height": 0.3,
        "offset_genes": {
            "membrane_glycoprotein_precursor_prM": {"offset_x": 0, "offset_y": -0.02, "fontsize": 8},
            "membrane_glycoprotein_M": {"offset_x": 0, "offset_y": 0.02, "fontsize": 8}
        }
    },
    "gene_coordinates": {
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
    "structural_genes": [
        "capsid_protein_C",
        "membrane_glycoprotein_precursor_prM",
        "membrane_glycoprotein_M",
        "envelope_protein_E"
    ],
    "nonstructural_genes": [
        "nonstructural_protein_NS1",
        "nonstructural_protein_NS2A", 
        "nonstructural_protein_NS2B",
        "nonstructural_protein_NS3",
        "nonstructural_protein_NS4A",
        "nonstructural_protein_NS4B",
        "RNA-dependent_RNA_polymerase_NS5"
    ],
    "quality_score": 8,
    "status": "curated",
    "last_updated": "2025-01-24",
    "source": "literature_based",
    "reference_publication": "Tick-adapted flavivirus genome architecture analysis"
}

# Update the database
db["families"]["Flaviviridae_tickborne"]["references"]["HM440560.1"] = powv_entry

# Save updated database
with open('virus_configs/reference_genomes_database.json', 'w') as f:
    json.dump(db, f, indent=2)

print("Updated POWV database entry with correct 10-protein structure")
