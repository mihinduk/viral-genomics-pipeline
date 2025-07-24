# Development Tools

This directory contains scripts and utilities for pipeline developers and curators. **End users do not need these tools.**

## Directory Structure

- database_management/ - Database update and maintenance scripts
- curation/ - Genome curation and positioning tools  
- README.md - This file

## Database Management

**Purpose**: Scripts for updating virus configuration databases

- update_powv_entry.py - Updates POWV entry with correct 10-protein structure
- Future: Scripts for adding new virus families, updating gene coordinates

## Curation Tools

**Purpose**: Tools for fine-tuning visualization positioning and quality

- add_powv_positioning.py - Adds POWV-specific gene label positioning
- Future: Label positioning optimization tools, color scheme generators

## Usage Notes

- These are **one-time** or **maintenance** scripts
- Run only when updating the pipeline or adding new viruses
- Not part of the standard analysis workflow
- Require developer knowledge of the pipeline internals

## For End Users

If you're running viral genome analysis, you only need:
- Scripts in viral_pipeline/analysis/ (main analysis)
- Scripts in viral_pipeline/visualization/ (plotting)
- README_MODULES.md (usage instructions)

The tools in this dev_tools/ directory are for pipeline developers only.
