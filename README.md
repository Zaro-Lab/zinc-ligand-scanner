# Zinc Ligand Scanner
Script to find PDB structures with Zinc, filter for human ones, and report ligands near the Zinc ions in an Excel file

![mermaid-diagram-2025-05-11-225715](https://github.com/user-attachments/assets/f8bb8452-3afc-4165-87a5-47a9bc053106)

## Components
- Searches RCSB PDB database for structures containing zinc
- Optional filtering for human-only proteins via `--human-only` flag
- Scans for ligands within specified radius of zinc atoms
- Excludes water molecules and other zinc ions
- Parallel processing with configurable worker count
- Outputs results in Excel format

## Usage

### Basic Usage
```
python zinc_ligand_scan.py
```

### Human-Only Structures
```
python zinc_ligand_scan.py --human-only
```

### Advanced Options
```
python zinc_ligand_scan.py \
    --radius 5.0 \
    --workers 8 \
    --human-only \
    --cache pdbs \
    --out zinc_hits.xlsx
```

## Command Arguments
- `--radius` cutoff distance in Ångstroms default 5.0
- `--human-only` restricts search to Homo sapiens entries
- `--workers` number of parallel processes default = CPU count
- `--cache` directory for downloaded CIFs default ./pdbs
- `--out` output Excel filename default zinc_ligand_hits.xlsx

## Technical Implementation

### PDB Query
- Uses RCSB Search API v2 to find structures with zinc
- For human-only adds taxonomy filter to query
- Paginates results to handle large datasets

### Structure Processing
1. Fetches list of PDB IDs containing zinc
2. Downloads CIF files with caching to avoid redundancy
3. Parses structures using Biopython MMCIFParser
4. Identifies zinc atoms in each structure
5. Scans for ligands within specified radius
6. Records ligand names and minimum distances
7. Exports sorted results to Excel

### Parallel Processing
Uses `ProcessPoolExecutor` to analyze multiple structures simultaneously significantly improving performance on multicore systems

## Output Data
Excel file with columns:
- PDB_ID identifier for structure
- LigandNames list of ligands with distances
- MinDistance closest approach in Ångstroms

## Dependencies
- Python ≥ 3.8
- requests
- pandas
- openpyxl
- biopython
- tqdm
