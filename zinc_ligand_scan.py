#!/usr/bin/env python3
"""
Example usage
-----

    python zinc_ligand_scan.py \
           --radius 5.0 \
           --workers 8 \
           --human-only \
           --out zinc_hits.xlsx

Dependencies
-----
Python ≥ 3.8  
pip install requests pandas openpyxl biopython tqdm
"""

from __future__ import annotations
import argparse
import gzip
import os
import shutil
from pathlib import Path
from typing import List, Tuple, Optional

import pandas as pd
import requests
from Bio.PDB import MMCIFParser
from concurrent.futures import ProcessPoolExecutor, as_completed

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def _post_search(payload: dict) -> dict:
    r = requests.post(SEARCH_URL, json=payload, timeout=60)
    r.raise_for_status()
    return r.json()


def fetch_pdb_ids_with_zinc(page_size: int = 1_000, human_only: bool = False) -> List[str]:
    
    zn_node = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
            "operator": "exact_match",
            "value": "ZN",
        },
    }

    if human_only:
        human_node = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                "operator": "exact_match",
                "value": "Homo sapiens",
            },
        }
        query_block = {
            "type": "group",
            "logical_operator": "and",
            "nodes": [zn_node, human_node]
        }
    else:
        query_block = zn_node

    payload = {
        "query": query_block,
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": page_size}},
    }

    ids, start = [], 0
    total_count_known = False
    current_total_count = 0

    while True:
        payload["request_options"]["paginate"]["start"] = start
        
        response_json = _post_search(payload)
        page = response_json.get("result_set") or []
        
        if not total_count_known:
            current_total_count = response_json.get("total_count", 0)
            if current_total_count > 0 :
                 total_count_known = True

        if not page:
            break
            
        ids.extend(hit["identifier"] for hit in page)
        
        if total_count_known and current_total_count > 0 and len(ids) >= current_total_count:
            break
        
        if len(page) < page_size:
            break
            
        start += len(page)

    return ids

PARSER = MMCIFParser(QUIET=True)
FILE_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"


def ensure_cif(pdb_id: str, cache_dir: Path) -> Path:
    pdb_id = pdb_id.lower()
    gz_path = cache_dir / f"{pdb_id}.cif.gz"
    cif_path = cache_dir / f"{pdb_id}.cif"
    if not cif_path.exists():
        if not gz_path.exists():
            url = FILE_URL.format(pdb=pdb_id)
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            gz_path.write_bytes(r.content)
        with gzip.open(gz_path, "rb") as fin, cif_path.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    return cif_path


def scan_structure_for_ligands_near_zn(
    cif_path: Path, radius: float = 5.0
) -> List[Tuple[str, float]]:
    structure = PARSER.get_structure(cif_path.stem, cif_path)
    hits: List[Tuple[str, float]] = []

    for model in structure:
        zn_atoms = [
            atm
            for res in model.get_residues()
            if res.id[0].startswith("H_") and res.get_resname() == "ZN"
            for atm in res.get_atoms()
        ]
        if not zn_atoms:
            continue

        for res in model.get_residues():
            het_flag, _, _ = res.id
            if not het_flag.startswith("H_"):
                continue
            if res.get_resname() in {"ZN", "HOH"}:
                continue

            for zn in zn_atoms:
                for lig_atom in res.get_atoms():
                    if (dist := zn - lig_atom) <= radius:
                        hits.append((res.get_resname(), round(dist, 2)))
                        break
                else:
                    continue
                break

    return hits

def worker(
    pdb_id: str, cache_dir: str, radius: float
) -> Optional[Tuple[str, str, float]]:
    try:
        cif = ensure_cif(pdb_id, Path(cache_dir))
        hits = scan_structure_for_ligands_near_zn(cif, radius)
        if not hits:
            return None
        lig_list = ", ".join(f"{lig}:{d}" for lig, d in hits)
        min_d = min(d for _, d in hits)
        return pdb_id.upper(), lig_list, min_d
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser(description="Find ligands near Zn ions in the PDB.")
    ap.add_argument("--radius", type=float, default=5.0, help="cut-off distance in Å.")
    ap.add_argument(
        "--human-only",
        action="store_true",
        help="restrict search to Homo sapiens entries",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 4,
        help="parallel processes (default = #CPUs)",
    )
    ap.add_argument(
        "--cache",
        default="pdbs",
        help="directory to store downloaded CIFs (default: ./pdbs)",
    )
    ap.add_argument("--out", default="zinc_ligand_hits.xlsx", help="output .xlsx file")
    args = ap.parse_args()

    cache_dir = Path(args.cache)
    cache_dir.mkdir(exist_ok=True)

    print("Fetching PDB IDs …")
    pdb_ids = fetch_pdb_ids_with_zinc(human_only=args.human_only)
    print(f"→ {len(pdb_ids):,} entries contain Zn²⁺")

    print(f"Scanning structures with {args.workers} processes …")
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(worker, pid, cache_dir, args.radius): pid for pid in pdb_ids
        }
        try:
            from tqdm import tqdm

            iterator = tqdm(as_completed(futures), total=len(futures))
        except ImportError:
            iterator = as_completed(futures)

        for fut in iterator:
            result = fut.result()
            if result:
                rows.append(
                    {"PDB_ID": result[0], "LigandNames": result[1], "MinDistance": result[2]}
                )

    print(f"✅ {len(rows):,} structures have at least one ligand within {args.radius} Å")

    df = pd.DataFrame(rows).sort_values("MinDistance")
    df.to_excel(args.out, index=False)
    print(f"Results written to {args.out}")


if __name__ == "__main__":
    main()
