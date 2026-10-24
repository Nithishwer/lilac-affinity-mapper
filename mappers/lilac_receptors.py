import csv
from pathlib import Path

from sources import pdb

INPUT = Path(__file__).resolve().parents[1] / "mappings" / "lilac_sar.csv"
OUTPUT = INPUT.with_name("lilac_sar_receptors.csv")
COLUMNS = ["pdb_protein_name", "pdb_species"]


def receptor_fields(pdb_code, uniprot_id):
    names = set()
    species = set()
    for protein in pdb.get_proteins(pdb_code):
        identifiers = protein.get("rcsb_polymer_entity_container_identifiers") or {}
        # Match the target accession instead of taking an unrelated protein in the entry.
        if uniprot_id not in (identifiers.get("uniprot_ids") or []):
            continue
        name = (protein.get("rcsb_polymer_entity") or {}).get("pdbx_description")
        if name:
            names.add(name)
        species.update(
            organism["ncbi_scientific_name"]
            for organism in protein.get("rcsb_entity_source_organism") or []
            if organism.get("ncbi_scientific_name")
        )
    return {
        "pdb_protein_name": "; ".join(sorted(names)),
        "pdb_species": "; ".join(sorted(species)),
    }


def main():
    with INPUT.open(newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames + COLUMNS
        rows = list(reader)

    for row in rows:
        row.update(receptor_fields(row["pdb_code"], row["uniprot_id"]))

    # Keep the previous output intact if a PDB request fails.
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
