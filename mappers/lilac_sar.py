import csv
import re
from collections import defaultdict
from functools import cache
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from mappers.lilac_affinity import COLUMNS as AFFINITY_COLUMNS
from mappers.lilac_affinity import affinity_record, canonical_smiles
from sources import bindingdb, chembl

INPUT = Path(__file__).resolve().parents[1] / "mappings" / "lilac_affinity.csv"
OUTPUT = INPUT.with_name("lilac_sar.csv")
MIN_SERIES_SIZE = 3
COLUMNS = AFFINITY_COLUMNS + ["similarity", "assay_id", "same_assay", "series_size"]
FINGERPRINT = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048,
    includeChirality=True,
)


def publication_doi(value):
    doi = re.sub(
        r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", (value or "").strip().lower()
    )
    return "" if doi.startswith("10.7270/") else doi


@cache
def structure(smiles):
    smiles = canonical_smiles(smiles)
    if not smiles:
        return None
    molecule = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return smiles, scaffold, FINGERPRINT.GetFingerprint(molecule)


def series_rows(seed, records):
    seed_structure = structure(seed["smiles"])
    if not seed_structure or not seed_structure[1]:
        return []
    seed_smiles, seed_scaffold, seed_fingerprint = seed_structure
    seed_assays = defaultdict(set)
    candidates = []
    for record in records:
        candidate = structure(record["smiles"])
        if candidate[0] == seed_smiles:
            if record["source"] == "chembl" and record["assay_id"]:
                seed_assays[record["affinity_measurement_type"]].add(record["assay_id"])
        elif candidate[1] == seed_scaffold:
            candidates.append(record)

    # Count structures, not measurements or IDs from different databases.
    series_size = len({record["smiles"] for record in candidates} | {seed_smiles})
    if series_size < MIN_SERIES_SIZE:
        return []
    rows = []
    for candidate in candidates:
        assays = seed_assays[candidate["affinity_measurement_type"]]
        same_assay = "unknown"
        if candidate["source"] == "chembl" and candidate["assay_id"] and assays:
            same_assay = str(candidate["assay_id"] in assays).lower()
        similarity = DataStructs.TanimotoSimilarity(
            seed_fingerprint,
            structure(candidate["smiles"])[2],
        )
        rows.append(
            {
                **{
                    key: seed[key]
                    for key in (
                        "pdb_code",
                        "pdb_ligand_id",
                        "pdb_ligand_instance",
                        "uniprot_id",
                    )
                },
                **candidate,
                "similarity": f"{similarity:.4f}",
                "same_assay": same_assay,
                "series_size": series_size,
            }
        )
    return sorted(rows, key=lambda row: float(row["similarity"]), reverse=True)


def main():
    studies = defaultdict(lambda: defaultdict(dict))
    with INPUT.open(newline="") as stream:
        for row in csv.DictReader(stream):
            doi = publication_doi(row["doi"])
            if not row["affinity_value"] or not doi:
                continue
            seed_structure = structure(row["smiles"])
            if not seed_structure or not seed_structure[1]:
                continue
            key = tuple(
                row[field]
                for field in ("pdb_code", "pdb_ligand_id", "pdb_ligand_instance")
            )
            studies[row["uniprot_id"]][doi][key] = {
                **row,
                "smiles": seed_structure[0],
                "doi": doi,
            }

    output_rows = []
    for uniprot, publications in studies.items():
        print(f"Querying bindingdb: {uniprot}", flush=True)
        bindingdb_studies = defaultdict(list)
        for record in bindingdb.get_binding_data(uniprot):
            doi = publication_doi(record.get("doi"))
            if doi in publications:
                bindingdb_studies[doi].append(record)
        for doi, seeds in publications.items():
            print(f"Querying chembl: {uniprot} / {doi}", flush=True)
            records = []
            for source, raw_records in (
                ("chembl", chembl.get_study_data(uniprot, doi)),
                ("bindingdb", bindingdb_studies[doi]),
            ):
                for raw in raw_records:
                    smiles = (
                        raw.get("canonical_smiles")
                        if source == "chembl"
                        else raw.get("smile")
                    )
                    candidate = structure(smiles)
                    if not candidate or not candidate[1]:
                        continue
                    affinity = affinity_record(raw, source)
                    if not affinity or publication_doi(affinity["doi"]) != doi:
                        continue
                    records.append(
                        {
                            **affinity,
                            "doi": doi,
                            "source": source,
                            "smiles": candidate[0],
                            "assay_id": raw.get("assay_chembl_id") or ""
                            if source == "chembl"
                            else "",
                        }
                    )
            for seed in seeds.values():
                output_rows.extend(series_rows(seed, records))

    # Keep the previous output intact if any source request fails.
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Saved {len(output_rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
