import csv
import re
from collections import defaultdict
from decimal import Decimal
from functools import cache
from pathlib import Path

from rdkit import Chem

from sources import bindingdb, chembl

LILAC_DB = Path("/biggin/b211/reub0138/Projects/fep-lipid-exposed/LILAC_DB")
OUTPUT = Path(__file__).resolve().parents[1] / "mappings" / "lilac_affinity.csv"
COLUMNS = [
    "pdb_code",
    "pdb_ligand_id",
    "pdb_ligand_instance",
    "uniprot_id",
    "source",
    "source_id",
    "smiles",
    "exp_type",
    "affinity_value",
    "affinity_measurement_type",
    "doi",
]
MEASUREMENT_TYPES = {
    name.lower(): name for name in ("IC50", "Ki", "Kd", "pIC50", "pKi", "pKd")
}
NM_FACTORS = {
    "M": "1e9",
    "mM": "1e6",
    "uM": "1e3",
    "nM": "1",
    "pM": "1e-3",
    "fM": "1e-6",
}
METHODS = {
    "SPR": r"\bspr\b|surface plasmon resonance",
    "ITC": r"\bitc\b|isothermal titration calorimetry",
    "MST": r"\bmst\b|microscale thermophoresis",
    "BLI": r"\bbli\b|bio.?layer interferometry",
    "radioligand_binding": r"radioligand|radiolabel|tritiated|\[3h\]|\[125i\]",
    "fluorescence_polarization": r"\bfp\b|fluorescence polarization",
}


@cache
def canonical_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles or "")
    if molecule is None or molecule.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def experimental_method(description):
    methods = [
        name
        for name, pattern in METHODS.items()
        if re.search(pattern, description or "", re.IGNORECASE)
    ]
    return methods[0] if len(methods) == 1 else ""


def measurement(kind, value, units, relation=""):
    kind = MEASUREMENT_TYPES.get((kind or "").lower())
    match = re.fullmatch(
        r"\s*([<>=~≤≥]*)\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*",
        str(value),
    )
    if not kind or not match:
        return None
    number = Decimal(match[2])
    if not kind.startswith("p"):
        factor = NM_FACTORS.get((units or "").replace("µ", "u").replace("μ", "u"))
        if factor is None:
            return None
        number *= Decimal(factor)
    qualifier = match[1] or relation or ""
    return {
        "affinity_value": ("" if qualifier == "=" else qualifier)
        + format(number.normalize(), "f"),
        "affinity_measurement_type": kind,
    }


def affinity_record(record, source):
    if source == "chembl":
        # Preserve reported logarithmic endpoints rather than converting to concentration.
        if (record.get("type") or "").lower() in {"pic50", "pki", "pkd"}:
            result = measurement(
                record["type"], record.get("value"), None, record.get("relation")
            )
        else:
            result = measurement(
                record.get("standard_type"),
                record.get("standard_value"),
                record.get("standard_units"),
                record.get("standard_relation"),
            )
        if result:
            result.update(
                source_id=record["molecule_chembl_id"],
                exp_type=experimental_method(record.get("assay_description")),
                doi=chembl.get_document_doi(record.get("document_chembl_id")),
            )
    else:
        # BindingDB's REST endpoint returns concentration affinities in nM.
        result = measurement(record.get("affinity_type"), record.get("affinity"), "nM")
        if result:
            doi = record.get("doi") or ""
            result.update(
                source_id=record["monomerid"],
                exp_type="",
                # A BindingDB dataset DOI is not the publication DOI.
                doi="" if "10.7270/" in doi.lower() else doi,
            )
    return result


def main():
    input_csv = LILAC_DB / "LILAC_DB_Non_Redundant" / "LILAC_DB_Non_Redundant.csv"
    complexes_by_uniprot = defaultdict(list)
    with input_csv.open(newline="") as stream:
        for complex_row in csv.DictReader(stream):
            smiles = canonical_smiles(complex_row["Ligand SMILES"])
            if smiles is None:
                raise ValueError(f"Invalid LILAC SMILES: {complex_row['pdb_id']}")
            for uniprot in dict.fromkeys(complex_row["UniProt Accessions"].split(",")):
                complexes_by_uniprot[uniprot.strip()].append(
                    {
                        "pdb_code": complex_row["pdb_id"],
                        "pdb_ligand_id": complex_row["ligand_id"],
                        "pdb_ligand_instance": complex_row["ligand_copy"],
                        "uniprot_id": uniprot.strip(),
                        "smiles": smiles,
                    }
                )

    output_rows = []
    for uniprot, complexes in complexes_by_uniprot.items():
        wanted_smiles = {row["smiles"] for row in complexes}
        for source_name, source in (("chembl", chembl), ("bindingdb", bindingdb)):
            print(f"Querying {source_name}: {uniprot}", flush=True)
            matches = defaultdict(list)
            for record in source.get_binding_data(uniprot):
                smiles = canonical_smiles(
                    record.get("canonical_smiles")
                    if source_name == "chembl"
                    else record.get("smile")
                )
                if smiles in wanted_smiles:
                    affinity = affinity_record(record, source_name)
                    if affinity:
                        matches[smiles].append(affinity)
            for complex_row in complexes:
                for affinity in matches[complex_row["smiles"]] or [{}]:
                    output_rows.append(
                        {**complex_row, "source": source_name, **affinity}
                    )

    # Write only after all requests succeed, so failures do not create a partial CSV.
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Saved {len(output_rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
