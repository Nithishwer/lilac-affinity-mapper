"""Build the offline LILAC SAR explorer from the mapper and local LILAC database."""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from functools import cache
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def normalize_doi(value):
    return re.sub(
        r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", value.strip().lower()
    )


def seed_key(row):
    return tuple(
        row[field]
        for field in ("pdb_code", "pdb_ligand_id", "pdb_ligand_instance", "uniprot_id")
    )


@cache
def canonical_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles or "")
    return Chem.MolToSmiles(molecule, isomericSmiles=True) if molecule else None


def assay_context(assay):
    """Summarize only systems explicitly named in the assay description/fields."""
    text = assay.get("description") or ""
    lower = text.lower()
    cell_patterns = [
        (r"hek[- ]?293t", "HEK293T"),
        (r"hek[- ]?293", "HEK293"),
        (r"\b293e\b", "293E"),
        (r"cho-rd-hga16", "CHO-RD-HGA16"),
        (r"cho-?k1", "CHO-K1"),
        (r"\bcho\b", "CHO"),
        (r"\bsf9\b", "Sf9"),
        (r"\ba549\b", "A549"),
        (r"\bcaco-?2\b", "Caco-2"),
        (r"\bhepg2\b", "HepG2"),
        (r"\bu937\b", "U937"),
        (r"\ba9 l\b", "A9 L"),
    ]
    cell = next(
        (name for pattern, name in cell_patterns if re.search(pattern, lower)), ""
    )
    cell = cell or assay.get("assay_cell_type") or ""
    fraction = assay.get("assay_subcellular_fraction") or ""
    prep = "Unspecified"
    status = "Reported system; membrane preparation unspecified"
    if re.search(r"\bwhole blood\b", lower):
        label, prep = "Human whole blood", "Whole blood"
        status = "Reported biological matrix"
    elif re.search(r"\bplasma\b", lower) and "plasma membrane" not in lower:
        label, prep = "Plasma", "Plasma"
        status = "Reported biological matrix"
    elif re.search(r"\bmicrosomes?\b", lower):
        label, prep = f"{cell or 'Unspecified source'} microsomes", "Microsomes"
        status = "Reported membrane fraction"
    elif re.search(r"liposome|proteoliposome|nanodisc|lipid suspension", lower):
        prep = re.search(
            r"proteoliposomes?|liposomes?|nanodiscs?|lipid suspension", lower
        )[0]
        label, status = prep.capitalize(), "Reported reconstituted system"
    elif fraction:
        label, prep = f"{cell + ' ' if cell else ''}{fraction}", fraction
        status = "Reported subcellular fraction"
    elif cell:
        label = f"{cell} {'expression system' if re.search(r'express|transfect|baculovirus', lower) else 'cell assay'}"
    elif "cell-free" in lower:
        label = "Cell-free assay"
    elif "cell based" in lower or "cell-based" in lower:
        label = "Cell-based assay; cell type unspecified"
    elif assay.get("assay_tissue") or re.search(
        r"canine (?:liver|intestine|adrenal)", lower
    ):
        tissue = assay.get("assay_tissue") or "liver"
        label = (
            f"{assay.get('assay_organism') or 'Unspecified organism'} {tissue.lower()}"
        )
    else:
        label, status = "Experimental system unspecified", "Unspecified in ChEMBL"
    return {"label": label, "preparation": prep, "status": status}


def structure_index(evidence):
    """Compare downloaded stereochemical descriptors locally for each target."""
    index = defaultdict(list)
    entries = evidence["entries"]
    for accession, identifiers in evidence["targets"].items():
        for identifier in identifiers:
            entry = entries[identifier]
            matches = {}
            for entity in entry.get("nonpolymer_entities") or []:
                descriptor = (entity.get("nonpolymer_comp") or {}).get(
                    "rcsb_chem_comp_descriptor"
                ) or {}
                smiles = canonical_smiles(descriptor.get("SMILES_stereo") or "")
                if smiles:
                    matches.setdefault(smiles, set()).add(
                        entity["pdbx_entity_nonpoly"]["comp_id"]
                    )
            for smiles, components in matches.items():
                index[(accession, smiles)].append(
                    {
                        "pdb": identifier,
                        "components": sorted(components),
                        "methods": [
                            item["method"] for item in entry.get("exptl") or []
                        ],
                    }
                )
    return index


def molecule_svg(molecule, highlights=(), color=(0.79, 0.72, 0.96)):
    molecule = Chem.Mol(molecule)
    rdDepictor.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(420, 240)
    drawer.drawOptions().clearBackground = False
    drawer.drawOptions().padding = 0.10
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=list(highlights),
        highlightAtomColors={index: color for index in highlights},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText().split("?>")[-1]


def seed_geometry(database, metadata):
    base = (
        database
        / metadata["pdb_id"].lower()
        / metadata["ligand_id"]
        / (metadata["ligand_id"] + "_" + metadata["ligand_copy"].replace(" ", ""))
    )
    atoms = read_csv(base.with_suffix(".csv"))
    molecule = Chem.SDMolSupplier(str(base.with_suffix(".sdf")), removeHs=False)[0]
    if molecule is None:
        raise ValueError(f"Cannot read {base}.sdf")

    # Check atom identity and coordinates before transferring exposure highlights.
    coordinates = molecule.GetConformer()
    for atom in atoms:
        index = int(atom["atom_idx"])
        point = coordinates.GetAtomPosition(index)
        assert molecule.GetAtomWithIdx(index).GetSymbol() == atom["element"]
        assert all(
            abs(float(atom[key]) - value) < 0.02
            for key, value in zip(
                ("x_coord", "y_coord", "z_coord"), (point.x, point.y, point.z)
            )
        ), f"SDF/CSV atom mismatch in {base}"

    heavy_atoms = [atom for atom in atoms if atom["element"] != "H"]
    heavy_indices = {int(atom["atom_idx"]) for atom in heavy_atoms}
    exposed_indices = {
        int(atom["atom_idx"])
        for atom in heavy_atoms
        if atom["touches_bilayer"].lower() == "true"
    }
    # RDKit retains this property when explicit hydrogens are removed.
    for atom in molecule.GetAtoms():
        atom.SetIntProp("original_index", atom.GetIdx())
    drawing_molecule = Chem.RemoveHs(molecule)
    highlights = [
        atom.GetIdx()
        for atom in drawing_molecule.GetAtoms()
        if atom.GetIntProp("original_index") in exposed_indices
    ]
    return {
        "svg": molecule_svg(drawing_molecule, highlights),
        "atoms": [
            {
                "index": int(atom["atom_idx"]),
                "element": atom["element"],
                "x": float(atom["x_coord"]),
                "y": float(atom["y_coord"]),
                "z": float(atom["z_coord"]),
                "exposed": int(atom["atom_idx"]) in exposed_indices,
            }
            for atom in heavy_atoms
        ],
        "bonds": [
            [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
            for bond in molecule.GetBonds()
            if bond.GetBeginAtomIdx() in heavy_indices
            and bond.GetEndAtomIdx() in heavy_indices
        ],
        "exposed": len(exposed_indices),
        "heavy": len(heavy_atoms),
        "halfThickness": float(metadata["1/2 Lipid Tail Thickness"]),
        "membrane": metadata["Membrane"],
        "classification": metadata["Protein Classifications"],
        "annotationPath": str(base.with_suffix(".csv")),
    }


def build(mapper, lilac, template, output, evidence_path):
    mappings = mapper / "mappings"
    sar = read_csv(mappings / "lilac_sar_receptors.csv")
    affinity = read_csv(mappings / "lilac_affinity.csv")
    evidence = json.loads(evidence_path.read_text())
    pdb_index = structure_index(evidence)
    for assay in evidence["assays"].values():
        assay["context"] = assay_context(assay)
    metadata = {
        (row["pdb_id"], row["ligand_id"], row["ligand_copy"]): row
        for row in read_csv(lilac / "LILAC_DB_Non_Redundant.csv")
    }
    groups = defaultdict(list)
    for row in sar:
        groups[(*seed_key(row), normalize_doi(row["doi"]))].append(row)
    seed_measurements = defaultdict(list)
    for row in affinity:
        if row["affinity_value"]:
            seed_measurements[(*seed_key(row), normalize_doi(row["doi"]))].append(row)

    seeds, series, structures = {}, [], {}
    for key, rows in groups.items():
        pdb, ligand, instance, uniprot, doi = key
        seed_id = "|".join(key[:4])
        if seed_id not in seeds:
            meta = metadata[(pdb, ligand, instance)]
            seeds[seed_id] = {
                "pdb": pdb,
                "ligand": ligand,
                "instance": instance,
                "uniprot": uniprot,
                "name": rows[0]["pdb_protein_name"] or meta["Protein Name"],
                "pdbName": rows[0]["pdb_protein_name"],
                "species": rows[0]["pdb_species"],
                "smiles": Chem.MolToSmiles(Chem.MolFromSmiles(meta["Ligand SMILES"])),
                "methods": [
                    item["method"]
                    for item in evidence["entries"][pdb].get("exptl") or []
                ],
                "pdbMatches": pdb_index[
                    (uniprot, canonical_smiles(meta["Ligand SMILES"]))
                ],
                **seed_geometry(lilac, meta),
            }
        seed = seeds[seed_id]
        scaffold = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(seed["smiles"]))
        compounds = {}
        for row in rows:
            smiles = row["smiles"]
            compound_id = hashlib.sha256(smiles.encode()).hexdigest()[:12]
            if smiles not in structures:
                structures[smiles] = compound_id
            if compound_id not in compounds:
                molecule = Chem.MolFromSmiles(smiles)
                match = set(molecule.GetSubstructMatch(scaffold))
                highlights = [
                    i for i in range(molecule.GetNumAtoms()) if i not in match
                ]
                compounds[compound_id] = {
                    "id": compound_id,
                    "smiles": smiles,
                    "svg": molecule_svg(molecule, highlights, (0.75, 0.86, 0.96)),
                    "similarity": float(row["similarity"]),
                    "records": [],
                    "pdbMatches": pdb_index[(uniprot, canonical_smiles(smiles))],
                }
            compounds[compound_id]["records"].append(row)
        series.append(
            {
                "id": f"series-{len(series) + 1:02}",
                "seedId": seed_id,
                "doi": doi,
                "reportedSize": int(rows[0]["series_size"]),
                "compounds": list(compounds.values()),
                "seedRecords": seed_measurements[key],
            }
        )

    payload = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "sources": [
            {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in (
                mappings / "lilac_sar_receptors.csv",
                mappings / "lilac_affinity.csv",
                lilac / "LILAC_DB_Non_Redundant.csv",
                evidence_path,
            )
        ],
        "rowCount": len(sar),
        "uniqueAnalogues": len(structures),
        "seeds": seeds,
        "series": series,
        "assays": evidence["assays"],
        "evidenceRetrieved": evidence["retrieved"],
        "structureSearch": {
            "targetEntryCounts": {
                key: len(values) for key, values in evidence["targets"].items()
            },
            "entriesRetrieved": len(evidence["entries"]),
            "matching": evidence["matching"],
            "source": evidence["structure_source"],
        },
    }
    # Escape the HTML script terminator even if future source strings contain it.
    encoded = json.dumps(payload, separators=(",", ":")).replace("<", "\\u003c")
    output.write_text(template.read_text().replace("__LILAC_DATA__", encoded))
    print(
        f"Wrote {output}: {len(series)} series, {len(seeds)} seeds, "
        f"{len(sar)} analogue records, {len(structures)} unique analogues"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapper",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument(
        "--lilac",
        type=Path,
        default=Path(
            "/biggin/b211/reub0138/Projects/fep-lipid-exposed/LILAC_DB/LILAC_DB_Non_Redundant"
        ),
    )
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    build(
        args.mapper,
        args.lilac,
        folder / "explorer.template.html",
        folder / "index.html",
        folder / "source_evidence.json",
    )
