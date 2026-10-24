# LILAC SAR explorer

Open `index.html` in a browser. It is a self-contained snapshot: chemical drawings,
experimental measurements and seed bilayer annotations are embedded. It works
offline without a server or JavaScript dependencies. Links to the source databases
and publications open externally when clicked.

- **Series explorer:** select a study in the sidebar, inspect seed exposure, and
  filter analogue structures and potency plots by endpoint, assay and source.
  Click a structure or plotted measurement to inspect its original records.
- **Campaign overview:** compare all 22 seed–UniProt–publication series, with
  same-assay coverage, seed exposure, exact experimental ranges, experimental
  membrane/system evidence and analogue structure availability. Click any column
  heading to sort, and click again to reverse the order. Choose a range endpoint
  and a structure method (X-ray, EM or any experimental method).
- **CSV export:** the series export includes all source records for that study,
  including the seed, regardless of display filters. The summary export follows
  the overview search, sorting, endpoint and structure-method settings, including
  the experimental context, analogue coverage, PDB IDs and evidence date.

The snapshot contains 16 seed complexes/accessions, 169 distinct analogue SMILES
and 514 analogue measurement rows. Series overlap; do not sum analogue counts as
independent compounds. PDB metadata identifies the seed protein/construct, not
necessarily the assay species or construct. The details panel explains how
bounds, duplicate records, assays and exposure annotations are handled.

There are no FEP predictions or paired membrane/no-membrane results in these
source tables. This page displays the experimental evidence for selecting and
curating the campaign; it cannot estimate a membrane effect.

## Refresh the snapshot

The evidence snapshot includes 32 ChEMBL assays and 346 experimental PDB entries
across the 12 target accessions, retrieved on 24 October 2026. Cell systems and
membrane fractions come from assay descriptions and structured fields. Missing
preparations remain unspecified: HEK293 expression does not establish whether
the measured sample contained intact cells or isolated membranes, and whole
blood is not plasma. Expand **Assay evidence** to review the original text.

Structure counts use exact RDKit canonical isomeric SMILES equality with RCSB
stereochemical ligand descriptors, within entries containing the same UniProt
accession. All chemical matching happens locally; only public identifiers are
sent to the APIs. Salts, protonation states, tautomers and unspecified stereo are
not merged. Zero means no exact match in the retrieved snapshot. Entries can
contain several proteins; the ligand's binding subunit, construct, mutations and
pose still require review. Counts are distinct PDB entries, not ligand copies.
Analogues matching the seed PDB ligand are flagged for seed-stereochemistry review.

To refresh public metadata before rebuilding, run this from the visualization
directory (network access required):

```bash
python3 fetch_evidence.py \
  --mapper /biggin/b211/reub0138/Projects/fep-lipid-exposed/lilac-affinity-mapper \
  --cache /tmp/lilac-evidence-cache
```

`source_evidence.json` is the reproducible input snapshot. The fetcher reuses its
HTTP cache; choose a new cache directory when you want a fresh upstream snapshot.

The HTML generator requires RDKit and uses the saved evidence without network
access. From this directory, run:

```bash
/biggin/b211/reub0138/miniforge3/envs/rdkit-env/bin/python build_visualization.py \
  --mapper /biggin/b211/reub0138/Projects/fep-lipid-exposed/lilac-affinity-mapper
```

`--lilac` can override the default local LILAC non-redundant database directory.
The generator reads `mappings/lilac_sar_receptors.csv`, joins seed measurements
from `mappings/lilac_affinity.csv`, and uses each seed's LILAC CSV/SDF files.
Normalized DOIs preserve joins across differences in capitalization. Source file
hashes and the build time are embedded in the page. Exposure highlights are
transferred only after checking atom elements and coordinates against the SDF.

Edit `explorer.template.html` to change the interface, then rerun the generator.
