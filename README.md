# lilac-affinity-mapper

This repository maps experimental affinities from ChEMBL and BindingDB to
[LILAC-DB](https://zenodo.org/records/14835079), introduced in
[Allison Pearl Barkdull et al.'s paper](https://doi.org/10.1038/s42004-025-01472-8).

The workflow is as follows:

- `mappers.lilac_affinity` matches LILAC ligands by exact stereochemical SMILES
  and UniProt accession, producing `mappings/lilac_affinity.csv`.
- `mappers.lilac_sar` retrieves measured analogues for the same target and
  publication, retaining the seed's nonempty Murcko scaffold and series with at
  least three distinct structures including the seed. It writes
  `mappings/lilac_sar.csv` with similarity scores and same-assay flags, excluding
  the seed.
- `mappers.lilac_receptors` adds PDB protein names and species, producing
  `mappings/lilac_sar_receptors.csv`.

We also provide a [website](visualization/index.html) for interactive viewing of
molecules, affinity measurements and bilayer exposure.
