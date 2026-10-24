# lilac-affinity-mapper

Map LILAC complexes to ChEMBL and BindingDB measurements using exact canonical
SMILES, including stereochemistry. Salts, protonation states and tautomers are not
merged. Each UniProt accession in a multi-protein entry is queried separately;
this does not establish which subunit or construct was measured.

```bash
python -m pip install -r requirements.txt
python -m mappers.lilac_affinity
```

Set `LILAC_DB` in `mappers/lilac_affinity.py` to the database folder. The script
reads its Non-Redundant CSV and writes `mappings/lilac_affinity.csv`.

Each measurement gets a row. A complex/source with no matching measurement has
empty measurement fields. Request failures stop the run before writing output.

ChEMBL includes binding (B) and functional (F) assays.
Retains IC50, Ki, Kd and their p variants; other endpoint types are omitted.
Concentration values are in nM; pIC50, pKi and pKd values remain logarithmic.
Inequality qualifiers stay in `affinity_value`. Experimental methods are assigned
from explicit assay-description keywords; missing or ambiguous methods stay blank.
BindingDB's current endpoint supplies no assay methods, and its dataset DOIs are
excluded from the publication DOI column. ChEMBL publication DOIs are fetched
from its document endpoint. Measurements are not averaged or deduplicated across
sources.

Find measured analogues from the same study:

```bash
python -m mappers.lilac_sar
```

Reads `mappings/lilac_affinity.csv` and writes `mappings/lilac_sar.csv`. Candidates
must share the seed's nonempty Murcko scaffold. `MIN_SERIES_SIZE` defaults to three
distinct measured structures including the seed; cross-source duplicates count
once. The seed itself is excluded from candidate rows. Ring-changing analogues
and acyclic scaffolds are not covered by this conservative screen.

The original PDB fields identify the seed; `source_id` and `smiles` describe the
candidate. Similarity is chirality-aware Morgan Tanimoto (radius 2, 2048 bits).
`same_assay` compares ChEMBL assay IDs for the same measurement type, and is
`unknown` when the seed's assay cannot be identified or for BindingDB records.
Different assays and endpoints remain separate measurement rows. `series_size`
counts measured structures in the study, not only those sharing an assay.
Only rows with measurements are written; missing DOI seeds are skipped. DOI
prefixes and case are normalized. Same-study series require review before being
treated as measurements from one experiment.
