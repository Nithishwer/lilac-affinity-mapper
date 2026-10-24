import json
from functools import cache
from urllib.parse import urlencode
from urllib.request import urlopen

BASE_URL = "https://data.rcsb.org/graphql"
QUERY = """
query ($pdb_id: String!) {
  entry(entry_id: $pdb_id) {
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { uniprot_ids }
      rcsb_entity_source_organism { ncbi_scientific_name }
    }
  }
}
"""


@cache
def get_proteins(pdb_code: str) -> list[dict]:
    """Return PDB polymer entities with protein names, accessions and source species."""
    query = urlencode(
        {"query": QUERY, "variables": json.dumps({"pdb_id": pdb_code.upper()})}
    )
    with urlopen(f"{BASE_URL}?{query}", timeout=60) as response:
        body = json.load(response)
    if body.get("errors"):
        raise ValueError(f"PDB query failed for {pdb_code}: {body['errors']}")
    entry = body["data"]["entry"]
    if entry is None:
        raise ValueError(f"PDB entry not found: {pdb_code}")
    return entry.get("polymer_entities") or []
