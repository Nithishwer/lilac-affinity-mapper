"""Download public assay/PDB metadata by identifier; never upload compound structures."""

import argparse
import csv
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ENTRY_FIELDS = """{
  rcsb_id exptl { method }
  nonpolymer_entities {
    pdbx_entity_nonpoly { comp_id }
    nonpolymer_comp { rcsb_chem_comp_descriptor { SMILES_stereo } }
  }
  polymer_entities {
    rcsb_polymer_entity_container_identifiers { uniprot_ids }
  }
}"""


def fetch_json(url, cache):
    path = cache / (hashlib.sha256(url.encode()).hexdigest() + ".json")
    if path.exists():
        return json.loads(path.read_text())
    for attempt in range(4):
        try:
            request = Request(url, headers={"User-Agent": "LILAC-SAR-evidence/1.0"})
            with urlopen(request, timeout=90) as response:
                result = {} if response.status == 204 else json.load(response)
            path.write_text(json.dumps(result))
            return result
        except (HTTPError, URLError, TimeoutError) as error:
            if isinstance(error, HTTPError) and error.code not in (
                429,
                500,
                502,
                503,
                504,
            ):
                raise
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))


def target_entries(accession, cache):
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "results_content_type": ["experimental"],
        },
    }
    data = fetch_json(
        "https://search.rcsb.org/rcsbsearch/v2/query?"
        + urlencode({"json": json.dumps(query)}),
        cache,
    )
    entries = [item["identifier"] for item in data.get("result_set", [])]
    assert len(entries) == data.get("total_count", 0), "Incomplete PDB search results"
    return accession, entries


def entry_batch(identifiers, cache):
    query = "{entries(entry_ids:" + json.dumps(identifiers) + ")" + ENTRY_FIELDS + "}"
    data = fetch_json(
        "https://data.rcsb.org/graphql?" + urlencode({"query": query}), cache
    )
    if data.get("errors"):
        raise ValueError(data["errors"])
    entries = data["data"]["entries"]
    if any(entry is None for entry in entries) or len(entries) != len(identifiers):
        raise ValueError("Missing PDB entries in metadata response")
    return entries


def main(mapper, output, cache):
    cache.mkdir(parents=True, exist_ok=True)
    with (mapper / "mappings" / "lilac_sar_receptors.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assay_ids = sorted({row["assay_id"] for row in rows if row["assay_id"]})
    accessions = sorted({row["uniprot_id"] for row in rows})

    # Every outbound request contains public identifiers only. SMILES stay local.
    def get_assay(identifier):
        url = f"https://www.ebi.ac.uk/chembl/api/data/assay/{identifier}.json"
        return identifier, {**fetch_json(url, cache), "source_url": url}

    with ThreadPoolExecutor(max_workers=4) as pool:
        assays = dict(pool.map(get_assay, assay_ids))
        print(f"Retrieved {len(assays)} public assay descriptions", flush=True)
        targets = dict(
            pool.map(lambda identifier: target_entries(identifier, cache), accessions)
        )
        for identifier, entries in targets.items():
            print(f"{identifier}: {len(entries)} experimental PDB entries", flush=True)
        identifiers = sorted(set().union(*map(set, targets.values())))
        batches = [identifiers[i : i + 40] for i in range(0, len(identifiers), 40)]
        entries = {}
        for batch in pool.map(lambda values: entry_batch(values, cache), batches):
            entries.update((entry["rcsb_id"], entry) for entry in batch)
            print(
                f"Retrieved structure metadata: {len(entries)}/{len(identifiers)}",
                flush=True,
            )

    result = {
        "retrieved": datetime.now(timezone.utc).isoformat(),
        "assays": assays,
        "targets": targets,
        "entries": entries,
        "structure_source": "https://search.rcsb.org/",
        "structure_metadata_source": "https://data.rcsb.org/",
        "matching": "Local RDKit canonical isomeric SMILES equality; no salt, tautomer or protonation normalization",
    }
    output.write_text(json.dumps(result, indent=2))
    print(f"Saved public evidence snapshot: {output}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapper", type=Path, default=Path(__file__).resolve().parent.parent
    )
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("source_evidence.json")
    )
    parser.add_argument(
        "--cache", type=Path, default=Path(__file__).parent / ".evidence-cache"
    )
    args = parser.parse_args()
    main(args.mapper, args.output, args.cache)
