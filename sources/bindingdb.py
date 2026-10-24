import json
from urllib.parse import urlencode
from urllib.request import urlopen

BASE_URL = "https://bindingdb.org/rest/getLigandsByUniprots"


def get_binding_data(uniprot_id: str) -> list[dict]:
    """Return raw BindingDB affinity records for a UniProt accession."""
    query = urlencode({"uniprot": uniprot_id, "response": "application/json"})
    with urlopen(f"{BASE_URL}?{query}", timeout=60) as response:
        body = response.read()
    if not body.strip():
        return []
    return json.loads(body)["getLindsByUniprotsResponse"]["affinities"]
