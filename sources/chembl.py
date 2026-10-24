import json
from functools import cache
from urllib.parse import urlencode, urljoin
from urllib.request import urlopen

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/"


@cache
def get_document_doi(document_id):
    if not document_id:
        return ""
    with urlopen(f"{BASE_URL}document/{document_id}.json", timeout=60) as response:
        return json.load(response).get("doi") or ""


def _fetch(resource, key, **filters):
    url = f"{BASE_URL}{resource}.json?{urlencode({'limit': 1000, **filters})}"
    records = []
    while url:
        with urlopen(url, timeout=60) as response:
            page = json.load(response)
        records.extend(page[key])
        next_page = page["page_meta"]["next"]
        url = urljoin(BASE_URL, next_page) if next_page else None
    return records


@cache
def _get_targets(uniprot_id):
    return _fetch(
        "target",
        "targets",
        target_components__accession=uniprot_id,
        target_type="SINGLE PROTEIN",
    )


def get_binding_data(uniprot_id: str) -> list[dict]:
    """Return raw binding and functional activities for matching single-protein targets."""
    activities = []
    for target in _get_targets(uniprot_id):
        activities.extend(
            _fetch(
                "activity",
                "activities",
                target_chembl_id=target["target_chembl_id"],
                assay_type__in="B,F",
            )
        )
    return activities


def get_study_data(uniprot_id: str, doi: str) -> list[dict]:
    """Return binding and functional activities for one target and publication."""
    targets = _get_targets(uniprot_id)
    if not targets:
        return []
    documents = _fetch("document", "documents", doi__iexact=doi)
    if not documents:
        return []
    return _fetch(
        "activity",
        "activities",
        target_chembl_id__in=",".join(target["target_chembl_id"] for target in targets),
        document_chembl_id__in=",".join(
            document["document_chembl_id"] for document in documents
        ),
        assay_type__in="B,F",
    )
