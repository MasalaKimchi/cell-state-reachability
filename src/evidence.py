"""Orthogonal external evidence + pathway interpretation for a nominated gene.

Evidence is SUPPORT, never proof. We surface prior knowledge so a human can judge a
hypothesis — we never let a literature hit inflate a claim of causality or rescue.

ALL SOURCES ARE OPEN AND REQUIRE NO AUTHENTICATION:
  - PubMed E-utilities (esearch)            — co-mention counts / example PMIDs
  - Open Targets GraphQL (public endpoint)  — target -> associated diseases
  - Enrichr via gseapy                      — pathway / GO over-representation
No CZI / Synapse / Wiley login is needed anywhere in this module. (This file also
absorbs the former pathways.py so external interpretation lives in one place.)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import pandas as pd
import requests

PUBMED = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"

DEFAULT_LIBRARIES = [
    "GO_Biological_Process_2021",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]


# --------------------------------------------------------------------------- #
# Per-gene evidence
# --------------------------------------------------------------------------- #
@dataclass
class GeneEvidence:
    gene: str
    pubmed_count: int
    example_pmids: list[str]
    opentargets_diseases: list[dict] = field(default_factory=list)

    def external_score(self) -> float:
        """Squash co-mention counts into [0,1]. Deliberately conservative:
        absence of literature is NOT evidence of absence, so 0 mentions -> 0.3, not 0.
        """
        base = 0.3 + 0.7 * (1 - math.exp(-self.pubmed_count / 10.0))
        return float(min(base, 1.0))


def pubmed_comention(gene: str, process: str, retmax: int = 5,
                     timeout: int = 20) -> tuple[int, list[str]]:
    """Count PubMed records co-mentioning a gene and the target process (no auth)."""
    term = f'("{gene}"[Title/Abstract]) AND ("{process}"[Title/Abstract])'
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax}
    r = requests.get(PUBMED, params=params, timeout=timeout)
    r.raise_for_status()
    res = r.json().get("esearchresult", {})
    return int(res.get("count", 0)), list(res.get("idlist", []))


def opentargets_diseases(gene_symbol: str, size: int = 10,
                         timeout: int = 20) -> list[dict]:
    """Top associated diseases for a target from Open Targets public GraphQL (no auth).

    Two steps: resolve the gene symbol to an Ensembl target id, then pull the ranked
    associatedDiseases. Returns [{disease, score}] sorted by association score. Never
    raises — the evidence layer must not crash the pipeline.
    """
    resolve = """
    query resolve($q: String!) {
      search(queryString: $q, entityNames: ["target"]) { hits { id name } }
    }"""
    assoc = """
    query assoc($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        associatedDiseases(page: {index: 0, size: $size}) {
          rows { disease { name } score }
        }
      }
    }"""
    try:
        r = requests.post(OPENTARGETS, json={"query": resolve,
                          "variables": {"q": gene_symbol}}, timeout=timeout)
        r.raise_for_status()
        hits = r.json().get("data", {}).get("search", {}).get("hits", [])
        if not hits:
            return []
        target_id = hits[0]["id"]
        r2 = requests.post(OPENTARGETS, json={"query": assoc,
                           "variables": {"id": target_id, "size": size}}, timeout=timeout)
        r2.raise_for_status()
        rows = (r2.json().get("data", {}).get("target", {})
                or {}).get("associatedDiseases", {}).get("rows", [])
        return [{"disease": row["disease"]["name"], "score": row["score"]} for row in rows]
    except Exception:
        return []


def collect(gene: str, process: str = "T cell polarization") -> GeneEvidence:
    try:
        count, pmids = pubmed_comention(gene, process)
    except Exception:
        count, pmids = 0, []
    return GeneEvidence(
        gene=gene, pubmed_count=count, example_pmids=pmids,
        opentargets_diseases=opentargets_diseases(gene),
    )


# --------------------------------------------------------------------------- #
# Pathway / gene-set interpretation (formerly pathways.py)
# --------------------------------------------------------------------------- #
def enrich(gene_list: list[str], libraries: list[str] | None = None) -> pd.DataFrame:
    """Over-representation analysis on a nominated gene set via Enrichr (no auth).

    Returns a tidy DataFrame (term, library, adjusted p, overlap genes). Empty and
    safe if gseapy or network is unavailable — interpretation is a nicety, not a gate.
    """
    libraries = libraries or DEFAULT_LIBRARIES
    if len(gene_list) < 2:
        return pd.DataFrame(columns=["Term", "Gene_set", "Adjusted P-value", "Genes"])
    try:
        import gseapy as gp

        res = gp.enrichr(gene_list=list(gene_list), gene_sets=libraries,
                         outdir=None, no_plot=True)
        df = res.results.sort_values("Adjusted P-value").head(25)
        return df[["Term", "Gene_set", "Adjusted P-value", "Combined Score", "Genes"]]
    except Exception as exc:  # pragma: no cover
        return pd.DataFrame({"error": [str(exc)]})
