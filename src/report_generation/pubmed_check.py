"""Check report citations against PubMed (NCBI E-utilities, no key required).

A citation is *verified* only when first author, publication year and journal
all match a PubMed record; surname + year alone matches almost anything. Any
network failure yields ``unchecked`` instead of an exception, so callers work
offline. Results are cached in memory and, optionally, in a JSON file.

Statuses: verified | mismatch (author+year exist, not in that journal) |
not_found | incomplete (no journal given, cannot verify) | unchecked.
"""

import json
import os
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Callable, Dict, List, Optional

from src.report_generation.faithfulness import (
    HEADING,
    extract_citations,
    score_report,
    split_sentences,
    unsupported_cell_types,
    unsupported_numbers,
)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "medgemma-spatial"
UNVERIFIED = {"mismatch", "not_found", "incomplete"}


def _http_get(url: str, timeout: float) -> dict:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


class PubMedChecker:
    """Rate-limited, cached PubMed lookups; inject ``fetch`` in tests."""

    def __init__(self, cache_path: Optional[str] = None, fetch: Optional[Callable] = None,
                 timeout: float = 6.0, api_key: Optional[str] = None):
        self.cache_path = Path(cache_path) if cache_path else None
        self.cache: Dict[str, dict] = {}
        if self.cache_path and self.cache_path.exists():
            try:
                self.cache = json.loads(self.cache_path.read_text())
            except (OSError, ValueError):
                self.cache = {}
        self.fetch = fetch or _http_get
        self.timeout = timeout
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self.min_interval = 0.11 if self.api_key else 0.4
        self.retries = 3
        self._last = 0.0
        self._lock = threading.Lock()

    def _get(self, endpoint: str, **params) -> dict:
        params.update(retmode="json", tool=TOOL)
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS}{endpoint}?{urllib.parse.urlencode(params)}"
        for attempt in range(self.retries):
            with self._lock:
                wait = self.min_interval - (time.monotonic() - self._last)
                if wait > 0:
                    time.sleep(wait)
                self._last = time.monotonic()
            try:
                return self.fetch(url, self.timeout)
            except urllib.error.HTTPError as e:
                if e.code not in (429, 500, 502, 503) or attempt == self.retries - 1:
                    raise
                time.sleep(1.0 * (attempt + 1))

    def search(self, term: str) -> List[str]:
        """PMIDs for an E-utilities term (cached)."""
        if term in self.cache:
            return self.cache[term]["ids"]
        ids = self._get("esearch.fcgi", db="pubmed", term=term, retmax=5)["esearchresult"]["idlist"]
        self.cache[term] = {"ids": ids}
        return ids

    def summary(self, pmid: str) -> dict:
        """Title / journal / year for one PMID (best effort)."""
        key = f"summary:{pmid}"
        if key not in self.cache:
            rec = self._get("esummary.fcgi", db="pubmed", id=pmid)["result"].get(pmid, {})
            self.cache[key] = {"title": rec.get("title"), "journal": rec.get("fulljournalname"),
                               "year": (rec.get("pubdate") or "")[:4]}
        return self.cache[key]

    def verify(self, surname: str, year: int, journal: Optional[str]) -> dict:
        """Status of one author-year-journal citation."""
        base = f"{surname}[1au] AND {year}[pdat]"
        try:
            if journal:
                ids = self.search(f'{base} AND "{journal}"[jour]')
                if ids:
                    return {"status": "verified", "pmid": ids[0], **self.summary(ids[0])}
            loose = self.search(base)
        except Exception as e:  # network, HTTP or malformed payload
            return {"status": "unchecked", "error": f"{type(e).__name__}: {e}"}
        if not loose:
            return {"status": "not_found"}
        return {"status": "mismatch" if journal else "incomplete", "candidates": loose[:3]}

    def save(self) -> None:
        if self.cache_path:
            try:
                self.cache_path.parent.mkdir(parents=True, exist_ok=True)
                self.cache_path.write_text(json.dumps(self.cache))
            except OSError:
                pass


def check_citations(text: str, checker: Optional[PubMedChecker] = None) -> List[dict]:
    """Verify every author-year citation found in a report."""
    checker = checker or PubMedChecker(os.environ.get("MEDGEMMA_PUBMED_CACHE"))
    results = [{**c, **checker.verify(c["surname"], c["year"], c["journal"])} for c in extract_citations(text)]
    checker.save()
    return results


def strip_unverified(text: str, results: List[dict]) -> str:
    """Remove sentences carrying a citation that PubMed could not confirm."""
    bad = [r["text"] for r in results if r["status"] in UNVERIFIED]
    if not bad:
        return text
    out = []
    for para in text.split("\n"):
        kept = [s for s in split_sentences(para) if not any(b in s for b in bad)] if para.strip() else []
        out.append(" ".join(kept) if para.strip() else "")
    return "\n".join(out)


def strip_ungrounded(text: str, prompt: str) -> str:
    """Remove sentences stating numbers or cell types that the prompt did not provide."""
    out = []
    for para in text.split("\n"):
        if not para.strip():
            out.append("")
            continue
        kept = [s for s in split_sentences(para)
                if not unsupported_numbers(s, prompt) and not unsupported_cell_types(s, prompt)]
        out.append(" ".join(kept))
    return _drop_empty_headings("\n".join(out))


def _drop_empty_headings(text: str) -> str:
    """Remove a bold heading line whose section has no text left."""
    lines = text.split("\n")
    keep = []
    for i, line in enumerate(lines):
        if HEADING.match(line):
            nxt = next((x for x in lines[i + 1:] if x.strip()), None)
            if nxt is None or HEADING.match(nxt):
                continue
        keep.append(line)
    return re.sub(r"\n{3,}", "\n\n", "\n".join(keep)).strip()


def review_report(text: str, prompt: str, verify_citations: bool = True,
                  checker: Optional[PubMedChecker] = None) -> dict:
    """Grounding review shown next to a generated report.

    Returns the score, PubMed results and ``clean_text``: the report without
    sentences that state numbers or cell types absent from the prompt, or cite
    papers PubMed could not confirm. With ``verify_citations=False`` (or
    offline) citations are left in place and reported as unchecked.
    """
    score = score_report(text, prompt)
    grounded = strip_ungrounded(text, prompt) if prompt else text
    if not score["citations"]:
        return {**score, "pubmed": [], "clean_text": grounded}
    if not verify_citations:
        results = [{**c, "status": "unchecked"} for c in extract_citations(text)]
    else:
        results = check_citations(text, checker)
    return {**score, "pubmed": results, "clean_text": strip_unverified(grounded, results)}
