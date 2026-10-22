#!/usr/bin/env python3
"""Snapshot OmniPath ligand-receptor pairs for the offline biology check.

Keeps interactions whose two partners are both annotated as an intercellular
ligand or receptor, from the ``omnipath`` and ``ligrecextra`` datasets,
restricted to resources that allow commercial use (``license=commercial``).
Complexes are split into their members. Writes ``data/knowledge/omnipath_lr.tsv.gz``
with columns a, b, directed, stimulation, inhibition, n_sources, and
``omnipath_roles.tsv.gz`` marking each gene as ligand and/or receptor.

    python scripts/fetch_omnipath.py
"""

import argparse
import csv
import datetime
import gzip
import io
import json
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BASE = "https://omnipathdb.org"
INTERACTIONS = (f"{BASE}/interactions?datasets=omnipath,ligrecextra&genesymbols=yes"
                "&license=commercial&fields=sources&format=tsv")
INTERCELL = f"{BASE}/intercell?categories=ligand,receptor&scope=generic&license=commercial&format=tsv"


def fetch_tsv(url):
    """Rows of a TSV endpoint as dicts."""
    with urllib.request.urlopen(url, timeout=600) as r:
        return list(csv.DictReader(io.StringIO(r.read().decode()), delimiter="\t"))


def members(symbol):
    """Gene symbols of a complex ('TGFBR1_TGFBR2') or a single protein."""
    return [s for s in symbol.split("_") if s]


def lr_pairs(interactions, intercell_genes):
    """Ligand-receptor pairs with signs, merged over complex members and duplicate records."""
    pairs = {}
    for row in interactions:
        for a in members(row["source_genesymbol"]):
            for b in members(row["target_genesymbol"]):
                if a == b or a not in intercell_genes or b not in intercell_genes:
                    continue
                rec = pairs.setdefault((a, b), {"directed": False, "stimulation": False,
                                                "inhibition": False, "sources": set()})
                for flag in ("directed", "stimulation", "inhibition"):
                    rec[flag] |= row[f"is_{flag}"] == "True"
                rec["sources"].update(s.split("_")[0] for s in row["sources"].split(";") if s)
    return pairs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(REPO / "data/knowledge/omnipath_lr.tsv.gz"))
    args = ap.parse_args(argv)

    intercell = fetch_tsv(INTERCELL)
    genes = {r["genesymbol"] for r in intercell}
    pairs = lr_pairs(fetch_tsv(INTERACTIONS), genes)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["a", "b", "directed", "stimulation", "inhibition", "n_sources"])
        for (a, b), rec in sorted(pairs.items()):
            w.writerow([a, b, int(rec["directed"]), int(rec["stimulation"]), int(rec["inhibition"]),
                        len(rec["sources"])])
    roles = {}
    for r in intercell:
        roles.setdefault(r["genesymbol"], set()).add(r["category"])
    with gzip.open(out.with_name("omnipath_roles.tsv.gz"), "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "ligand", "receptor"])
        for g, cats in sorted(roles.items()):
            w.writerow([g, int("ligand" in cats), int("receptor" in cats)])
    meta = {"fetched": datetime.date.today().isoformat(), "interactions": INTERACTIONS,
            "intercell": INTERCELL, "pairs": len(pairs), "intercell_genes": len(genes)}
    out.with_name("omnipath_lr.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"{len(pairs)} ligand-receptor pairs over {len(genes)} intercell genes -> {out}")


if __name__ == "__main__":
    main()
