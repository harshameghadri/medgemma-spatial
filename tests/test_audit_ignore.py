"""Accepted pip-audit advisories must stay scoped to packages still pinned below the fix."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_pillow_ignores_only_while_pillow_below_13():
    entries = [line.split() for line in (ROOT / "audit-ignore.txt").read_text().splitlines()
               if re.match(r"(GHSA|PYSEC|CVE)-", line)]
    assert all(len(e) >= 2 for e in entries)
    pinned = re.search(r"^pillow==(\d+)", (ROOT / "requirements.txt").read_text(), re.M)
    if pinned and int(pinned.group(1)) >= 13:
        assert not [e for e in entries if e[1] == "pillow"], "pillow fixed: drop its audit-ignore entries"
