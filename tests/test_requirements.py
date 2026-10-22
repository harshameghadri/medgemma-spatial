"""Training pins must match the runtime lock wherever they overlap."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIN = re.compile(r"^([A-Za-z0-9_.\-]+)==([^\s;]+)")


def pins(path):
    out = {}
    for line in (ROOT / path).read_text().splitlines():
        m = PIN.match(line.strip())
        if m:
            out[m.group(1).lower()] = m.group(2)
    return out


def test_training_pins_match_lock():
    lock, train = pins("requirements.txt"), pins("requirements-train.txt")
    mismatched = {k: (train[k], lock[k]) for k in train.keys() & lock.keys() if train[k] != lock[k]}
    assert not mismatched, f"requirements-train.txt vs requirements.txt: {mismatched}"
