"""Notebooks must be valid nbformat and must not embed credentials."""

import re
from pathlib import Path

import nbformat
import pytest

NOTEBOOKS = sorted(Path("notebooks").glob("*.ipynb"))
SECRET_PATTERNS = [
    re.compile(r"hf_[A-Za-z0-9]{30,}"),
    re.compile(r"sk-ant-[A-Za-z0-9_-]{20,}"),
]


@pytest.mark.parametrize("path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_is_valid(path):
    nb = nbformat.read(path, as_version=4)
    nbformat.validate(nb)


@pytest.mark.parametrize("path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_has_no_secrets(path):
    text = path.read_text(encoding="utf-8")
    for pattern in SECRET_PATTERNS:
        assert not pattern.search(text), f"{path} contains a credential-like string"
