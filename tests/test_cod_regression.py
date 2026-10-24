"""Real COD CIFs (tests/data/cod/) against their stored output (expected.json).

A failure means a change altered the output for a real CIF: if the change is intended, run dev/update_expected.py,
read which fields moved, and commit the new expected.json with it.
"""

import json
from pathlib import Path

import pytest

from cod_summary import differences, summarize

DATA = Path(__file__).parent / "data" / "cod"
EXPECTED = json.loads((DATA / "expected.json").read_text(encoding="utf-8"))


def test_every_cif_has_an_expected_entry():
    assert sorted(EXPECTED) == sorted(p.stem for p in DATA.glob("*.cif"))


@pytest.mark.parametrize("cod_id", sorted(EXPECTED))
def test_cod_cif(cod_id):
    diff = differences(EXPECTED[cod_id], summarize(DATA / f"{cod_id}.cif"))
    assert not diff, "\n".join(f"{field}: {text}" for field, text in diff.items())
