# Plan

## Purpose and Scope

This is a **minimalist shim** to help migrate codebases from pyhgvs to biocommons HGVS. Nothing more.

`vg_dump/` is a code reference to check in case something was missed when extracting from VariantGrid. It is NOT a feature target — do not implement everything in there. Delete it once we're confident nothing has been missed.

---

## Analysis

### The Core Approach is Right

The abstract `HGVSConverter`/`HGVSVariant` pattern is the correct architecture. The VG dump confirms it worked in production. Keep it.

---

### Active Bug: `all_equal` is Undefined

`hgvs_converter_combo.py:42` calls `all_equal(results)` but this function is never defined or imported. This is a `NameError` at runtime. Fix with a simple local helper:

```python
def _all_equal(iterable):
    items = [x for x in iterable if not isinstance(x, Exception)]
    return len(set(items)) <= 1
```

---

### Hidden pyhgvs Dependency in the Biocommons Converter

`hgvs_converter_biocommons.py` imports `from pyhgvs import get_refseq_type` to detect RNA transcripts. This means the biocommons-only path has an undeclared pyhgvs dependency. Fix by inlining:

```python
def _get_refseq_type(transcript_accession: str) -> str:
    if transcript_accession.startswith('NR_'):
        return 'RNA'
    if transcript_accession.startswith('NM_'):
        return 'mRNA'
    return 'unknown'
```

---

### `variant_coordinate_to_c_hgvs` Signature Mismatch

`PyHGVSConverter` takes `transcript_accession: str`; `BioCommonsHGVSConverter` takes a `transcript_version` object with `.accession`, `.strand`, `.is_coding`, `.gene_symbol`. The biocommons converter needs strand to correctly normalize intronic variants on minus-strand genes before the g→c mapping.

Fix: define a `TranscriptInfo` dataclass in `hgvs_converter.py` that both implementations accept. The pyhgvs implementation just uses `.accession`.

```python
# src/hgvs_shim/hgvs_converter.py
from dataclasses import dataclass

@dataclass
class TranscriptInfo:
    accession: str           # e.g. "NM_001637.3"
    strand: str              # '+' or '-'
    is_coding: bool          # True for NM_, False for NR_
    gene_symbol: str | None = None
```

---

### The pyhgvs Tests Cannot Currently Run

`tests/pyhgvs_tests/` uses relative imports expecting to be inside the pyhgvs package:

```python
from .. import CDNACoord, HGVSName, parse_hgvs_name, format_hgvs_name
from ..utils import read_transcripts
from ..variants import normalize_variant
```

There are three levels of things tested:

| Level | What | Action |
|-------|------|--------|
| 1 | `CDNACoord` parsing, `genomic_to_cdna_coord` | pyhgvs-internal types with no biocommons equivalent. **Skip.** |
| 2 | `normalize_variant` | pyhgvs VCF normalizer, not equivalent to biocommons normalization. **Skip.** |
| 3 | Round-trip `name↔variant` (`_name_variants` data) | Exactly what the shim should test. **Port the data.** |

The right move: **delete `tests/pyhgvs_tests/`** and port the `_name_variants` test data into new pytest-parametrized tests against `BioCommonsHGVSConverter`. The test data is the valuable thing; the test infrastructure is not.

---

### Test Infrastructure: Transcript Data for Biocommons

The biocommons converter needs a transcript data provider (`hdp`, implementing `hgvs.dataproviders.interface.Interface`). Use `cdot`'s `JSONDataProvider` pointing at a local JSON file.

Fetch the required transcripts from the cdot REST API (e.g. `https://cdot.cc/transcript/NM_000352.3`) and store the result locally in `tests/data/`. The transcripts needed are:
`NM_000352.3, NM_007294.3, NM_000016.4, NM_000492.3, NM_000038.5, NM_000018.3, NM_000019.3, NM_000022.2, NM_000023.2, NM_000030.2, ENST00000357654`

This JSON will be small and committable.

---

### Parser Performance

The biocommons `Parser()` constructor was slow (~1.5s) because it compiled a parsley grammar at startup. A fix was merged into biocommons `main` and will be in the next release. Installing from GitHub gives correct performance. No `ParserSingleton` needed in the shim.

---

## Plan

### Phase 1: Fix Bugs

1. Add `_all_equal` helper to `hgvs_converter_combo.py` and use it.
2. Remove `from pyhgvs import get_refseq_type` from `hgvs_converter_biocommons.py`; inline the logic.

### Phase 2: Dependencies

Make `pyhgvs` and `pysam` optional. This gives two meaningful install targets:

```toml
[project.optional-dependencies]
pyhgvs = ["pyhgvs", "pysam"]   # full migration layer: run both backends side-by-side
# default install (no extra) = biocommons only: for codebases that have already migrated
```

`PyHGVSConverter` should raise a clear `ImportError` at instantiation time if pyhgvs/pysam are not installed, not at import time.

### Phase 3: API Cleanup

1. **Add `TranscriptInfo` dataclass** to `hgvs_converter.py` and update both converter implementations to accept it for `variant_coordinate_to_c_hgvs`.

2. **Add `BioCommonsHGVSConverter.from_local_file()` classmethod** so users don't need to know how to construct a cdot `JSONDataProvider`:

```python
@classmethod
def from_local_file(cls, assembly_name: str, transcript_json_path: str) -> 'BioCommonsHGVSConverter':
    from cdot.hgvs.dataproviders import JSONDataProvider
    hdp = JSONDataProvider([transcript_json_path])
    return cls(assembly_name, hdp)
```

3. **Fix `format()` parameter naming** on `HGVSVariant`. Currently `allow_inv=False` means "convert inv to delins", which is backwards-reading. Rename to `use_delins_for_inv: bool = True`.

4. **Clean up `__init__.py`** to export the full public API explicitly:

```python
from hgvs_shim.hgvs_converter import (
    HGVSConverter, HGVSConverterType, TranscriptInfo,
    HGVSException, HGVSNomenclatureException, HGVSImplementationException,
)
from hgvs_shim.hgvs_variant import HGVSVariant
from hgvs_shim.hgvs_converter_biocommons import BioCommonsHGVSConverter
from hgvs_shim.hgvs_converter_pyhgvs import PyHGVSConverter
from hgvs_shim.hgvs_converter_combo import ComboCheckerHGVSConverter
```

### Phase 4: Restructure Tests

Delete `tests/pyhgvs_tests/`. New structure:

```
tests/
  conftest.py                          # top-level: marks, shared config
  unit/
    test_hgvs_string_validation.py     # HGVSConverter._hgvs_string_validation — no external deps
    test_combo_checker.py              # ComboCheckerHGVSConverter logic with mocks
  integration/
    conftest.py                        # fixtures: biocommons_converter_grch37, pyhgvs_converter_grch37
    data/
      test_transcripts.json            # cdot JSON for the test transcripts (committed)
    test_biocommons_name_to_variant.py # ported _name_variants from pyhgvs test suite
    test_biocommons_variant_to_name.py # ported _name_variants (reverse direction)
    test_pyhgvs_fields.py              # existing test, moved here
    test_combo_agreement.py            # run both converters and assert they agree
```

Example ported test:

```python
# tests/integration/test_biocommons_name_to_variant.py
import pytest

# Ported from pyhgvs test_hgvs_names.py _name_variants (name_canonical=True, var_canonical=True only)
_name_variants = [
    ('NM_000352.3:c.215A>G',            ('chr11', 17496508, 'T', 'C')),
    ('NM_000352.3:c.1630+1G>A',         ('chr11', 17464266, 'C', 'T')),  # intronic
    ('NM_000492.3:c.935_937delTCT',     ('chr7',  117180210, 'CCTT', 'C')),  # alignment test
    ('NM_000492.3:c.1155_1156dupTA',    ('chr7',  117182104, 'A', 'AAT')),
    # ... all cases where name_canonical=True and var_canonical=True
]

@pytest.mark.parametrize("hgvs_name,expected", _name_variants)
def test_name_to_variant(biocommons_converter_grch37, hgvs_name, expected):
    result = biocommons_converter_grch37.hgvs_to_variant_coordinate(hgvs_name)
    assert result == expected
```

The `integration/conftest.py` fixture:

```python
@pytest.fixture(scope="session")
def biocommons_converter_grch37():
    return BioCommonsHGVSConverter.from_local_file(
        'GRCh37',
        str(Path(__file__).parent / 'data' / 'test_transcripts.json')
    )
```

---

## Not In Scope

- `CHGVS` / `VariantCoordinate` utility classes — VG-specific, not needed for migration
- `clean_hgvs` — input cleaning is the caller's responsibility
- `HgvsMatchRefAllele` / `HgvsOriginallyNormalized` — application-level concerns, not shim concerns
- `ParserSingleton` — fixed upstream in biocommons
- Anything else in `vg_dump/` not already present in `src/hgvs_shim/`
