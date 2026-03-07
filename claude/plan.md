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

### Hidden pyhgvs Dependencies in the Biocommons Converter

`hgvs_converter_biocommons.py` has two undeclared dependencies on pyhgvs:

**1. `from pyhgvs import get_refseq_type`** — used to detect whether a c. HGVS transcript is actually non-coding (NR_/XR_) so it can be converted as `n.` rather than `c.`. Fix by inlining using the same 3-char prefix approach pyhgvs uses (from its `REFSEQ_PREFIXES` table):

```python
_REFSEQ_TRANSCRIPT_TYPES = {
    'NM_': 'mRNA', 'XM_': 'mRNA',
    'NR_': 'RNA',  'XR_': 'RNA',
}

def _get_refseq_type(accession: str) -> str:
    return _REFSEQ_TRANSCRIPT_TYPES.get(accession[:3], None)
```

**2. `looks_like_transcript`** — called at lines 100 and 219 but never defined or imported anywhere in the shim. This is a live `NameError`. Based on pyhgvs's `REFSEQ_PREFIXES` (mRNA types: NM_, XM_; RNA types: NR_, XR_) and Ensembl (ENST prefix), the correct implementation is:

```python
def _looks_like_transcript(accession: str) -> bool:
    return accession.startswith(('NM_', 'NR_', 'XM_', 'XR_', 'ENST'))
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

The biocommons converter needs a transcript data provider (`hdp`, implementing `hgvs.dataproviders.interface.Interface`). cdot provides two implementations:

- **`JSONDataProvider(["path/to/file.json.gz"])`** — fast (500–1000 tx/sec), requires a local `.json.gz` file
- **`RESTDataProvider()`** — no local file needed, uses cdot.cc by default, slow (10 tx/sec)

For integration tests, use `RESTDataProvider()` to fetch and locally cache the specific transcripts needed, then store the result as a small `.json.gz` in `tests/integration/data/`. The transcripts needed are:
`NM_000352.3, NM_007294.3, NM_000016.4, NM_000492.3, NM_000038.5, NM_000018.3, NM_000019.3, NM_000022.2, NM_000023.2, NM_000030.2, ENST00000357654`

The cdot release JSON.gz files (~15MB for GRCh37 RefSeq) are too large to commit, but a subset file covering only these 11 transcripts will be tiny.

For pyhgvs integration tests, use `cdot.pyhgvs.pyhgvs_transcript.JSONPyHGVSTranscriptFactory` pointing at the same or a separate local JSON file (replacing the `RESTPyHGVSTranscriptFactory` currently in `conftest.py`).

---

### Parser Performance

The biocommons `Parser()` constructor was slow (~1.5s) because it compiled a parsley grammar at startup. The fix pre-generates the grammar as Python code (`hgvs.generated.hgvs_grammar.createParserClass`) so `Parser()` construction is now fast. This is merged into biocommons `main`.

The shim's current usage — `from hgvs import parser` then `parser.parse_hgvs_variant(hgvs_string)` — uses the module-level convenience function and is fine as-is. No `ParserSingleton` needed.

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

2. **Add `BioCommonsHGVSConverter` factory classmethods** so users don't need to know how to construct a cdot data provider:

```python
@classmethod
def from_local_file(cls, assembly_name: str, transcript_json_gz_path: str) -> 'BioCommonsHGVSConverter':
    """Use a local cdot JSON.gz file (fast: 500-1000 tx/sec)."""
    from cdot.hgvs.dataproviders import JSONDataProvider
    hdp = JSONDataProvider([transcript_json_gz_path])
    return cls(assembly_name, hdp)

@classmethod
def from_rest_api(cls, assembly_name: str, url: str = None) -> 'BioCommonsHGVSConverter':
    """Use the cdot REST API (default: cdot.cc). Slower but no local file needed."""
    from cdot.hgvs.dataproviders import RESTDataProvider
    hdp = RESTDataProvider(url) if url else RESTDataProvider()
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
        str(Path(__file__).parent / 'data' / 'test_transcripts.json.gz')
    )

@pytest.fixture(scope="session")
def pyhgvs_transcript_factory_grch37():
    from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory
    return JSONPyHGVSTranscriptFactory(
        str(Path(__file__).parent / 'data' / 'test_transcripts.json.gz')
    )
```

---

---

## Todo

### Phase 1: Fix Bugs ✅

- [x] `hgvs_converter_combo.py`: define `_all_equal()` helper at module level
- [x] `hgvs_converter_combo.py`: replace call to undefined `all_equal(results)` with `_all_equal(results)`
- [x] `hgvs_converter_biocommons.py`: add `_looks_like_transcript()` — checks NM_, NR_, XM_, XR_, ENST prefixes
- [x] `hgvs_converter_biocommons.py`: remove `from pyhgvs import get_refseq_type`; inline via `TranscriptInfo.is_coding`

### Phase 2: Dependencies ✅

- [x] `pyproject.toml`: remove `pyhgvs` and `pysam` from core `dependencies`
- [x] `pyproject.toml`: add `[project.optional-dependencies]` with `pyhgvs = ["pyhgvs", "pysam"]`
- [x] `pyproject.toml`: update `dev` extras to include the `pyhgvs` extra
- [x] `hgvs_converter_pyhgvs.py`: wrap top-level imports in try/except
- [x] `hgvs_converter_pyhgvs.py`: raise `ImportError` with install instructions if pyhgvs/pysam absent
- [x] `CLAUDE.md`: updated

### Phase 3: API Cleanup ✅

- [x] `hgvs_converter.py`: add `TranscriptInfo` dataclass
- [x] `hgvs_converter.py`: update abstract `variant_coordinate_to_c_hgvs` to take `TranscriptInfo`
- [x] `hgvs_converter_pyhgvs.py`: updated to accept `TranscriptInfo`
- [x] `hgvs_converter_biocommons.py`: updated to accept `TranscriptInfo`
- [x] `hgvs_converter_biocommons.py`: add `from_local_file()` and `from_rest_api()` classmethods
- [x] `hgvs_variant.py`: rename `allow_inv` → `use_delins_for_inv`
- [x] `src/hgvs_shim/__init__.py`: explicit imports and `__all__`

### Phase 4: Restructure Tests ✅

- [x] Create `tests/unit/` and `tests/integration/` directories
- [x] Generate `tests/integration/data/test_transcripts.json.gz`
- [x] Create `tests/integration/conftest.py`
- [x] Create `tests/unit/test_hgvs_string_validation.py`
- [x] Create `tests/unit/test_combo_checker.py`
- [x] Create `tests/integration/test_biocommons_name_to_variant.py` (34 tests pass)
- [x] Create `tests/integration/test_biocommons_variant_to_name.py` (34 tests pass)
- [x] Move and update `tests/integration/test_pyhgvs_fields.py`
- [x] Delete `tests/pyhgvs_tests/` and old `tests/test_pyhgvs_fields.py`
- [x] `pyproject.toml`: add pytest markers; integration tests excluded by default

Note: `test_combo_agreement.py` was removed — pyhgvs returns NCBI accessions (NC_000001.10) while biocommons returns UCSC names (chr1), making direct comparison infeasible without a full contig mapping.

### Phase 5: Cleanup ✅

- [x] Reviewed `vg_dump/` — all VG-specific Django code, nothing missing from shim
- [x] Deleted `vg_dump/`
- [x] Updated `CLAUDE.md`

---

## Not In Scope

- `CHGVS` / `VariantCoordinate` utility classes — VG-specific, not needed for migration
- `clean_hgvs` — input cleaning is the caller's responsibility
- `HgvsMatchRefAllele` / `HgvsOriginallyNormalized` — application-level concerns, not shim concerns
- `ParserSingleton` — fixed upstream in biocommons
- Anything else in `vg_dump/` not already present in `src/hgvs_shim/`
