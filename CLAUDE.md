# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install for development (biocommons + cdot + pyhgvs + test deps)
pip install -e ".[dev]"

# Install biocommons-only (use with UTA or any other hdp)
pip install -e "."

# Run unit tests (fast, no external dependencies)
pytest

# Run a single test
pytest tests/unit/test_hgvs_string_validation.py::test_valid_hgvs_passes

# Run integration tests (requires local FASTA at /data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz)
pytest -m integration

# Run with coverage
pytest --cov
```

## Test Structure

```
tests/
  unit/                              # No external dependencies
    test_hgvs_string_validation.py   # HGVSConverter._hgvs_string_validation
    test_combo_checker.py            # ComboCheckerHGVSConverter with mocks
  integration/                       # Marked integration, excluded by default
    conftest.py                      # Session fixtures: biocommons_converter_grch37, pyhgvs_converter_grch37
    data/test_transcripts.json.gz    # Committed cdot transcript data for test transcripts
    test_biocommons_name_to_variant.py
    test_biocommons_variant_to_name.py
    test_pyhgvs_fields.py
```

Integration tests skip if the FASTA file is absent. pyhgvs tests skip if pyhgvs is not installed.

## Architecture

This library is a compatibility shim for migrating code from [pyhgvs (Counsyl)](https://github.com/counsyl/hgvs) to [Biocommons HGVS](https://github.com/biocommons/hgvs). It exposes a single abstract API that both libraries implement.

### Core Abstractions (`src/hgvs_shim/`)

- **`hgvs_converter.py`** — `HGVSConverter` (ABC): the primary interface. Defines operations: `create_hgvs_variant`, `normalize`, `variant_coordinate_to_g_hgvs`, `variant_coordinate_to_c_hgvs`, `hgvs_to_variant_coordinate`, `c_hgvs_remove_gene_symbol`, `get_transcript_accession`. Also defines `TranscriptInfo` dataclass (fields: `accession`, `strand`, `is_coding`, `gene_symbol`) and exception hierarchy: `HGVSException` → `HGVSNomenclatureException` (user-fixable) and `HGVSImplementationException` (library bug).

- **`hgvs_variant.py`** — `HGVSVariant` (ABC): wraps the library-specific variant object (`HGVSName` for pyhgvs, `SequenceVariant` for biocommons). Properties: `transcript`, `gene`, `kind`, `ref_allele`, `mutation_type`, `length`, `contig_accession`.

### Implementations

- **`hgvs_converter_biocommons.py`** — `BioCommonsHGVSConverter` / `BioCommonsHGVSVariant`. Constructor takes `(assembly_name, hdp)` where `hdp` is any `hgvs.dataproviders.interface.Interface` (e.g. from cdot or biocommons' UTA). Chromosome naming: biocommons uses NCBI style (`'11'`, `'MT'`); the converter maps to/from UCSC (`'chr11'`, `'chrM'`). Biocommons modern HGVS format omits ref sequence after `del`/`dup` (e.g. `c.203del` not `c.203delA`).

- **`hgvs_converter_pyhgvs.py`** — `PyHGVSConverter` / `PyHGVSVariant`. Requires pyhgvs+pysam extras (`pip install 'hgvs_shim[pyhgvs]'`). Raises `ImportError` at instantiation if not installed. Takes a `pysam.FastaFile` and optional `get_transcript` callable. `normalize()` is a no-op. Returns NCBI-style chromosome accessions from `hgvs_to_variant_coordinate`.

- **`hgvs_converter_combo.py`** — `ComboCheckerHGVSConverter`: runs multiple converters in parallel and raises `ValueError` (or logs an error) if they return different results. Used during migration to validate that both libraries agree.

### Key Design Decisions

- Variant coordinates use `(chrom, position, ref, alt)` — VCF-style, 1-based, UCSC chromosome names.
- `variant_coordinate_to_c_hgvs` takes a `TranscriptInfo` dataclass (not a bare accession string).
- `HGVSVariant.format(use_delins_for_inv=True)` converts `inv` edits to `delins` by default.
- `HGVSVariant.__eq__` compares via `format()` string output.
- `PyHGVSVariant.format(max_ref_length=N)` strips long ref alleles from del/dup/delins per HGVS spec.
- `BioCommonsHGVSConverter.__init__` accepts any `hgvs.dataproviders.interface.Interface` — the `from_local_file`/`from_rest_api` classmethods use cdot but are purely optional. The plain constructor works with biocommons' default UTA provider (`hgvs.dataproviders.uta.connect()`).
- Transcripts for integration tests are provided via a committed `test_transcripts.json.gz` (cdot format).
