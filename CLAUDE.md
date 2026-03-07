# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install for development (includes pysam and cdot)
pip install -e ".[dev]"

# Run all tests
pytest

# Run a single test
pytest tests/test_pyhgvs_fields.py::test_variant_coordinate_to_g

# Run with coverage
pytest --cov
```

**Note:** Tests in `tests/test_pyhgvs_fields.py` require an external FASTA file at `/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz` and a running cdot REST service. The `tests/pyhgvs_tests/` subdirectory uses local mock genome files and runs without external data.

## Project arrangement

vg_dump - a dump of hgvs related source control from VariantGrid if you want to see original implementation
tests/pyhgvs - a dump of pyhgvs tests - a goal should be to get ALL of the pyhgvs tests to run, but with the biocommons implementation

## Architecture

This library is a compatibility shim for migrating code from [pyhgvs (Counsyl)](https://github.com/counsyl/hgvs) to [Biocommons HGVS](https://github.com/biocommons/hgvs). It exposes a single abstract API that both libraries implement.

### Core Abstractions (`src/hgvs_shim/`)

- **`hgvs_converter.py`** — `HGVSConverter` (ABC): the primary interface. Defines operations: `create_hgvs_variant`, `normalize`, `variant_coordinate_to_g_hgvs`, `variant_coordinate_to_c_hgvs`, `hgvs_to_variant_coordinate`, `c_hgvs_remove_gene_symbol`, `get_transcript_accession`. Also defines the exception hierarchy: `HGVSException` → `HGVSNomenclatureException` (user-fixable) and `HGVSImplementationException` (library bug).

- **`hgvs_variant.py`** — `HGVSVariant` (ABC): wraps the library-specific variant object (`HGVSName` for pyhgvs, `SequenceVariant` for biocommons). Properties: `transcript`, `gene`, `kind`, `ref_allele`, `mutation_type`, `length`, `contig_accession`.

### Implementations

- **`hgvs_converter_pyhgvs.py`** — `PyHGVSConverter` / `PyHGVSVariant`: wraps pyhgvs. Requires a `pysam.FastaFile` and an optional `get_transcript` callable for c.HGVS conversions. `normalize()` is a no-op (pyhgvs limitation).

- **`hgvs_converter_biocommons.py`** — `BioCommonsHGVSConverter` / `BioCommonsHGVSVariant`: wraps biocommons hgvs. Requires an `assembly_name` and a `hdp` (data provider implementing `hgvs.dataproviders.interface.Interface`, typically supplied via cdot). Maps biocommons exceptions to the shim's exception hierarchy.

- **`hgvs_converter_combo.py`** — `ComboCheckerHGVSConverter`: runs multiple converters in parallel and raises `ValueError` (or logs an error) if they return different results. Used during migration to validate that both libraries agree.

### Key Design Decisions

- Variant coordinates use `(chrom, position, ref, alt)` — VCF-style, 1-based.
- `HGVSVariant.__eq__` compares via `format()` string output.
- `PyHGVSVariant.format(max_ref_length=N)` strips long ref alleles from del/dup/delins per HGVS spec.
- `BioCommonsHGVSConverter` converts `inv` edits to `delins` on format (unless `allow_inv=True`).
- Transcripts for biocommons are provided via [cdot](https://github.com/SACGF/cdot); pyhgvs transcripts via `cdot.pyhgvs.pyhgvs_transcript.RESTPyHGVSTranscriptFactory`.
