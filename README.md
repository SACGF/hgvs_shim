# hgvs_shim

[![pypi_rel](https://img.shields.io/pypi/v/hgvs_shim.svg)](https://pypi.org/project/hgvs_shim/)

hgvs-shim is a small compatibility layer that simplifies migrating code from the [Counsyl pyhgvs](https://github.com/counsyl/hgvs) library to [Biocommons HGVS](https://github.com/biocommons/hgvs).

# Background

There are 2 Python HGVS libraries:

- [pyhgvs](https://github.com/counsyl/hgvs) — simpler, but appears abandoned (no activity for years)
- [Biocommons HGVS](https://github.com/biocommons/hgvs) — actively developed, more features (alignment gaps, inversions, uncertain coordinates)

# Motivation

In the VariantGrid project, we initially started with pyhgvs. We eventually decided Biocommons HGVS was worth the added complexity and wrote [cdot](https://github.com/SACGF/cdot) for providing transcripts to both libraries.

Migrating was risky: pyhgvs classes were spread throughout our codebase, and we needed to validate that Biocommons would produce identical results before switching. The transition steps were:

1. Abstract away pyhgvs classes and methods behind our own interfaces
2. Write implementations for both pyhgvs and Biocommons
3. Run both in parallel during testing, raising errors if they disagreed
4. Switch to Biocommons-only once confident

hgvs_shim is a cleaned-up library from this work.

# Installation

```bash
# Biocommons only (pair with any hdp: UTA, cdot, etc.)
pip install hgvs_shim

# With pyhgvs support (for running both backends side-by-side during migration)
pip install 'hgvs_shim[pyhgvs]'
```

# Usage

## Setting up the converter

`BioCommonsHGVSConverter` takes any biocommons data provider (`hgvs.dataproviders.interface.Interface`). The simplest option is the default biocommons UTA database connection:

```python
import hgvs.dataproviders.uta as uta
from hgvs_shim import BioCommonsHGVSConverter

hdp = uta.connect()  # connects to biocommons' public UTA PostgreSQL instance
converter = BioCommonsHGVSConverter('GRCh38', hdp)
```

For faster or offline transcript data, [cdot](https://github.com/SACGF/cdot) provides JSON.gz files and a REST API (`pip install 'hgvs_shim[cdot]'`):

```python
from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider
from hgvs_shim import BioCommonsHGVSConverter

# Local JSON.gz file (fast: 500-1000 tx/sec)
hdp = JSONDataProvider(['/path/to/cdot-0.2.x.refseq.grch37.json.gz'])
converter = BioCommonsHGVSConverter('GRCh37', hdp)

# cdot.cc REST API (no local file needed, ~10 tx/sec)
hdp = RESTDataProvider()
converter = BioCommonsHGVSConverter('GRCh37', hdp)
```

## HGVS string → variant coordinate

```python
chrom, pos, ref, alt = converter.hgvs_to_variant_coordinate('NM_000492.3:c.1521_1523delCTT')
# ('chr7', 117548628, 'ACTT', 'A')
```

## Variant coordinate → c.HGVS

```python
from hgvs_shim import TranscriptInfo

transcript_info = TranscriptInfo(
    accession='NM_000492.3',
    strand='+',
    is_coding=True,
    gene_symbol='CFTR',  # optional
)

variant = converter.variant_coordinate_to_c_hgvs('chr7', 117548628, 'ACTT', 'A', transcript_info)
print(variant.format())
# NM_000492.3:c.1521_1523del
```

## Normalize a variant

```python
variant = converter.create_hgvs_variant('NM_000492.3:c.1521_1523delCTT')
normalized = converter.normalize(variant)
print(normalized.format())
```

## Running both backends in parallel (migration validation)

```python
from cdot.hgvs.dataproviders import JSONDataProvider
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory
from pysam import FastaFile
from hgvs_shim import BioCommonsHGVSConverter, PyHGVSConverter, ComboCheckerHGVSConverter

hdp = JSONDataProvider(['/path/to/cdot.json.gz'])
biocommons = BioCommonsHGVSConverter('GRCh37', hdp)

factory = JSONPyHGVSTranscriptFactory(['/path/to/cdot.json.gz'])
pyhgvs_conv = PyHGVSConverter(FastaFile('/path/to/GRCh37.fa'), factory.get_transcript_grch37)

# Raises ValueError if the two converters return different results
combo = ComboCheckerHGVSConverter([biocommons, pyhgvs_conv], die_on_error=True)
result = combo.hgvs_to_variant_coordinate('NM_000352.3:c.215A>G')
```

## Exception handling

```python
from hgvs_shim import HGVSNomenclatureException, HGVSImplementationException

try:
    result = converter.hgvs_to_variant_coordinate(hgvs_string)
except HGVSNomenclatureException:
    # Bad HGVS string — user error, fixable
    ...
except HGVSImplementationException:
    # Library failure — not user-fixable
    ...
```

## Format differences vs pyhgvs

Biocommons follows the modern HGVS specification, which differs from pyhgvs in two ways:

**del/dup notation** — biocommons omits the deleted/duplicated sequence:
```
pyhgvs:     NM_000492.3:c.442delA
biocommons: NM_000492.3:c.442del
```

**delins VCF representation** — biocommons uses minimal representation (no anchor base):
```
pyhgvs:     ('chr7', 117171119, 'CA', 'C')   # anchor base included
biocommons: ('chr7', 117171120, 'A',  'C')   # no anchor base
```
