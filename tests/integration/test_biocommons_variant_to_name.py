"""
Integration tests: variant coordinate → HGVS name (via variant_coordinate_to_c_hgvs).
Ported from pyhgvs test_hgvs_names.py _name_variants (name_canonical=True entries only).
"""
import pytest

from hgvs_shim.hgvs_converter import TranscriptInfo

pytestmark = pytest.mark.integration

# TranscriptInfo for each tested transcript (strand from pyhgvs _refgene data)
_TRANSCRIPT_INFO = {
    'NM_000016.4': TranscriptInfo('NM_000016.4', '+', True),
    'NM_000018.3': TranscriptInfo('NM_000018.3', '+', True),
    'NM_000019.3': TranscriptInfo('NM_000019.3', '+', True),
    'NM_000022.2': TranscriptInfo('NM_000022.2', '-', True),
    'NM_000023.2': TranscriptInfo('NM_000023.2', '+', True),
    'NM_000030.2': TranscriptInfo('NM_000030.2', '+', True),
    'NM_000038.5': TranscriptInfo('NM_000038.5', '+', True),
    'NM_000352.3': TranscriptInfo('NM_000352.3', '-', True),
    'NM_000492.3': TranscriptInfo('NM_000492.3', '+', True),
    'NM_007294.3': TranscriptInfo('NM_007294.3', '-', True),
}

# (expected_hgvs_name, variant_coordinate) for name_canonical=True entries
_variant_to_name = [
    # Simple SNPs
    ('NM_000352.3:c.215A>G',            ('chr11', 17496508, 'T', 'C')),
    ('NM_000352.3:c.72C>A',             ('chr11', 17498252, 'G', 'T')),
    ('NM_000352.3:c.3885C>G',           ('chr11', 17418843, 'G', 'C')),

    # SNPs within introns
    ('NM_000352.3:c.1630+1G>A',         ('chr11', 17464266, 'C', 'T')),
    ('NM_000352.3:c.1630+1G>C',         ('chr11', 17464266, 'C', 'G')),
    ('NM_000352.3:c.1630+1G>T',         ('chr11', 17464266, 'C', 'A')),
    ('NM_000352.3:c.1672-20A>G',        ('chr11', 17452526, 'T', 'C')),
    ('NM_000352.3:c.1672-74G>A',        ('chr11', 17452580, 'C', 'T')),
    ('NM_000352.3:c.1923+5G>T',         ('chr11', 17450107, 'C', 'A')),
    ('NM_000352.3:c.2041-21G>A',        ('chr11', 17449510, 'C', 'T')),
    ('NM_000352.3:c.2116+3A>G',         ('chr11', 17449411, 'T', 'C')),
    ('NM_000352.3:c.2116+1G>T',         ('chr11', 17449413, 'C', 'A')),
    ('NM_000352.3:c.2116+2T>C',         ('chr11', 17449412, 'A', 'G')),

    # Indels — biocommons uses delins notation without ref sequence
    ('NM_000018.3:c.1077_1077+1delinsCAC',    ('chr17', 7126183, 'CGG', 'CCAC')),
    ('NM_000023.2:c.585-2_585-1delinsT',       ('chr17', 48246450, 'CAG', 'CT')),
    ('NM_000030.2:c.2_3delinsAT',              ('chr2', 241808283, 'ATG', 'AAT')),

    # Single letter del and insert
    ('NM_000016.4:c.945+4delinsGC',    ('chr1', 76216234, 'AA', 'AGC')),
    ('NM_000038.5:c.1617delinsGAA',    ('chr5', 112163693, 'AC', 'AGAA')),
    ('NM_000038.5:c.4256delinsCC',     ('chr5', 112175546, 'AG', 'ACC')),

    # Delete region — biocommons omits ref sequence after del
    ('NM_000016.4:c.430_432del',        ('chr1', 76200511, 'AAAG', 'A')),

    # Single letter insert, delete, duplication — biocommons omits ref after del/dup
    ('NM_000016.4:c.203del',            ('chr1', 76198412, 'GA', 'G')),
    ('NM_000016.4:c.244dup',            ('chr1', 76198564, 'C', 'CT')),
    ('NM_000016.4:c.387+1del',          ('chr1', 76199309, 'TG', 'T')),
    ('NM_000016.4:c.475del',            ('chr1', 76205669, 'AT', 'A')),
    ('NM_000016.4:c.1189dup',           ('chr1', 76227049, 'C', 'CT')),
    ('NM_000016.4:c.1191del',           ('chr1', 76227051, 'AT', 'A')),
    ('NM_000016.4:c.306_307insG',       ('chr1', 76199232, 'T', 'TG')),

    # Alignment tests
    ('NM_000492.3:c.935_937del',        ('chr7', 117180210, 'CCTT', 'C')),
    ('NM_000492.3:c.442del',            ('chr7', 117171120, 'CA', 'C')),
    ('NM_000492.3:c.805_806del',        ('chr7', 117176660, 'AAT', 'A')),
    ('NM_000492.3:c.1155_1156dup',      ('chr7', 117182104, 'A', 'AAT')),
    ('NM_000492.3:c.3889dup',           ('chr7', 117292905, 'A', 'AT')),

    # After stop codon
    ('NM_000492.3:c.*3A>C',             ('chr7', 117307165, 'A', 'C')),

    # Non-canonical variant (name_canonical=True, var_canonical=False)
    ('NM_000492.3:c.1210-7_1210-6dup',  ('chr7', 117188682, 'GTT', 'GTTTT')),
]


def _get_transcript_accession(hgvs_name: str) -> str:
    return hgvs_name.split(':')[0]


@pytest.mark.parametrize("expected_name,variant", _variant_to_name)
def test_variant_to_name(biocommons_converter_grch37, expected_name, variant):
    chrom, position, ref, alt = variant
    accession = _get_transcript_accession(expected_name)
    transcript_info = _TRANSCRIPT_INFO[accession]
    result = biocommons_converter_grch37.variant_coordinate_to_c_hgvs(
        chrom, position, ref, alt, transcript_info
    )
    assert result.format() == expected_name
