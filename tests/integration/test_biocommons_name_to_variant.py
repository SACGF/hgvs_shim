"""
Integration tests: HGVS name → variant coordinate (via hgvs_to_variant_coordinate).
Ported from pyhgvs test_hgvs_names.py _name_variants (var_canonical=True entries only).
"""
import pytest

pytestmark = pytest.mark.integration

# (hgvs_name, expected_variant) for var_canonical=True entries
# Excludes ENST00000357654 (not in cdot.cc REST data)
_name_to_variant = [
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

    # Indels (name_canonical=True, var_canonical=True)
    # Note: biocommons uses minimal VCF representation (no anchor base) vs pyhgvs which adds one
    ('NM_000018.3:c.1077_1077+1delGGinsCAC',        ('chr17', 7126184, 'GG', 'CAC')),
    ('NM_000023.2:c.585-2_585-1delAGinsT',           ('chr17', 48246451, 'AG', 'T')),
    ('NM_000030.2:c.2_3delTGinsAT',                  ('chr2', 241808284, 'TG', 'AT')),

    # Single letter del and insert
    ('NM_000016.4:c.945+4delAinsGC',    ('chr1', 76216235, 'A', 'GC')),
    ('NM_000038.5:c.1617delCinsGAA',    ('chr5', 112163694, 'C', 'GAA')),
    ('NM_000038.5:c.4256delGinsCC',     ('chr5', 112175547, 'G', 'CC')),

    # Delete region
    ('NM_000016.4:c.430_432delAAG',     ('chr1', 76200511, 'AAAG', 'A')),

    # Single letter insert, delete, duplication
    ('NM_000016.4:c.203delA',           ('chr1', 76198412, 'GA', 'G')),
    ('NM_000016.4:c.244dupT',           ('chr1', 76198564, 'C', 'CT')),
    ('NM_000016.4:c.387+1delG',         ('chr1', 76199309, 'TG', 'T')),
    ('NM_000016.4:c.475delT',           ('chr1', 76205669, 'AT', 'A')),
    ('NM_000016.4:c.1189dupT',          ('chr1', 76227049, 'C', 'CT')),
    ('NM_000016.4:c.1191delT',          ('chr1', 76227051, 'AT', 'A')),
    ('NM_000016.4:c.306_307insG',       ('chr1', 76199232, 'T', 'TG')),

    # Alignment tests for HGVS 3' and VCF left-alignment
    ('NM_000492.3:c.935_937delTCT',     ('chr7', 117180210, 'CCTT', 'C')),
    ('NM_000492.3:c.442delA',           ('chr7', 117171120, 'CA', 'C')),
    ('NM_000492.3:c.805_806delAT',      ('chr7', 117176660, 'AAT', 'A')),
    ('NM_000492.3:c.1155_1156dupTA',    ('chr7', 117182104, 'A', 'AAT')),
    ('NM_000492.3:c.3889dupT',          ('chr7', 117292905, 'A', 'AT')),

    # Transcript prefix (var_canonical=True)
    ('NM_007294.3:c.2207A>C',           ('chr17', 41245341, 'T', 'G')),

    # After stop codon
    ('NM_000492.3:c.*3A>C',             ('chr7', 117307165, 'A', 'C')),
]


@pytest.mark.parametrize("hgvs_name,expected", _name_to_variant)
def test_name_to_variant(biocommons_converter_grch37, hgvs_name, expected):
    result = biocommons_converter_grch37.hgvs_to_variant_coordinate(hgvs_name)
    assert result == expected
