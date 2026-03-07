import pytest

from hgvs_shim.hgvs_converter import HGVSConverter, HGVSNomenclatureException


@pytest.mark.parametrize("hgvs_string", [
    "NM_000352.3:c.215A>G",
    "NM_007294.3:c.2207A>C",
    "NC_000007.13:g.117180210A>G",
])
def test_valid_hgvs_passes(hgvs_string):
    HGVSConverter._hgvs_string_validation(hgvs_string)


@pytest.mark.parametrize("hgvs_string,expected_msg", [
    ("NM_000352.3:c.215ins", "Insertions require inserted sequence"),
    ("NM_000352.3:c.215ins5", "Insertions require inserted sequence, not an integer length"),
])
def test_ins_without_sequence_raises(hgvs_string, expected_msg):
    with pytest.raises(HGVSNomenclatureException, match=expected_msg):
        HGVSConverter._hgvs_string_validation(hgvs_string)


@pytest.mark.parametrize("hgvs_string", [
    "NM_000352.3c.215A>G",  # no colon
    "c.215A>G",             # no colon + missing ref
    "g.117180210A>G",       # no colon + missing ref
    ".215A>G",              # no colon + missing ref
])
def test_missing_colon_raises(hgvs_string):
    with pytest.raises(HGVSNomenclatureException):
        HGVSConverter._hgvs_string_validation(hgvs_string)


@pytest.mark.parametrize("hgvs_string", [
    ":NM_000352.3:c.215A>G",   # starts with ':'
])
def test_missing_reference_raises(hgvs_string):
    with pytest.raises(HGVSNomenclatureException, match="Missing reference sequence"):
        HGVSConverter._hgvs_string_validation(hgvs_string)
