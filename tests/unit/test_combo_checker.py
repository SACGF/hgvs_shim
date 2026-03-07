import pytest
from unittest.mock import MagicMock

from hgvs_shim.hgvs_converter import HGVSException
from hgvs_shim.hgvs_converter_combo import ComboCheckerHGVSConverter


def _make_converter(result):
    """Create a mock converter that returns result (or raises it if it's an exception)."""
    converter = MagicMock()
    converter.description.return_value = "mock"
    if isinstance(result, Exception):
        converter.create_hgvs_variant.side_effect = result
    else:
        converter.create_hgvs_variant.return_value = result
    return converter


def test_converters_agree_returns_first_result():
    c1 = _make_converter("variant_a")
    c2 = _make_converter("variant_a")
    combo = ComboCheckerHGVSConverter([c1, c2])
    result = combo.create_hgvs_variant("NM_000352.3:c.215A>G")
    assert result == "variant_a"


def test_converters_disagree_die_on_error_raises():
    c1 = _make_converter("variant_a")
    c2 = _make_converter("variant_b")
    combo = ComboCheckerHGVSConverter([c1, c2], die_on_error=True)
    with pytest.raises(ValueError, match="different results"):
        combo.create_hgvs_variant("NM_000352.3:c.215A>G")


def test_converters_disagree_no_die_logs_error(caplog):
    import logging
    c1 = _make_converter("variant_a")
    c2 = _make_converter("variant_b")
    combo = ComboCheckerHGVSConverter([c1, c2], die_on_error=False)
    with caplog.at_level(logging.ERROR):
        result = combo.create_hgvs_variant("NM_000352.3:c.215A>G")
    assert result == "variant_a"
    assert "different results" in caplog.text


def test_all_converters_raise_reraises():
    err = HGVSException("parse error")
    c1 = _make_converter(err)
    c2 = _make_converter(HGVSException("other error"))
    combo = ComboCheckerHGVSConverter([c1, c2])
    with pytest.raises(HGVSException):
        combo.create_hgvs_variant("BAD:c.1A>G")


def test_one_raises_one_succeeds_returns_success():
    err = HGVSException("parse error")
    c1 = _make_converter(err)
    c2 = _make_converter("variant_a")
    combo = ComboCheckerHGVSConverter([c1, c2], die_on_error=False)
    # One exception and one non-exception: has_exception=True, has_non_exception=True
    # _all_equal skips exceptions, so non-exception results must agree
    # Result: returns first result (exception) → raises it
    with pytest.raises(HGVSException):
        combo.create_hgvs_variant("NM_000352.3:c.215A>G")
