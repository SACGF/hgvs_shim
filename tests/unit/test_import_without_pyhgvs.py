"""pyhgvs is an optional extra - importing hgvs_shim without it must work."""
import os
import subprocess
import sys
import textwrap

import pytest

_SCRIPT = textwrap.dedent("""
    import sys


    class BlockOptionalExtras:
        def find_spec(self, name, path=None, target=None):
            if name.split(".")[0] in ("pyhgvs", "pysam"):
                raise ImportError(f"No module named '{name}'")
            return None


    sys.meta_path.insert(0, BlockOptionalExtras())

    import hgvs_shim

    assert hgvs_shim.BioCommonsHGVSConverter is not None

    try:
        hgvs_shim.PyHGVSConverter(fasta_file=None)
    except ImportError as e:
        assert "hgvs_shim[pyhgvs]" in str(e), e
    else:
        raise AssertionError("PyHGVSConverter should raise ImportError without the pyhgvs extra")
""")


def test_import_and_converter_error_without_pyhgvs():
    """Run in a subprocess so the block is in place before hgvs_shim is first imported."""
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p))
    result = subprocess.run([sys.executable, "-c", _SCRIPT], capture_output=True, text=True, env=env)
    if result.returncode:
        pytest.fail(f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")
