import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / 'data'
TEST_TRANSCRIPTS_JSON_GZ = str(DATA_DIR / 'test_transcripts.json.gz')

# Standard GRCh37 FASTA path — required for biocommons sequence lookups
GRCH37_FASTA = "/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"


def _fasta_path_or_skip():
    """Return FASTA path if present, otherwise skip the test."""
    if not Path(GRCH37_FASTA).exists():
        pytest.skip(f"Local FASTA not found: {GRCH37_FASTA}")
    return GRCH37_FASTA


@pytest.fixture(scope="session")
def biocommons_converter_grch37():
    from cdot.hgvs.dataproviders import JSONDataProvider
    from cdot.hgvs.dataproviders.fasta_seqfetcher import FastaSeqFetcher
    from hgvs_shim.hgvs_converter_biocommons import BioCommonsHGVSConverter
    fasta_path = _fasta_path_or_skip()
    hdp = JSONDataProvider([TEST_TRANSCRIPTS_JSON_GZ], seqfetcher=FastaSeqFetcher(fasta_path))
    return BioCommonsHGVSConverter('GRCh37', hdp)


@pytest.fixture(scope="session")
def pyhgvs_transcript_factory_grch37():
    pytest.importorskip('pyhgvs', reason="pyhgvs not installed (pip install 'hgvs_shim[pyhgvs]')")
    from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory
    return JSONPyHGVSTranscriptFactory([TEST_TRANSCRIPTS_JSON_GZ])


@pytest.fixture(scope="session")
def pyhgvs_converter_grch37(pyhgvs_transcript_factory_grch37):
    pytest.importorskip('pysam', reason="pysam not installed (pip install 'hgvs_shim[pyhgvs]')")
    fasta_path = _fasta_path_or_skip()
    from pysam import FastaFile
    from hgvs_shim.hgvs_converter_pyhgvs import PyHGVSConverter
    fasta = FastaFile(fasta_path)
    return PyHGVSConverter(fasta, pyhgvs_transcript_factory_grch37.get_transcript_grch37)
