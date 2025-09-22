import pytest
from pysam import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import RESTPyHGVSTranscriptFactory

@pytest.fixture(scope="session")
def fasta_file37():
    return FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")

@pytest.fixture(scope="session")
def pyhgvs_transcript_factory():
    return RESTPyHGVSTranscriptFactory()
