"""
Integration tests for PyHGVSConverter field-level correctness.
Requires pyhgvs extra and a local FASTA file.
"""
import pytest

pytestmark = [pytest.mark.integration, pytest.mark.pyhgvs]


def test_variant_coordinate_to_g(pyhgvs_converter_grch37):
    pytest.importorskip('pyhgvs')
    import pyhgvs
    from pysam import FastaFile

    fasta = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
    chrom = "NC_000007.13"
    position = 36561663
    ref = "A"
    alt = "T"

    hgvs_name_g = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt, fasta, None)
    hgvs_variant_g = pyhgvs_converter_grch37.variant_coordinate_to_g_hgvs(chrom, position, ref, alt)
    _assert_fields_equal(hgvs_name_g, hgvs_variant_g)


def test_variant_coordinate_to_c(pyhgvs_converter_grch37, pyhgvs_transcript_factory_grch37):
    pytest.importorskip('pyhgvs')
    import pyhgvs
    from pysam import FastaFile
    from hgvs_shim.hgvs_converter import TranscriptInfo

    fasta = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
    transcript_accession = "NM_001637.3"
    chrom = "NC_000007.13"
    position = 36561663
    ref = "A"
    alt = "T"

    transcript = pyhgvs_transcript_factory_grch37.get_transcript(transcript_accession, "GRCh37")
    hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt, fasta, transcript)

    transcript_info = TranscriptInfo(accession=transcript_accession, strand='+', is_coding=True)
    hgvs_variant = pyhgvs_converter_grch37.variant_coordinate_to_c_hgvs(
        chrom, position, ref, alt, transcript_info
    )
    _assert_fields_equal(hgvs_name, hgvs_variant)


def _assert_fields_equal(hgvs_name, other_hgvs_variant):
    assert hgvs_name.transcript == other_hgvs_variant.transcript
    assert hgvs_name.gene == other_hgvs_variant.gene
    assert hgvs_name.kind == other_hgvs_variant.kind
    assert hgvs_name.mutation_type == other_hgvs_variant.mutation_type
    assert hgvs_name.format() == other_hgvs_variant.format()
    assert hgvs_name.get_ref_alt() == other_hgvs_variant.get_ref_alt()
