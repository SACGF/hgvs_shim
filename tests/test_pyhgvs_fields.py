import pyhgvs

from hgvs_shim.hgvs_converter_pyhgvs import PyHGVSConverter

def test_variant_coordinate_to_g(fasta_file37):
    chrom = "NC_000007.13"
    position = 36561663
    ref = "A"
    alt = "T"

    hgvs_name_g = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt, fasta_file37, None)
    hgvs_converter = PyHGVSConverter(fasta_file37)
    hgvs_variant_g = hgvs_converter.variant_coordinate_to_g_hgvs(chrom, position, ref, alt)
    _assert_fields_equal(hgvs_name_g, hgvs_variant_g)


def test_variant_coordinate_to_c(fasta_file37, pyhgvs_transcript_factory):
    transcript_accession = "NM_001637.3"
    chrom = "NC_000007.13"
    position = 36561663
    ref = "A"
    alt = "T"

    transcript = pyhgvs_transcript_factory.get_transcript(transcript_accession, "GRCh37")
    hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, position, ref, alt, fasta_file37, transcript)

    hgvs_converter = PyHGVSConverter(fasta_file37, pyhgvs_transcript_factory.get_transcript_grch37)
    hgvs_variant = hgvs_converter.variant_coordinate_to_c_hgvs(chrom, position, ref, alt, transcript_accession)
    _assert_fields_equal(hgvs_name, hgvs_variant, transcript)


def _assert_fields_equal(hgvs_name, other_hgvs_variant, transcript: str | None = None):
    # assert hgvs_name.name == other_hgvs_variant.name
    # assert hgvs_name.prefix == other_hgvs_variant.prefix
    # assert hgvs_name.chrom == other_hgvs_variant.chrom
    assert hgvs_name.transcript == other_hgvs_variant.transcript
    assert hgvs_name.gene == other_hgvs_variant.gene
    assert hgvs_name.kind == other_hgvs_variant.kind
    assert hgvs_name.mutation_type == other_hgvs_variant.mutation_type
    # assert hgvs_name.start == other_hgvs_variant.start
    # assert hgvs_name.end == other_hgvs_variant.end
    # assert hgvs_name.ref_allele == other_hgvs_variant.ref_allele
    # assert hgvs_name.ref2_allele == other_hgvs_variant.ref2_allele
    # assert hgvs_name.alt_allele == other_hgvs_variant.alt_allele
    # assert hgvs_name.cdna_start == other_hgvs_variant.cdna_start
    # assert hgvs_name.cdna_end == other_hgvs_variant.cdna_end

    # assert str(hgvs_name) == str(other_hgvs_variant)
    assert hgvs_name.format() == other_hgvs_variant.format()
    assert hgvs_name.get_ref_alt() == other_hgvs_variant.get_ref_alt()
    if transcript:
        # assert hgvs_name.get_vcf_coords(transcript) == other_hgvs_variant.get_vcf_coords(transcript)
        pass


