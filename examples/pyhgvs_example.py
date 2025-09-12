import pyhgvs
from pysam.libcfaidx import FastaFile
from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, RESTPyHGVSTranscriptFactory

def print_hgvs(hgvs_name):
    print(f"""
    {hgvs_name=}
    {hgvs_name.name=}
    {hgvs_name.prefix=}
    {hgvs_name.chrom=}
    {hgvs_name.transcript=}
    {hgvs_name.gene=}
    {hgvs_name.kind=}
    {hgvs_name.mutation_type=}
    {hgvs_name.start=}
    {hgvs_name.end=}
    {hgvs_name.ref_allele=}
    {hgvs_name.ref2_allele=}
    {hgvs_name.alt_allele=}   
    {hgvs_name.cdna_start=}
    {hgvs_name.cdna_end=}
    {hgvs_name.get_ref_alt()}
    {hgvs_name.get_vcf_coords(transcript)}
    {hgvs_name.format()}
    """)


genome = FastaFile("/data/annotation/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz")
factory = RESTPyHGVSTranscriptFactory()
# factory = JSONPyHGVSTranscriptFactory(["./cdot-0.2.14.refseq.grch37.json.gz"])  # Uses local JSON file
(chrom, start, ref, alt) = pyhgvs.parse_hgvs_name('NM_001637.3:c.1582G>A', genome, get_transcript=factory.get_transcript_grch37)

transcript = factory.get_transcript("NM_001637.3", "GRCh37")

hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, start, ref, alt, genome, transcript)
print_hgvs(hgvs_name)

