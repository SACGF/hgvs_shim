import os.path
import re
from importlib import metadata
from typing import Optional

from bioutils.sequences import reverse_complement
from hgvs.parser import Parser as _HGVSParser
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.interface import Interface
from hgvs.edit import NARefAlt
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSError, HGVSInvalidVariantError, HGVSInvalidIntervalError, \
    HGVSNormalizationError, HGVSParseError, HGVSUnsupportedOperationError, HGVSInternalError, HGVSUsageError, \
    HGVSVerifyFailedError
from hgvs.extras.babelfish import Babelfish
from hgvs.normalizer import Normalizer
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator
from hgvs.variantmapper import VariantMapper

from hgvs_shim.hgvs_converter import HGVSConverter, HGVSNomenclatureException, HGVSImplementationException, \
    HGVSException, HGVSConverterType, TranscriptInfo
from hgvs_shim.hgvs_variant import HGVSVariant


_hgvs_parser = _HGVSParser()


def _looks_like_transcript(accession: str) -> bool:
    return accession.startswith(('NM_', 'NR_', 'XM_', 'XR_', 'ENST'))


def _chrom_to_ncbi(chrom: str) -> str:
    """Convert UCSC-style 'chr1'/'chrM' to NCBI-style '1'/'MT' for babelfish input."""
    if chrom == 'chrM':
        return 'MT'
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom


def _chrom_to_ucsc(chrom: str) -> str:
    """Convert NCBI-style '1'/'MT' to UCSC-style 'chr1'/'chrM' for VCF output."""
    if chrom == 'MT':
        return 'chrM'
    if not chrom.startswith('chr'):
        return 'chr' + chrom
    return chrom


class BioCommonsHGVSVariant(HGVSVariant):
    def __init__(self, sequence_variant: SequenceVariant):
        self._sequence_variant = sequence_variant

    def _get_length(self) -> int:
        return self._sequence_variant.posedit.length_change()

    def _get_contig_accession(self) -> str:
        return self._sequence_variant.ac

    def _get_gene(self) -> str:
        return self._sequence_variant.gene

    def _set_gene(self, value):
        self._sequence_variant.gene = value

    def _get_transcript(self) -> str:
        return self._sequence_variant.ac

    def _set_transcript(self, value):
        self._sequence_variant.ac = value

    def _get_kind(self) -> str:
        return self._sequence_variant.type

    def _set_kind(self, value):
        self._sequence_variant.type = value

    def _get_ref_allele(self):
        return self._sequence_variant.posedit.edit.ref

    def _set_ref_allele(self, value):
        self._sequence_variant.posedit.edit.ref = value

    def _get_mutation_type(self):
        biocommons_type = self._sequence_variant.posedit.edit.type
        if biocommons_type == "sub":
            mutation_type = ">"
        else:
            mutation_type = biocommons_type
        return mutation_type

    def get_ref_alt(self):
        edit = self._sequence_variant.posedit.edit
        ref = edit.ref or ''
        alt = edit.alt or ''
        return ref, alt

    def get_cdna_coords(self) -> str:
        return str(self._sequence_variant.posedit.pos.start)

    def format(self, use_delins_for_inv: bool = True, max_ref_length=None):
        sv: SequenceVariant = self._sequence_variant
        if use_delins_for_inv:
            if sv.posedit.edit.type == "inv":
                ref = sv.posedit.edit.ref
                sv.posedit.edit = NARefAlt(ref=ref, alt=reverse_complement(ref))

        conf = {}
        if max_ref_length:
            conf["max_ref_length"] = max_ref_length
        return sv.format(conf)

    def get_gene_symbol_if_no_transcript(self) -> Optional[str]:
        # Biocommons SequenceVariant works like:
        # NM_001145661.2:c.1113dup          transcript=NM_001145661.2, gene=None
        # NM_001145661.2(GATA2):c.1113dup   transcript=NM_001145661.2, gene=GATA2
        # GATA2:c.1113dup                   transcript=GATA2
        gene_symbol = None
        gene = self._get_gene()
        transcript = self._get_transcript()
        if not (gene and transcript):
            # Will always be transcript
            if not _looks_like_transcript(transcript):
                gene_symbol = transcript
        return gene_symbol


class BioCommonsHGVSConverter(HGVSConverter):
    hgvs_span_trailing_int_length_pattern = re.compile(r"(.*(?:del|dup|inv))(\d+)$")

    def __init__(self, assembly_name: str, hdp: Interface):
        self.hdp = hdp
        self.babelfish = Babelfish(self.hdp, assembly_name)
        self.am = AssemblyMapper(self.hdp,
                                 assembly_name=assembly_name,
                                 alt_aln_method='splign',
                                 replace_reference=True)
        self.ev = ExtrinsicValidator(self.hdp)
        self.norm_5p = Normalizer(self.hdp, shuffle_direction=5)
        self.no_validate_mapper = VariantMapper(self.hdp, replace_reference=True, prevalidation_level="NONE")
        self.no_validate_normalizer = Normalizer(self.hdp, cross_boundaries=True, validate=False,
                                                 variantmapper=self.no_validate_mapper)

    @staticmethod
    def _parser_hgvs(hgvs_string: str) -> SequenceVariant:
        """ All calls to parsing go through here """

        HGVSConverter._hgvs_string_validation(hgvs_string)

        # Biocommons HGVS doesn't accept integers on the end of indels - ie NM_001354689.1(RAF1):c.1_2dup3
        # We want to strip these and raise an error if the span is wrong
        provided_span_length = None
        if m := BioCommonsHGVSConverter.hgvs_span_trailing_int_length_pattern.match(hgvs_string):
            hgvs_string, provided_span_length = m.groups()
            provided_span_length = int(provided_span_length)

        try:
            sequence_variant = _hgvs_parser.parse(hgvs_string)
        except HGVSError as hgvs_error:
            klass = BioCommonsHGVSConverter._get_exception_class(hgvs_error)
            raise klass(hgvs_error) from hgvs_error

        if provided_span_length is not None:
            if sequence_variant.posedit.edit.type == 'inv':
                # HGVS is 0 based
                coord_span = (sequence_variant.posedit.pos.end - sequence_variant.posedit.pos.start) + 1
            else:
                coord_span = abs(sequence_variant.posedit.length_change())
            if coord_span != provided_span_length:
                raise HGVSNomenclatureException(f"coordinate span ({coord_span}) not equal to provided ref length {provided_span_length}")
        return sequence_variant

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:

        try:
            sequence_variant = self._parser_hgvs(hgvs_string)
            return BioCommonsHGVSVariant(sequence_variant)
        except HGVSError as e:
            raise HGVSNomenclatureException from e

    def normalize(self, hgvs_variant: BioCommonsHGVSVariant) -> HGVSVariant:
        sv = hgvs_variant._sequence_variant
        sv_normalized = self.no_validate_normalizer.normalize(sv)
        return BioCommonsHGVSVariant(sv_normalized)

    def _variant_coordinate_to_sequence_variant(self, chrom: str, position: int, ref: str, alt: str) -> SequenceVariant:
        return self.babelfish.vcf_to_g_hgvs(_chrom_to_ncbi(chrom), position, ref, alt)

    def variant_coordinate_to_g_hgvs(self, chrom: str, position: int, ref: str, alt: str) -> HGVSVariant:
        var_g = self._variant_coordinate_to_sequence_variant(chrom, position, ref, alt)
        return BioCommonsHGVSVariant(var_g)

    def variant_coordinate_to_c_hgvs(self, chrom: str, position: int, ref: str, alt: str, transcript_info: TranscriptInfo) -> HGVSVariant:
        """ In VG we call non-coding "c.HGVS" as well - so have to handle that """
        try:
            var_g = self._variant_coordinate_to_sequence_variant(chrom, position, ref, alt)  # returns normalized (default HGVS 3')
            # Biocommons HGVS doesn't normalize introns as it works with transcript sequences so doesn't have introns
            # workaround is to normalize on genome sequence first, so if it can't norm it's correct
            if transcript_info.strand == '-':
                var_g = self.norm_5p.normalize(var_g)

            if transcript_info.is_coding:
                var_c = self.am.g_to_c(var_g, transcript_info.accession)
            else:
                var_c = self.am.g_to_n(var_g, transcript_info.accession)
        except HGVSError as e:  # Can be out of bounds etc
            klass = self._get_exception_class(e)
            raise klass(e) from e

        if gene_symbol := transcript_info.gene_symbol:
            var_c.gene = gene_symbol
        return BioCommonsHGVSVariant(var_c)

    def hgvs_to_variant_coordinate(self, hgvs_string: str) -> tuple[str, int, str, str]:
        try:
            var_g = self._hgvs_to_g_hgvs(hgvs_string)
            (chrom, position, ref, alt, _typ) = self.babelfish.hgvs_to_vcf(var_g)
            if alt == '.':
                alt = ref
            chrom = _chrom_to_ucsc(chrom)
        except HGVSError as hgvs_error:
            klass = self._get_exception_class(hgvs_error)
            raise klass(hgvs_error) from hgvs_error
        return chrom, position, ref, alt

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        sequence_variant = self._parser_hgvs(hgvs_string)
        sequence_variant.gene = None
        return sequence_variant.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        """ Only returns anything if c. HGVS """
        transcript_accession = ''
        if hgvs_string is not None:
            sequence_variant = self._parser_hgvs(hgvs_string)
            transcript_accession = self._get_transcript_accession_from_sequence_variant(sequence_variant)
        return transcript_accession

    @staticmethod
    def _get_transcript_accession_from_sequence_variant(sequence_variant: SequenceVariant) -> str:
        transcript_accession = ''
        if sequence_variant.type != 'g':
            if _looks_like_transcript(sequence_variant.ac):
                transcript_accession = sequence_variant.ac
        return transcript_accession

    @staticmethod
    def _strip_common_prefix(ref: str, alt: str) -> tuple[str, str]:
        if common_prefix := os.path.commonprefix((ref, alt)):
            i = len(common_prefix)
            ref = ref[i:]
            alt = alt[i:]
        return ref, alt

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.BIOCOMMONS_HGVS

    def get_version(self) -> str:
        return metadata.version('hgvs')

    @staticmethod
    def _m_to_g(var_m):
        # mito is basically the same as genomic except for the letter
        var_m.type = 'g'
        return var_m

    def _hgvs_to_g_hgvs(self, hgvs_string: str) -> SequenceVariant:
        CONVERT_TO_G = {
            'c': self.am.c_to_g,
            'n': self.am.n_to_g,
            'm': self._m_to_g,
        }

        var_x = self._parser_hgvs(hgvs_string)
        if converter := CONVERT_TO_G.get(var_x.type):
            var_g = converter(var_x)
        else:
            var_g = var_x
        return var_g

    @staticmethod
    def _get_exception_class(hgvs_error: HGVSError) -> type:
        """ Convert from HGVS to our generic errors """

        exception_mappings = {
            HGVSNomenclatureException: {
                HGVSInvalidIntervalError,
                HGVSInvalidVariantError,
                HGVSNormalizationError,
                HGVSParseError,
                HGVSUnsupportedOperationError
            },
            HGVSImplementationException: {
                HGVSDataNotAvailableError,
                HGVSInternalError,
                HGVSUsageError,
                HGVSVerifyFailedError,
            },
        }

        for our_ex, biocommons_hgvs_exceptions in exception_mappings.items():
            for hgvs_ex in biocommons_hgvs_exceptions:
                if isinstance(hgvs_error, hgvs_ex):
                    return our_ex
        return HGVSException  # General one...
