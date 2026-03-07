import abc
import re
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from hgvs_shim.hgvs_variant import HGVSVariant


class HGVSException(Exception):
    """ A wrapper for pyhgvs and Biocommons HGVS Exceptions to allow library independent code """


class HGVSNomenclatureException(HGVSException):
    """ HGVSException subclass for when problem is with HGVS string (users can fix) """


class HGVSImplementationException(HGVSException):
    """ HGVSException subclass for when problem is with the library (users can NOT fix) """


class HGVSConverterType(Enum):
    PYHGVS = 1
    BIOCOMMONS_HGVS = 2
    COMBO = 3


@dataclass
class TranscriptInfo:
    accession: str
    strand: str              # '+' or '-'
    is_coding: bool          # True for NM_/XM_, False for NR_/XR_
    gene_symbol: Optional[str] = None


class HGVSConverter(abc.ABC):
    """ This is the base object for PyHGVS and BioCommons HGVS
        implementations """

    @staticmethod
    def _hgvs_string_validation(hgvs_string: str):
        """ raise exceptions on any errors """

        if "ins" in hgvs_string:
            if re.match(r".*ins\d+$", hgvs_string):
                raise HGVSNomenclatureException("Insertions require inserted sequence, not an integer length")
            if re.match(".*ins$", hgvs_string):
                raise HGVSNomenclatureException("Insertions require inserted sequence")
        if ":" not in hgvs_string:
            raise HGVSNomenclatureException("No colon (':') provided")
        for char in [":", "c", "g", "."]:
            if hgvs_string.startswith(char):
                raise HGVSNomenclatureException("Missing reference sequence")

    @abc.abstractmethod
    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        pass

    @abc.abstractmethod
    def normalize(self, hgvs_variant: HGVSVariant) -> HGVSVariant:
        pass

    @abc.abstractmethod
    def variant_coordinate_to_g_hgvs(self, chrom: str, position: int, ref: str, alt: str) -> HGVSVariant:
        pass

    @abc.abstractmethod
    def variant_coordinate_to_c_hgvs(self, chrom: str, position: int, ref: str, alt: str, transcript_info: 'TranscriptInfo') -> HGVSVariant:
        pass

    @abc.abstractmethod
    def hgvs_to_variant_coordinate(self, hgvs_string: str) -> tuple[str, int, str, str]:
        pass

    @abc.abstractmethod
    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        pass

    @abc.abstractmethod
    def get_transcript_accession(self, hgvs_string: str) -> str:
        pass

    @abc.abstractmethod
    def get_hgvs_converter_type(self) -> HGVSConverterType:
        pass

    @abc.abstractmethod
    def get_version(self) -> str:
        pass

    def description(self) -> str:
        hgvs_converter_type = self.get_hgvs_converter_type()
        version = self.get_version()
        desc = f"{hgvs_converter_type.name} {version}"
        return desc
