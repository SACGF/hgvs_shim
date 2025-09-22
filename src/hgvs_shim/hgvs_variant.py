import abc
from typing import Optional

class HGVSVariant(abc.ABC):
    """ This class wraps pyhgvs HGVSName and BioCommons SequenceVariant functionality,
        to allow library independent code """

    @property
    def contig_accession(self) -> str:
        _genomic_kinds = ('g', 'm')
        if self.kind not in _genomic_kinds:
            raise ValueError(f"'{self}' can only request contig for genomic kinds '{','.join(_genomic_kinds)}'")
        return self._get_contig_accession()

    @abc.abstractmethod
    def _get_contig_accession(self) -> str:
        pass

    @property
    def gene(self) -> str:
        return self._get_gene()

    @gene.setter
    def gene(self, value):
        self._set_gene(value)

    @abc.abstractmethod
    def _get_gene(self):
        pass

    @abc.abstractmethod
    def _set_gene(self, value):
        pass

    @property
    def length(self) -> int:
        return self._get_length()

    @abc.abstractmethod
    def _get_length(self) -> int:
        pass

    @property
    def transcript(self) -> str:
        return self._get_transcript()

    @transcript.setter
    def transcript(self, value):
        self._set_transcript(value)

    @abc.abstractmethod
    def _get_transcript(self):
        pass

    @abc.abstractmethod
    def _set_transcript(self, value):
        pass

    @property
    def kind(self) -> str:
        return self._get_kind()

    @kind.setter
    def kind(self, value):
        self._set_kind(value)

    @property
    def ref_allele(self) -> str:
        return self._get_ref_allele()

    @ref_allele.setter
    def ref_allele(self, value):
        self._set_ref_allele(value)

    @abc.abstractmethod
    def _get_kind(self):
        pass

    @abc.abstractmethod
    def _set_kind(self, value):
        pass

    @property
    def mutation_type(self) -> str:
        return self._get_mutation_type()

    @abc.abstractmethod
    def _get_mutation_type(self):
        pass

    @abc.abstractmethod
    def get_ref_alt(self):
        pass

    @abc.abstractmethod
    def get_cdna_coords(self) -> str:
        pass

    @abc.abstractmethod
    def _get_ref_allele(self):
        pass

    @abc.abstractmethod
    def _set_ref_allele(self, value):
        pass

    @abc.abstractmethod
    def format(self, allow_inv=False, max_ref_length=None):
        pass

    @abc.abstractmethod
    def get_gene_symbol_if_no_transcript(self) -> Optional[str]:
        pass

    def __repr__(self):
        return self.format()

    def __eq__(self, other):
        if isinstance(other, HGVSVariant):
            return self.format() == other.format()
        return False
