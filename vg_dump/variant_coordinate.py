class VariantCoordinate(FormerTuple, pydantic.BaseModel):
    """ This stores coordinates, when you want to use it, be sure to call either:
            * as_external_explicit() - if you want to interface with eg VCF
            * as_internal_symbolic() - to store in internal Database
    """
    chrom: str
    position: int
    ref: str
    alt: str
    svlen: Optional[int] = None

    @field_validator('svlen')
    @classmethod
    def validate_svlen(cls, value):
        if value is not None and settings.VARIANT_SYMBOLIC_ALT_ENABLED is False:
            raise ValueError('Symbolic variants disabled via settings.')
        return value

    @property
    def as_tuple(self) -> tuple:
        return self.chrom, self.position, self.ref, self.alt, self.svlen

    @property
    def end(self) -> int:
        """
            This corresponds to VCF INFO["END"] which is defined in spec as:

            END position of the longest variant described in this record. The END of each allele is defined as:
                * Non-symbolic alleles: POS + length of REF allele − 1.
                * <INS> symbolic structural variant alleles: POS + length of REF allele − 1.
                * <DEL>, <DUP>, <INV>, and <CNV> symbolic structural variant alleles:, POS + SVLEN
        """
        if self.is_symbolic:
            # We don't support <INV> so don't need to worry about it
            return self.position + abs(self.svlen)
        return self.position + len(self.ref) - 1

    def __lt__(self, other):
        return self.as_tuple < other.as_tuple

    def __str__(self):
        return self.format()

    def __hash__(self):
        return self.as_tuple.__hash__()

    def format(self):
        if Sequence.allele_is_symbolic(self.alt):
            range_end = self.position + abs(self.svlen)
            return f"{self.chrom}:{self.position}-{range_end} {self.alt}"
        else:
            return f"{self.chrom}:{self.position} {self.ref}>{self.alt}"

    @staticmethod
    def from_variant_match(match, genome_build: Optional[GenomeBuild] = None):
        chrom = match.group(1)
        if genome_build:
            chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        position = int(match.group(2))
        ref = match.group(3)
        alt = match.group(4)
        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref, alt)
        return vc.as_internal_canonical_form(genome_build)

    @staticmethod
    def from_symbolic_match(match, genome_build):
        chrom, start, range_end, alt = match.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        alt = alt.upper()
        start = int(start)
        range_end = int(range_end)
        if range_end < start + 1:
            raise ValueError("End coordinate must be at least 1 base after start")

        if alt == "DEL":
            svlen = start - range_end
        else:
            svlen = range_end - start
        alt = "<" + alt + ">"  # captured what was inside of brackets ie <()>
        contig_sequence = genome_build.genome_fasta.fasta[chrom]
        # 0-based
        ref = contig_sequence[start-1:start].upper()
        vc = VariantCoordinate(chrom=chrom, position=start, ref=ref, alt=alt, svlen=svlen)
        return vc.as_internal_canonical_form(genome_build)

    @staticmethod
    def from_string(variant_string: str, genome_build):
        """ Pass in genome build to be able to set REF from symbolic (will be N otherwise) """
        if full_match := VARIANT_PATTERN.fullmatch(variant_string):
            return VariantCoordinate.from_variant_match(full_match, genome_build)
        elif full_match := VARIANT_SYMBOLIC_PATTERN.fullmatch(variant_string):
            return VariantCoordinate.from_symbolic_match(full_match, genome_build)
        regex_patterns = ", ".join((str(s) for s in (VARIANT_PATTERN, VARIANT_SYMBOLIC_PATTERN)))
        raise ValueError(f"{variant_string=} did not match against {regex_patterns=}")

    @staticmethod
    def from_explicit_no_svlen(chrom: str, position: int, ref: str, alt: str) -> 'VariantCoordinate':
        """ Initialise w/o providing an end """

        ref = ref.strip().upper()
        alt = alt.strip().upper()
        if Sequence.allele_is_symbolic(alt):
            raise ValueError("Must pass 'svlen' when using symbolic alt")
        return VariantCoordinate(chrom=chrom, position=position, ref=ref, alt=alt)

    @property
    def is_symbolic(self):
        return Sequence.allele_is_symbolic(self.alt)

    def calculated_reference(self, genome_build) -> str:
        contig_sequence = genome_build.genome_fasta.fasta[self.chrom]
        # reference sequence is 0-based
        ref_sequence = contig_sequence[self.position - 1:self.end]
        return ref_sequence.upper()

    def as_external_explicit(self, genome_build) -> 'VariantCoordinate':
        """ explicit ref/alt """
        if self.is_symbolic:
            if self.svlen is None:
                raise ValueError(f"{self} has 'svlen' = None")

            ref_sequence = self.calculated_reference(genome_build)
            if self.alt == VCFSymbolicAllele.DEL:
                ref = ref_sequence
                alt = ref_sequence[0]
            elif self.alt == VCFSymbolicAllele.DUP:
                ref = ref_sequence[0]
                alt = ref_sequence
            elif self.alt == VCFSymbolicAllele.INV:
                ref = ref_sequence
                alt = reverse_complement(ref_sequence)
            else:
                raise ValueError(f"Unknown symbolic alt of '{self.alt}'")

            vc = VariantCoordinate(chrom=self.chrom, position=self.position, ref=ref, alt=alt)
            # print(f"Symbolic = {self} -> {vc}")
            return vc

        if self.alt == Variant.REFERENCE_ALT:
            # Convert from our internal format (alt='=' for ref) to explicit
            alt = self.ref
        else:
            alt = self.alt
        return VariantCoordinate(chrom=self.chrom, position=self.position, ref=self.ref, alt=alt)

    @property
    def max_sequence_length(self) -> int:
        if self.is_symbolic:
            return abs(self.svlen)
        else:
            return max(len(self.ref), len(self.alt))

    def as_internal_symbolic(self, genome_build: GenomeBuild,
                             min_symbolic_alt_size=settings.VARIANT_SYMBOLIC_ALT_SIZE) -> 'VariantCoordinate':
        """ Internal format - alt can be <DEL> or <DUP>
            Uses our internal reference representation
        """
        if self.is_symbolic:
            return self

        ref = self.ref
        alt = self.alt  # default, only change if symbolic
        svlen = None
        if self.alt in (Variant.REFERENCE_ALT, self.ref):
            alt = Variant.REFERENCE_ALT
        else:
            # Could be symbolic
            ref_length = len(ref)
            alt_length = len(self.alt)
            diff = alt_length - ref_length
            if ref_length == 1:
                if alt_length >= min_symbolic_alt_size:
                    # Possible dup
                    # TODO: Can probably remove HGVS dependency from here - just look directly at sequence
                    from genes.hgvs import HGVSMatcher
                    matcher = HGVSMatcher(genome_build)
                    hgvs_variant = matcher.variant_coordinate_to_hgvs_variant(self)
                    if hgvs_variant.mutation_type == 'dup':
                        ref = self.ref[0]
                        alt = VCFSymbolicAllele.DUP
                        svlen = diff
            elif alt_length == 1 and self.alt == self.ref[0]:
                if ref_length >= min_symbolic_alt_size:
                    ref = self.ref[0]
                    alt = VCFSymbolicAllele.DEL
                    svlen = diff

            elif ref_length == alt_length:
                if ref_length > min_symbolic_alt_size:
                    if self.ref == reverse_complement(self.alt):
                        ref = self.ref[0]
                        alt = VCFSymbolicAllele.INV
                        svlen = len(self.ref) - 1  # explicit inv had same length ref/alt, now we have len(ref) == 1
        return VariantCoordinate(chrom=self.chrom, position=self.position, ref=ref, alt=alt, svlen=svlen)

    def as_internal_canonical_form(self, genome_build: GenomeBuild) -> 'VariantCoordinate':
        """ Make sure we only have 1 representation for a variant """

        # Easiest way is to just convert to symbolic then check svlen
        vc_symbolic = self.as_internal_symbolic(genome_build)
        if vc_symbolic.svlen and abs(vc_symbolic.svlen) >= settings.VARIANT_SYMBOLIC_ALT_SIZE:
            vc = vc_symbolic
        elif self.is_symbolic:
            vc = self.as_external_explicit(genome_build)
        else:
            vc = self

        if vc.alt in (Variant.REFERENCE_ALT, vc.ref):
            vc.alt = Variant.REFERENCE_ALT
        return vc

    def as_contig_accession(self, genome_build: GenomeBuild) -> 'VariantCoordinate':
        contig = genome_build.chrom_contig_mappings[self.chrom]
        return VariantCoordinate(chrom=contig.refseq_accession,
                                 position=self.position, ref=self.ref, alt=self.alt, svlen=self.svlen)

