The purpose of "vg_dump" is to act as a code repository for Claude to examine in case I forgot to copy any code over
to the shim code. You do not want to implement all features in there into the shim - which is JUST to help migrate
from pyhgvs to biocommons.

The goal is to delete vg_Dump once done. 


---                                                                                                                          
  Analysis                                                                                                                     
                                                                                                                               
  The Core Approach is Right                                                                                                   
                                                                  
  The abstract HGVSConverter/HGVSVariant shim pattern is the correct architecture for this problem. The VG dump confirms it
  actually worked in production. Keep it.

  ---
  Active Bug: all_equal is Undefined

  hgvs_converter_combo.py:42 calls all_equal(results) but this function is never defined or imported anywhere in the shim. This
   is a NameError waiting to happen. In the VG code it was from library.utils import all_equal. Needs a simple local
  definition.

  ---
  The pyhgvs Test Goal is Harder Than It Looks

  The tests in tests/pyhgvs_tests/ cannot currently run — they use relative imports expecting to be inside the pyhgvs package:

  from .. import CDNACoord, HGVSName, parse_hgvs_name, format_hgvs_name
  from ..utils import read_transcripts
  from ..variants import normalize_variant

  There are three levels of things being tested:

  ┌───────┬──────────────────────────────────────┬─────────────────────────────────────────────────────────────────────────┐
  │ Level │                 What                 │                    Feasibility with biocommons shim                     │
  ├───────┼──────────────────────────────────────┼─────────────────────────────────────────────────────────────────────────┤
  │ 1     │ CDNACoord parsing,                   │ These are pyhgvs-internal types. biocommons has equivalent logic but no │
  │       │ genomic_to_cdna_coord                │  CDNACoord type to shim. Very hard.                                     │
  ├───────┼──────────────────────────────────────┼─────────────────────────────────────────────────────────────────────────┤
  │ 2     │ normalize_variant                    │ pyhgvs internal VCF normalizer. Biocommons does different normalization │
  │       │                                      │  via HGVS. Not directly equivalent.                                     │
  ├───────┼──────────────────────────────────────┼─────────────────────────────────────────────────────────────────────────┤
  │ 3     │ Round-trip name↔variant              │ This is exactly what the shim should test. The data is gold.            │
  │       │ (_name_variants data)                │                                                                         │
  └───────┴──────────────────────────────────────┴─────────────────────────────────────────────────────────────────────────┘

  The right plan is NOT to make pyhgvs's internal function signatures work via biocommons (deep compat shim) — that would mean
  reimplementing CDNACoord, transcript.cdna_to_genomic_coord, etc. in terms of biocommons internals. Fragile, huge effort.

  The right plan IS: port the test DATA, not the test infrastructure. The _name_variants list in test_hgvs_names.py contains
  ~50 real HGVS↔variant coordinate pairs with real transcripts that are the actual acceptance criteria for biocommons. These
  should drive new tests against BioCommonsHGVSConverter.

  ---
  What the Shim is Missing (from VG Dump)

  Performance — ParserSingleton: The biocommons Parser() constructor takes ~1.5 seconds. <--- Dave Note - this was slow
  as it compiled the class from parsley upon startup. I made a fix to biocommons to compile this class so that it is
  almost instant - this is not in the package but is merged into biocommons main/master so installing from github performance is ok
    they will release a new version soon.


Current shim calls parser.parse_hgvs_variant() using the module-level from hgvs import parser — this is fine, but if someone creates multiple
  BioCommonsHGVSConverter instances they don't share state. The VG version used a class-level singleton.

  looks_like_transcript: Used internally by BioCommonsHGVSConverter in get_gene_symbol_if_no_transcript — it's already calling
  looks_like_transcript but this is imported from pyhgvs (from pyhgvs import get_refseq_type). This is an undocumented pyhgvs
  internal dependency. Should be a standalone utility.

  CHGVS utility class: The VG CHGVS class is a genuinely useful standalone utility for parsing and manipulating c.HGVS strings
  without needing a full converter. It handles: extracting transcript, gene, c_dot; diff()ing two c.HGVS strings; normalizing
  format. It doesn't require a FASTA or data provider. Would help migration codebases that do a lot of string manipulation.

  clean_hgvs: The VG HGVSMatcher.clean_hgvs() handles real-world messy HGVS (missing c., doubled colons, wrong case,
  gene/transcript swapped). This is very useful for any production system accepting user input.

  HgvsMatchRefAllele: Tells you if the ref allele you provided matches the genome. Important for user-facing warnings during
  migration (when pyhgvs was more lenient than biocommons about ref mismatches).

  ---
  The Test Infrastructure Gap

  To test the biocommons converter properly (especially the _name_variants data), you need a transcript data provider (hdp).
  The current conftest.py only provides a FASTA file for pyhgvs. For biocommons you need something implementing
  hgvs.dataproviders.interface.Interface.

  cdot provides this via cdot.hgvs.dataproviders.JSONDataProvider pointing at a local JSON file. The test transcripts
  (NM_000352.3, NM_007294.3, NM_000016.4, etc.) would need to be in that JSON file.

    You can also use CDOT rest implementation to retrieve that off the internet.
We could download just what we need to build the tests then store that json locally.

  ---
  Plan

  Phase 1: Fix and Stabilize

  1. Fix the all_equal bug

  # src/hgvs_shim/hgvs_converter_combo.py

  def _all_equal(iterable):
      """Return True if all elements are equal (ignores exceptions)."""
      items = [x for x in iterable if not isinstance(x, Exception)]
      return len(set(items)) <= 1

  (Or define it as a module-level helper. The point is it's missing entirely.)

  2. Fix the hidden pyhgvs dependency in hgvs_converter_biocommons.py

  Currently from pyhgvs import get_refseq_type is used to detect RNA transcripts. This creates a pyhgvs dependency that's not
  stated as required for the biocommons converter path:

  # current - bad: biocommons converter importing from pyhgvs
  from pyhgvs import get_refseq_type

  # fix: inline the simple regex
  def _get_refseq_type(transcript_accession: str) -> str:
      """Determine refseq type from accession prefix."""
      if transcript_accession.startswith('NR_'):
          return 'RNA'
      if transcript_accession.startswith('NM_'):
          return 'mRNA'
      return 'unknown'

  3. Clarify the BioCommonsHGVSConverter constructor

  Currently the constructor takes assembly_name: str and hdp: Interface with no guidance on how to create them. Add a factory
  classmethod:

  @classmethod
  def from_local_file(cls, assembly_name: str, transcript_json_path: str) -> 'BioCommonsHGVSConverter':
      """Convenience constructor using cdot JSONDataProvider."""
      from cdot.hgvs.dataproviders import JSONDataProvider
      hdp = JSONDataProvider([transcript_json_path])
      return cls(assembly_name, hdp)

  ---
  Phase 2: Restructure the Tests

  2a. Fix the pyhgvs test imports

  The tests/pyhgvs_tests/ files currently import from .. expecting pyhgvs's package. Two options:

  - Option A: Delete tests/pyhgvs_tests/ and port the test DATA to new files
  - Option B: Keep them as reference only (add a note they're not executable), port DATA

  I recommend Option A. The value is in the test data, not the test infrastructure.

  2b. New test structure

  tests/
    conftest.py                        # shared fixtures
    unit/
      test_hgvs_string_validation.py   # tests HGVSConverter._hgvs_string_validation (no deps, runs anywhere)
      test_combo_checker.py            # tests ComboCheckerHGVSConverter logic with mocks
    integration/
      conftest.py                      # fixtures requiring cdot JSON / FASTA files
      test_biocommons_name_to_variant.py  # the ported _name_variants data
      test_biocommons_variant_to_name.py
      test_pyhgvs_name_to_variant.py   # same data, pyhgvs backend
      test_combo_agreement.py          # run both and compare

  2c. The ported test looks like this

  # tests/integration/test_biocommons_name_to_variant.py
  import pytest
  from hgvs_shim.hgvs_converter_biocommons import BioCommonsHGVSConverter

  # Ported directly from pyhgvs test_hgvs_names.py _name_variants
  # format: (hgvs_name, (chrom, position, ref, alt))
  _name_variants = [
      ('NM_000352.3:c.215A>G',  ('chr11', 17496508, 'T', 'C')),
      ('NM_000352.3:c.72C>A',   ('chr11', 17498252, 'G', 'T')),
      ('NM_000352.3:c.3885C>G', ('chr11', 17418843, 'G', 'C')),
      # intronic
      ('NM_000352.3:c.1630+1G>A', ('chr11', 17464266, 'C', 'T')),
      # CFTR deletions (alignment tests)
      ('NM_000492.3:c.935_937delTCT', ('chr7', 117180210, 'CCTT', 'C')),
      ('NM_000492.3:c.1155_1156dupTA', ('chr7', 117182104, 'A', 'AAT')),
      # ... all 50+ cases
  ]

  @pytest.mark.parametrize("hgvs_name,expected", _name_variants)
  def test_name_to_variant(biocommons_converter_grch37, hgvs_name, expected):
      chrom, pos, ref, alt = biocommons_converter_grch37.hgvs_to_variant_coordinate(hgvs_name)
      assert (chrom, pos, ref, alt) == expected

  The biocommons_converter_grch37 fixture needs a cdot data provider with those transcripts. This is the big infrastructure
  question: can we bundle a small cdot JSON file for testing? The transcripts involved are: NM_000352.3, NM_007294.3,
  NM_000016.4, NM_000492.3, NM_000038.5, NM_000018.3, NM_000019.3, NM_000022.2, NM_000023.2, NM_000030.2, ENST00000357654. A
  cdot JSON covering just these would be tiny and committable.

  ---
  Phase 3: Add Useful Utilities

  3a. CHGVS class — a standalone parser for c.HGVS strings

  The VG version has Django/DB dependencies we can strip. A standalone version:

  # src/hgvs_shim/chgvs.py
  import re

  class CHGVS:
      """Parse and manipulate a c.HGVS string without needing a converter.
   
      Examples:
          >>> c = CHGVS("NM_001637.3(ABCA4):c.5882G>A")
          >>> c.transcript
          'NM_001637.3'
          >>> c.gene
          'ABCA4'
          >>> c.c_dot
          'c.5882G>A'
      """
      _REGEX = re.compile(r'^(.*?)(?:\((.+?)\))?:([a-z]\..*)$')

      def __init__(self, full_hgvs: str):
          self.full_hgvs = full_hgvs or ""
          self.transcript = ""
          self.gene = None
          self.c_dot = None
          if m := self._REGEX.match(full_hgvs):
              self.transcript = m.group(1)
              self.gene = m.group(2)
              self.c_dot = m.group(3)

      def with_gene(self, gene: str) -> 'CHGVS':
          return CHGVS(f"{self.transcript}({gene}):{self.c_dot}")

      def without_gene(self) -> 'CHGVS':
          return CHGVS(f"{self.transcript}:{self.c_dot}")

      def __str__(self):
          return self.full_hgvs

  3b. clean_hgvs — accept real-world malformed input

  Extracted from VG's HGVSMatcher.clean_hgvs(), stripped of Django dependencies:

  # src/hgvs_shim/hgvs_cleaner.py

  def clean_hgvs(hgvs_string: str) -> tuple[str, list[str]]:
      """
      Clean common formatting errors from a user-provided HGVS string.
   
      Returns (cleaned_string, list_of_messages_about_what_was_fixed)
   
      Handles:
      - Whitespace removal
      - Double colons "::" → ":"
      - Uppercase mutation types (DEL → del, DUP → dup)
      - Missing c. prefix on transcript HGVS
      - Missing g. prefix on contig HGVS
      - Gene/transcript swap: "BRCA1(NM_007294.3):c..." → "NM_007294.3(BRCA1):c..."
      - Lowercase transcript accessions
      - Unbalanced brackets
      """
      ...

  This is one of the most valuable things the VG codebase developed — a large, real-world test suite of malformed HGVS strings.
   This should be extracted and have its own tests.

  ---
  Phase 4: Public API Cleanup

  4a. The format parameter naming

  In the VG version, format(use_compat=False) converts inv to delins. In the current shim it's format(allow_inv=False) but the
  logic is inverted (allow_inv=False strips inv → delins, allow_inv=True keeps inv). The parameter name is confusing. Consider:

  def format(self, use_delins_for_inv: bool = True, max_ref_length: int | None = None) -> str:

  4b. variant_coordinate_to_c_hgvs signature mismatch

  PyHGVS version takes transcript_accession: str, biocommons takes transcript_version (an object with .accession, .strand,
  .is_coding, .gene_symbol). This is the biggest incompatibility between the two backends - the biocommons converter needs
  richer transcript info to do the g→c mapping correctly (especially strand for normalization).

  The current shim has this mismatch and it's a real problem. Options:
  - Define a TranscriptInfo dataclass in the shim that both implementations accept
  - Document that the biocommons path requires a richer object

  A shim TranscriptInfo:

  # src/hgvs_shim/hgvs_converter.py
  from dataclasses import dataclass

  @dataclass
  class TranscriptInfo:
      accession: str          # e.g. "NM_001637.3"
      strand: str             # '+' or '-'
      is_coding: bool         # True for NM_, False for NR_
      gene_symbol: str | None = None

  Then both implementations accept TranscriptInfo and the pyhgvs one just uses accession.

  4c. Export the public API clearly

  # src/hgvs_shim/__init__.py
  from hgvs_shim.hgvs_converter import (
      HGVSConverter,
      HGVSConverterType,
      HGVSException,
      HGVSNomenclatureException,
      HGVSImplementationException,
      TranscriptInfo,
  )
  from hgvs_shim.hgvs_variant import HGVSVariant
  from hgvs_shim.hgvs_converter_biocommons import BioCommonsHGVSConverter
  from hgvs_shim.hgvs_converter_pyhgvs import PyHGVSConverter
  from hgvs_shim.hgvs_converter_combo import ComboCheckerHGVSConverter

  __all__ = [...]

  ---
  Clarifying Questions for You

  1. Test infrastructure: Do you have access to a cdot JSON file for the transcripts in the pyhgvs test suite (NM_000352.3,
  NM_007294.3, etc.), or would generating one be part of the work? The integration tests can't work without it.

Use cdot rest implementation - or just call eg "https://cdot.cc/transcript/NM_001001890.3" - subbing in the transcript version
store the local JSON

  2. pyhgvs dependency: Should PyHGVSConverter remain in the main package (requiring pyhgvs + pysam as dependencies), or should
   it move to an optional extra? For someone who has already migrated, they don't want pyhgvs as a dependency anymore.

Yeah, make them optional. Then the two installation targets are like full compatability layer then just the interfaces and biocommons hgvs once fully migrated

  3. CHGVS class scope: Do you want the string-parsing utility class, or is the converter API sufficient? The VG version uses
  it extensively for string manipulation outside the conversion context (comparing two c.HGVS values, sorting, etc).

    There is that class, and the VariantCoordinate class - if you think they are worth having then use them

    Don't only include them as a help - this is supposed to be a minimalist shim to help migration only. 

Just because it is in the VG dir doesn't mean you should use it

9. clean_hgvs scope: Is cleaning/normalizing dirty input strings in scope for this library, or should it stay out-of-scope
  (the library expects valid HGVS input)?

The cleaning is just legacy VG code I copied over - this should JUST be a shim libray - the only goal is to help migrate from
pyhgtvs to biocommons

  5. The HgvsMatchRefAllele question: The VG version tracks whether the reference allele you provided matches the genome
  sequence. This is important for clinical data import (users submit HGVS with wrong ref alleles). Is this something you want
  in the shim, or is it too application-specific?

This is just copying the VG code over - do not implement all the features in here, this was just in case I forgot to move something over
into the shim code.

Eventually I will delete the vg_dump