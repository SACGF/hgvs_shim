# hgvs_shim

hgvs-shim is a small compatibility layer that simplifies migrating code from Counsyl HGVS python library
to Biocommons HGVS.

# Background

There are 2 Python HGVS libraries:

https://github.com/counsyl/hgvs - pyhvs
https://github.com/biocommons/hgvs Biocommons HGVS

Pyhgvs is simpler, but appears to have been abandoned (no comments on issues or commits for years)
Biocommons HGVS under active development, and has more features (such as projecting across alignment gaps)
inversions and uncertain coordinates.

# Motivation

In the VariantGrid project, we intially started with PyHGVS as it was simpler, so we found it easier to hook in
transcript provision, and use our own transcript versions (from GTF files). We eventually decided the complexity
of Biocommons HGVS was worth it, and wrote cdot for providing transcripts for both libraries.

We had made a number of fixes to PyHGVS, and workarounds for bug fixes, and the PyHGVS clasess were spread throuhgout
our code base. We needed to perform a lot of validation to ensure that Biocommons HGVS would meet all our requirements
during a long testing phase. The transition steps were:

* Abstract away PyHGVS classes and methods to our own classes
* Write implementations for both PyHGVS and biocommons
* Be able to change which is used from simple configuration hidden in a factory method
* Be able to eg calculate both, and throw an exception if they were different

HGVS shim is a cleaned up, library from this work.

