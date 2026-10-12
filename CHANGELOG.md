# Changelog

## 0.1.0 (2026-10-12)

First public release.

- VNIR (S-detector) summary parameters of Viviano-Beck et al. (2014), Table 2,
  computable from map-projected VNIR tiles as well as full-range cubes.
- Every non-retired IR (L-detector) Table 2 parameter between 1.0 and 3.9 um.
- Nine visible/near-infrared slope and ratio parameters adapted from the
  Curiosity-team conventions.
- The 16 standard browse-product recipes (Table 4) with percentile stretching.
- Readers for ENVI (TER3, MTR3, VRDR) and PDS3 (TRR3, WA CDR, DDR) products,
  an ENVI writer, and the `pycrism` command-line tool.
- Tests on synthetic cubes; continuous integration on Python 3.10 and 3.12.
