# pycrism

[![tests](https://github.com/euchas/pycrism/actions/workflows/tests.yml/badge.svg)](https://github.com/euchas/pycrism/actions/workflows/tests.yml)

CRISM summary parameters, cube readers and browse products in pure NumPy.
No ENVI, IDL or CRISM Analysis Toolkit (CAT) needed.

CRISM is the imaging spectrometer on Mars Reconnaissance Orbiter (Murchie et
al., 2007). Its standard "summary parameters" (band depths, slopes, ratios and
indices tuned to specific minerals) are defined in Viviano-Beck et al. (2014),
Table 2. This package implements them as plain functions of a reflectance
array and a wavelength axis, so they can be run on any CRISM product, on a
laboratory spectrum, or on a single pixel.

What it does:

* **69 parameters**: the 14 VNIR (S-detector) parameters that fit inside the
  map-projected VNIR tiles (371-1017 nm), every non-retired IR (L-detector)
  parameter of Table 2 between 1.0 and 3.9 um (46), and 9 visible/near-infrared
  slope and ratio parameters adapted from the Curiosity rover conventions.
* **16 browse recipes** (TRU, FEM, FM2, MAF, PHY, CAR, ...) as percentile-stretched
  RGB images, after Table 4 of the same paper.
* **Readers** for ENVI products (TER3, MTR3, map-projected VNIR) and PDS3 products
  (TRR3 with the WA wavelength CDR, DDR geometry), an ENVI writer, and a
  command-line tool.
* **Bounded memory**: parameters work on a window of a few bands at a time and the
  whole-spectrum products run in pixel blocks, so a 4 GB map tile is never
  converted to float64 as a whole.

## Install

```
pip install git+https://github.com/euchas/pycrism
pip install "pycrism[plot] @ git+https://github.com/euchas/pycrism"   # PNG browse output
```

Requires Python 3.10 or newer and NumPy. Matplotlib is optional.

## Quickstart, command line

```
pycrism list                                   # parameters, band positions, recipes
pycrism params tile.img                        # summary statistics of every parameter the axis covers
pycrism params tile.img -p BD530_2 BD860_2 RPEAK1 -o tile_params.img   # ENVI BSQ with band names
pycrism browse tile.img FM2 -o tile_FM2.png    # browse image (PNG needs matplotlib; .img writes ENVI)
pycrism params frt0000xxxx_07_if166l_trr3.img --wa cdr/wa0000000000_01_wa0000l_cdr.img
```

## Quickstart, Python

```python
import pycrism

cube, wl, meta = pycrism.read_cube("tile.img")          # (lines, samples, bands) float32, nm
maps = pycrism.compute(cube, wl, ["BD530_2", "SH600_2", "BDI1000VIS"])
bd = pycrism.band_depth(cube, wl, 2250.0, 2290.0, 2350.0, kernels=(5, 5, 5))  # your own band
rgb = pycrism.browse_rgb(cube, wl, "FM2")              # (lines, samples, 3) in 0..1
stack, names = pycrism.stack_maps(maps)
pycrism.write_envi("tile_params.img", stack, band_names=names)
```

Spectra are bands-last: one spectrum `(bands,)`, a table `(n, bands)` or a
cube `(lines, samples, bands)`; a `(bands, lines, samples)` cube is accepted.
`wl` is the band-centre axis in nm (um is converted). Results keep the leading
shape and are NaN wherever an input is NaN or a needed wavelength is missing.

## Inputs and outputs

| Product | File pair | Reader | Notes |
|---|---|---|---|
| TER3, MTR3 (targeted and map-projected, S+L joined) | `.img` + `.hdr` | `read_envi` | wavelengths and bad-band list from the header |
| Map-projected VNIR tiles (VRDR) | `.img` + `.hdr` | `read_envi` | VNIR only; the two dropped-band gaps are respected |
| TRR3 (raw I/F, per detector) | `.img` + `.lbl` + WA CDR | `read_pds_cube`, `read_cube(..., wa=)` | descending L-detector axis is flipped; per-column wavelengths averaged |
| DDR (geometry) | `.img` + `.lbl` | `read_ddr` | named planes: lat, lon, incidence, emission, ... |
| Parameter maps, browse images | `.img` + `.hdr` | `write_envi` | float32 BSQ, NaN written as 65535 |

The null value 65535 and non-finite samples become NaN. Nothing here applies
photometric or atmospheric corrections: feed corrected I/F (TER/MTR products
or your own correction of TRR3).

## Parameters

The full list with band positions, kernel widths and what each definition
was checked against is in [docs/PARAMETERS.md](docs/PARAMETERS.md), and
`pycrism list` prints it. In short:

| Group | Parameters |
|---|---|
| VNIR, Table 2 | `R440` `R530` `R600` `R770` `RBR` `IRR1` `BD530_2` `SH600_2` `SH770` `BD640_2` `BD860_2` `BD920_2` `RPEAK1` `BDI1000VIS` |
| VNIR slopes (Curiosity conventions) | `VIS_SLOPE` `S600_840` `S750_840` `RED_BLUE` `RATIO_600_840` `NIR_RATIO` `BD535` `BD867W` `BD867N` |
| IR mafic and albedo | `BDI1000IR` `R1330` `BD1300` `OLINDEX3` `LCPINDEX2` `HCPINDEX2` `VAR` `ISLOPE1` `BDI2000` `R1080` `R1506` `R2529` `R3920` |
| IR hydration, sulfates, Al-OH, Fe/Mg-OH | `BD1400` `BD1750_2` `BD1900_2` `BD1900r2` `BD2100_2` `BD2165` `BD2190` `MIN2200` `BD2210_2` `D2200` `BD2230` `BD2250` `MIN2250` `BD2265` `BD2290` `D2300` `BD2355` `SINDEX2` `BD3000` |
| IR carbonates | `MIN2295_2480` `MIN2345_2537` `BD2500_2` `BD3400_2` `CINDEX2` |
| IR ices and atmosphere | `BD1435` `BD1500_2` `ICER1_2` `ICER2_2` `BD3100` `BD3200` `BD2600` `IRR2` `IRR3` |

Verification status:

* All VNIR Table 2 parameters except `RPEAK1` and `BDI1000VIS`, and all 46 IR
  parameters, use the wavelengths, continuum weights and median-kernel widths
  as printed in Table 2. Five IR parameters (`BDI1000IR`, `BDI2000`, `VAR`,
  `BD1500_2`, `BD1900_2`) use the CAT v7.4 hyperspectral form where the
  printed table describes the multispectral one; the docstrings say which.
  The IR functions were first written for the author's Jezero rim analysis
  ([jezero-rim-mapping](https://github.com/euchas/jezero-rim-mapping)) and
  were cross-checked there against the CRISM Data Product SIS, the CAT IDL
  code and two open ports (PyHAT, HyPyRameter).
* `RPEAK1` and `BDI1000VIS` follow the published method (fifth-order
  polynomial peak, integrated 1 um depth against that peak) but fit every
  valid channel in 440-1000 nm and integrate 833-1023 nm, instead of the
  specific channel lists of Table 2. Values are comparable, not identical, to
  the CAT products.
* The nine `VIS_SLOPE`...`BD867N` parameters are adapted from the slope and
  band-depth conventions used by the Curiosity team for ChemCam passive and
  Mastcam spectra (Johnson et al., 2015; Wellington et al., 2017; Fraeman et
  al., 2020; Horgan et al., 2020). Anchors that fall in the VNIR dropped-band
  gaps were moved to the nearest clean channel, so they are not identical to
  the published rover parameters. They were not re-checked against those
  papers' tables for this release.

Median kernels: Table 2 gives widths in channels for CRISM's 6.55 nm
sampling. A kernel of k channels is applied as the channels within
(k // 2 + 0.5) x 6.55 nm of the channel nearest the requested wavelength,
which equals the k-channel window on a regular grid but stops at a
dropped-band gap instead of reaching across it. A wavelength more than 8 nm
from any channel gives NaN.

## Known limitations

| Topic | State |
|---|---|
| Spectral smile | TRR3 per-column wavelengths are averaged to one axis; parameters are not evaluated per column |
| Corrections | no photometric, atmospheric or thermal correction; no column ratioing or denoising |
| Multispectral products (MSP, HSP, MRDR) | not tested; the kernel logic assumes hyperspectral sampling |
| Table 2 variants | retired parameters (`BD530`, `OLINDEX`, `BD1750`, ...) are not provided |
| Validation | tests use synthetic cubes with known absorptions; no comparison against released CAT parameter products yet |

## Citing

If you use this package, please cite Viviano-Beck et al. (2014) for the
parameter definitions and this repository (see `CITATION.cff`).

* Viviano-Beck, C. E., et al. (2014). Revised CRISM spectral parameters and summary products based on the currently detected mineral diversity on Mars. Journal of Geophysical Research: Planets, 119, 1403-1431. https://doi.org/10.1002/2014JE004627
* Pelkey, S. M., et al. (2007). CRISM multispectral summary products: Parameterizing mineral diversity on Mars from reflectance. Journal of Geophysical Research, 112, E08S14. https://doi.org/10.1029/2006JE002831
* Murchie, S., et al. (2007). Compact Reconnaissance Imaging Spectrometer for Mars (CRISM) on Mars Reconnaissance Orbiter (MRO). Journal of Geophysical Research, 112, E05S03. https://doi.org/10.1029/2006JE002682
* Murchie, S. L., et al. (2009). Compact Reconnaissance Imaging Spectrometer for Mars investigation and data set from the Mars Reconnaissance Orbiter's primary science phase. Journal of Geophysical Research, 114, E00D07. https://doi.org/10.1029/2009JE003344
* Johnson, J. R., et al. (2015). ChemCam passive reflectance spectroscopy of surface materials at the Curiosity landing site, Mars. Icarus, 249, 74-92. https://doi.org/10.1016/j.icarus.2014.02.028
* Wellington, D. F., et al. (2017). Visible to near-infrared MSL/Mastcam multispectral imaging: Initial results from select high-interest science targets within Gale Crater, Mars. American Mineralogist, 102(6), 1202-1217. https://doi.org/10.2138/am-2017-5760CCBY
* Fraeman, A. A., et al. (2020). Evidence for a diagenetic origin of Vera Rubin ridge, Gale crater, Mars: Summary and synthesis of Curiosity's exploration campaign. Journal of Geophysical Research: Planets, 125, e2020JE006527. https://doi.org/10.1029/2020JE006527
* Horgan, B. H. N., et al. (2020). Diagenesis of Vera Rubin ridge, Gale crater, Mars, from Mastcam multispectral images. Journal of Geophysical Research: Planets, 125, e2019JE006322. https://doi.org/10.1029/2019JE006322

## Acknowledgement

Much of the code was written with an AI coding assistant (Claude, Anthropic)
under the direction of the author, who reviewed the code and results.

## License

MIT, see `LICENSE`.
