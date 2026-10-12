"""CRISM summary parameters after Viviano-Beck et al. (2014), Table 2.

Pure NumPy. Every wavelength, continuum weight and median-kernel width is
written out in the function that uses it, so each definition can be checked
line by line against the paper: Viviano-Beck, C. E., et al. (2014), Revised
CRISM spectral parameters and summary products based on the currently
detected mineral diversity on Mars, J. Geophys. Res. Planets, 119,
1403-1431, doi:10.1002/2014JE004627 (Table 2, "Updated summary products").

Two families are implemented:

* the VNIR (S-detector) set, computable from the map-projected VNIR
  products (VRDR tiles, 371-1017 nm) as well as from full-range cubes;
* the IR (L-detector) set of every non-retired Table 2 parameter between
  1.0 and 3.9 um, computable from TRR3, TER3 and MTR3 products.

A third, small group (``MSL_PARAMS``) holds visible/near-infrared slope and
ratio parameters adapted from the conventions used by the Curiosity rover
team for ChemCam passive and Mastcam multispectral data. They are not part
of Viviano-Beck et al. (2014); see the docstrings for the band positions.

Layout and units
----------------
* Spectra are arrays with the band axis LAST: one spectrum (bands,), a list
  of spectra (n, bands) or a cube (lines, samples, bands). A cube given as
  (bands, lines, samples) is accepted and moved to bands-last.
* ``wl`` is a 1-D array of band-centre wavelengths in nanometres, ascending.
  Null wavelengths must be NaN. A ``wl`` given in micrometres is converted.
* Reflectance is corrected I/F (0 to about 1). Nothing here applies
  photometric or atmospheric corrections.
* Results keep the leading shape of the input, in the input's float type,
  with NaN wherever an input value is NaN.
* The cube is never converted to float64 as a whole: each parameter works on
  a window of a few bands, and the whole-spectrum products (RPEAK1,
  BDI1000VIS, VAR, BDI1000IR, BDI2000) run in blocks of pixels.

Median kernels
--------------
Table 2 gives each kernel width in channels for CRISM's nominal 6.55 nm
sampling. A kernel of k channels is applied here as the channels whose
centre lies within (k // 2 + 0.5) * 6.55 nm of the channel nearest to the
requested wavelength. On a regular grid that is exactly the k-channel window;
across a dropped-band gap (the VNIR gaps near 390-404 and 644-684 nm) the
window stops at the gap instead of reaching across it. If the nearest channel
is more than ``DEFAULT_MAX_GAP_NM`` from the requested wavelength, the term
is NaN for that product.

Continuum conventions (Table 2)
-------------------------------
  band depth       BD = 1 - R_c / (a R_s + b R_l),  b = (W_c - W_s)/(W_l - W_s),  a = 1 - b
  shoulder height  SH = 1 - (a R_s + b R_l) / R_c       (inverted band depth)
  slope-anchored   a straight line through two distant anchors, evaluated at
                   each band: RB = (RC - R)/|RC| (relative band depth) or
                   CR = R/RC (continuum removed)
"""
from __future__ import annotations

import warnings
from typing import Callable

import numpy as np

CRISM_NM_PER_CHANNEL = 6.55     # nominal CRISM spectral sampling
DEFAULT_MAX_GAP_NM = 8.0        # farther than this from any channel -> NaN
_CHUNK = 65536                  # pixels per block for whole-spectrum products
_trapz = getattr(np, "trapezoid", None) or np.trapz


# --------------------------------------------------------------- primitives


def to_nm(wl) -> np.ndarray:
    """Wavelength axis as float64 nm (converted from um if all values < 100)."""
    wl = np.array(wl, dtype=np.float64)
    finite = wl[np.isfinite(wl)]
    if finite.size and finite.max() < 100.0:
        wl = wl * 1000.0
    return wl


def as_bands_last(spec, wl):
    """Return (spec, wl) with the band axis last and ``wl`` in nm."""
    wl = to_nm(wl)
    spec = np.asarray(spec)
    if spec.ndim == 0:
        raise ValueError("spec must have at least one axis")
    if spec.shape[-1] == wl.size:
        return spec, wl
    if spec.ndim == 3 and spec.shape[0] == wl.size:
        return np.moveaxis(spec, 0, -1), wl
    raise ValueError(f"no axis of spec {spec.shape} matches the {wl.size} wavelengths")


def nearest_channel(wl, target_nm: float, max_gap_nm: float = DEFAULT_MAX_GAP_NM) -> int:
    """Index of the channel nearest to ``target_nm``, or -1 if none is within ``max_gap_nm``."""
    wl = np.asarray(wl, dtype=np.float64)
    d = np.abs(wl - float(target_nm))
    d[~np.isfinite(d)] = np.inf
    i = int(np.argmin(d))
    return i if d[i] <= max_gap_nm else -1


def _nan_plane(spec: np.ndarray) -> np.ndarray:
    dtype = spec.dtype if np.issubdtype(spec.dtype, np.floating) else np.float32
    return np.full(spec.shape[:-1], np.nan, dtype=dtype)


def reflectance(spec, wl, target_nm: float, kernel: int = 1,
                max_gap_nm: float = DEFAULT_MAX_GAP_NM) -> np.ndarray:
    """Median reflectance over a ``kernel``-channel window centred on ``target_nm``.

    ``spec`` is bands-last, ``wl`` in nm. Returns an array with the leading
    shape of ``spec``; all-NaN if the wavelength is not covered.
    """
    spec = np.asarray(spec)
    wl = np.asarray(wl, dtype=np.float64)
    i = nearest_channel(wl, target_nm, max_gap_nm)
    if i < 0:
        return _nan_plane(spec)
    h = max(int(kernel) // 2, 0)
    if h == 0:
        return spec[..., i]
    tol = (h + 0.5) * CRISM_NM_PER_CHANNEL
    sel = np.flatnonzero(np.abs(wl - wl[i]) <= tol)
    if sel.size <= 1:
        return spec[..., i]
    if sel[-1] - sel[0] + 1 == sel.size:
        window = spec[..., sel[0]:sel[-1] + 1]
    else:
        window = spec[..., sel]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(window, axis=-1)


def _safe_div(num, den):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den != 0.0, num / den, np.nan)


def continuum_weights(ws: float, wc: float, wl_: float) -> tuple[float, float]:
    """Lever-rule weights (a on the short shoulder, b on the long shoulder)."""
    b = (wc - ws) / (wl_ - ws)
    return 1.0 - b, b


def band_depth(spec, wl, ws: float, wc: float, wl_: float,
               kernels: tuple[int, int, int] = (5, 5, 5)) -> np.ndarray:
    """Band depth, Table 2 convention: BD = 1 - R_c / (a R_s + b R_l)."""
    a, b = continuum_weights(ws, wc, wl_)
    rs = reflectance(spec, wl, ws, kernels[0])
    rc = reflectance(spec, wl, wc, kernels[1])
    rl = reflectance(spec, wl, wl_, kernels[2])
    return 1.0 - _safe_div(rc, a * rs + b * rl)


def shoulder_height(spec, wl, ws: float, wc: float, wl_: float,
                    kernels: tuple[int, int, int] = (5, 5, 5)) -> np.ndarray:
    """Shoulder height (inverted band depth): SH = 1 - (a R_s + b R_l) / R_c."""
    a, b = continuum_weights(ws, wc, wl_)
    rs = reflectance(spec, wl, ws, kernels[0])
    rc = reflectance(spec, wl, wc, kernels[1])
    rl = reflectance(spec, wl, wl_, kernels[2])
    return 1.0 - _safe_div(a * rs + b * rl, rc)


def continuum_line(spec, wl, anchor_lo: float, anchor_hi: float,
                   kernels: tuple[int, int] = (5, 5)) -> Callable[[float], np.ndarray]:
    """Straight continuum through two anchors -> callable of wavelength (nm)."""
    r_lo = reflectance(spec, wl, anchor_lo, kernels[0])
    r_hi = reflectance(spec, wl, anchor_hi, kernels[1])
    slope = (r_hi - r_lo) / (anchor_hi - anchor_lo)

    def at(w: float) -> np.ndarray:
        return r_lo + slope * (w - anchor_lo)
    return at


def relative_band_depth(spec, wl, cont, w: float, kernel: int = 7) -> np.ndarray:
    """RB = (RC - R)/|RC| against a slope-anchored continuum ``cont``."""
    rc = cont(w)
    return _safe_div(rc - reflectance(spec, wl, w, kernel), np.abs(rc))


def continuum_removed(spec, wl, cont, w: float, kernel: int = 5) -> np.ndarray:
    """CR = R/RC against a slope-anchored continuum ``cont``."""
    return _safe_div(reflectance(spec, wl, w, kernel), cont(w))


def spectral_slope(spec, wl, ws: float, wl_: float,
                   kernels: tuple[int, int] = (5, 5)) -> np.ndarray:
    """(R_long - R_short)/(W_long - W_short) in reflectance per nm (nominal span)."""
    rs = reflectance(spec, wl, ws, kernels[0])
    rl = reflectance(spec, wl, wl_, kernels[1])
    return (rl - rs) / (float(wl_) - float(ws))


def spectral_ratio(spec, wl, num_nm: float, den_nm: float,
                   kernels: tuple[int, int] = (5, 5)) -> np.ndarray:
    """R(num_nm) / R(den_nm)."""
    return _safe_div(reflectance(spec, wl, num_nm, kernels[0]),
                     reflectance(spec, wl, den_nm, kernels[1]))


def _flat(spec):
    spec = np.asarray(spec)
    return spec.reshape(-1, spec.shape[-1]), spec.shape[:-1]


# ------------------------------------------------ VNIR (S detector), Table 2


def r440(spec, wl):
    """0.44 um reflectance (TRU blue; clouds and haze)."""
    return reflectance(spec, wl, 440.0, 5)


def r530(spec, wl):
    """0.53 um reflectance (TRU green)."""
    return reflectance(spec, wl, 530.0, 5)


def r600(spec, wl):
    """0.60 um reflectance (TRU red)."""
    return reflectance(spec, wl, 600.0, 5)


def r770(spec, wl):
    """0.77 um reflectance (VNA, VNIR albedo)."""
    return reflectance(spec, wl, 770.0, 5)


def rbr(spec, wl):
    """Red/blue ratio R770/R440; nanophase ferric oxide."""
    return spectral_ratio(spec, wl, 770.0, 440.0)


def irr1(spec, wl):
    """IR ratio 1, R800/R997; aphelion ice clouds versus dust."""
    return spectral_ratio(spec, wl, 800.0, 997.0)


def bd530_2(spec, wl):
    """0.53 um band depth; fine-grained crystalline hematite."""
    return band_depth(spec, wl, 440.0, 530.0, 614.0, (5, 5, 5))


def sh600_2(spec, wl):
    """0.60 um shoulder height; ferric minerals, textured surfaces."""
    return shoulder_height(spec, wl, 533.0, 600.0, 716.0, (5, 5, 3))


def sh770(spec, wl):
    """0.77 um shoulder height (centre 775 nm); ferric minerals."""
    return shoulder_height(spec, wl, 716.0, 775.0, 860.0, (3, 5, 5))


def bd640_2(spec, wl):
    """0.64 um band depth (centre 624 nm); maghemite. Sensitive to VNIR artefacts."""
    return band_depth(spec, wl, 600.0, 624.0, 760.0, (5, 3, 5))


def bd860_2(spec, wl):
    """0.86 um band depth; crystalline ferric minerals, hematite."""
    return band_depth(spec, wl, 755.0, 860.0, 977.0, (5, 5, 5))


def bd920_2(spec, wl):
    """0.92 um band depth; crystalline ferric minerals, low-Ca pyroxene."""
    return band_depth(spec, wl, 807.0, 920.0, 984.0, (5, 5, 5))


def vnir_polynomial_products(spec, wl, fit_range=(440.0, 1000.0),
                             band_range=(833.0, 1023.0), order: int = 5) -> dict[str, np.ndarray]:
    """RPEAK1 and BDI1000VIS from a per-pixel polynomial fit of the VNIR spectrum.

    Method (after Pelkey et al., 2007, and Viviano-Beck et al., 2014):

    1. fit an ``order``-degree polynomial to reflectance against wavelength
       over all valid channels in ``fit_range`` (least squares per pixel);
    2. RPEAK1 = wavelength (nm) of the maximum of the fitted curve;
    3. the 1 um band continuum is taken as the fitted peak reflectance, a
       flat line at R(RPEAK1);
    4. BDI1000VIS = trapezoidal integral of (1 - R/R_peak) over
       ``band_range``, in micrometres.

    Table 2 lists specific channels for the fit and the integral; this
    implementation uses every valid channel in the two ranges instead, so
    values are comparable but not identical to the CAT products. On a VRDR
    tile the integral stops at the last channel near 1017 nm. Channels that
    are NaN in every pixel (dropped bands) are ignored; pixels with a NaN in
    any remaining channel are NaN. Runs in blocks of ``_CHUNK`` pixels.
    """
    spec, wl = as_bands_last(spec, wl)
    flat, lead = _flat(spec)
    n = flat.shape[0]
    rpeak = np.full(n, np.nan, dtype=np.float32)
    bdi = np.full(n, np.nan, dtype=np.float32)
    out = {"RPEAK1": rpeak.reshape(lead), "BDI1000VIS": bdi.reshape(lead)}
    fit_idx = np.flatnonzero((wl >= fit_range[0]) & (wl <= fit_range[1]))
    band_idx = np.flatnonzero((wl >= band_range[0]) & (wl <= band_range[1]))
    min_fit = order + 3
    if fit_idx.size < min_fit or band_idx.size < 3:
        return out
    # Drop channels that are NaN for every pixel (dropped or bad bands).
    keep_fit = np.zeros(fit_idx.size, dtype=bool)
    keep_band = np.zeros(band_idx.size, dtype=bool)
    for s in range(0, n, _CHUNK):
        blk = flat[s:s + _CHUNK]
        keep_fit |= np.isfinite(blk[:, fit_idx]).any(axis=0)
        keep_band |= np.isfinite(blk[:, band_idx]).any(axis=0)
    fit_idx, band_idx = fit_idx[keep_fit], band_idx[keep_band]
    if fit_idx.size < min_fit or band_idx.size < 3:
        return out
    w_fit, w_band = wl[fit_idx], wl[band_idx]
    w0 = float(w_fit.mean())
    vand = np.vander(w_fit - w0, N=order + 1, increasing=True)
    w_dense = np.linspace(w_fit.min(), w_fit.max(), 400)
    vand_dense = np.vander(w_dense - w0, N=order + 1, increasing=True)
    um_band = w_band / 1000.0
    for s in range(0, n, _CHUNK):
        blk = flat[s:s + _CHUNK]
        r_fit = blk[:, fit_idx].astype(np.float64)
        r_band = blk[:, band_idx].astype(np.float64)
        ok = np.isfinite(r_fit).all(axis=1) & np.isfinite(r_band).all(axis=1)
        if not ok.any():
            continue
        coeffs, *_ = np.linalg.lstsq(vand, r_fit[ok].T, rcond=None)
        fit_dense = vand_dense @ coeffs                       # (400, m)
        peak_idx = np.argmax(fit_dense, axis=0)
        peak_val = fit_dense[peak_idx, np.arange(peak_idx.size)]
        safe = peak_val > 1e-6
        with np.errstate(divide="ignore", invalid="ignore"):
            curve = 1.0 - r_band[ok] / peak_val[:, None]
        integ = _trapz(curve, um_band, axis=-1)
        idx = s + np.flatnonzero(ok)
        rpeak[idx] = np.where(safe, w_dense[peak_idx], np.nan)
        bdi[idx] = np.where(safe, integ, np.nan)
    return out


def rpeak1(spec, wl):
    """Wavelength (nm) of the VNIR reflectance peak from the polynomial fit."""
    return vnir_polynomial_products(spec, wl)["RPEAK1"]


def bdi1000vis(spec, wl):
    """1 um integrated band depth, VNIR side; olivine, pyroxene, Fe glass."""
    return vnir_polynomial_products(spec, wl)["BDI1000VIS"]


# ------------------------------- VNIR slopes and ratios (Curiosity conventions)
# Adapted from the visible/near-infrared slope and band-depth conventions used
# by the Curiosity team for ChemCam passive spectra and Mastcam multispectral
# images (Johnson et al., 2015; Wellington et al., 2017; Fraeman et al., 2020;
# Horgan et al., 2020). Anchors that fall into the VNIR dropped-band gaps of
# the map-projected products (390-404 and 644-684 nm) were moved to the
# nearest clean channel, so magnitudes are not identical to the published
# rover parameters. slope = (R_long - R_short)/(W_long - W_short) in
# reflectance per nm; band depth = 1 - R_c/continuum.


def vis_slope(spec, wl):
    """Visible slope 442 -> 639 nm (gap-safe form of the 440-670 nm slope); oxidation, dust."""
    return spectral_slope(spec, wl, 442.0, 639.0)


def s600_840(spec, wl):
    """Red-to-NIR slope 600 -> 840 nm; broad ferric trend."""
    return spectral_slope(spec, wl, 600.0, 840.0)


def s750_840(spec, wl):
    """NIR downturn slope 750 -> 840 nm; negative for red crystalline hematite."""
    return spectral_slope(spec, wl, 750.0, 840.0)


def red_blue(spec, wl):
    """Red/blue ratio R639/R527 (gap-safe form of R640/R527); above 1 is redder."""
    return spectral_ratio(spec, wl, 639.0, 527.0)


def ratio_600_840(spec, wl):
    """Oxidation ratio R600/R840."""
    return spectral_ratio(spec, wl, 600.0, 840.0)


def nir_ratio(spec, wl):
    """NIR ratio R1013/R751; below 1 ferrous (rising NIR), above 1 ferric."""
    return spectral_ratio(spec, wl, 1013.0, 751.0)


def bd535(spec, wl):
    """0.535 um ferric band depth, continuum 500/600 nm."""
    return band_depth(spec, wl, 500.0, 535.0, 600.0, (5, 5, 5))


def bd867w(spec, wl):
    """0.867 um crystalline hematite band depth, wide continuum 751/1012 nm."""
    return band_depth(spec, wl, 751.0, 867.0, 1012.0, (5, 5, 5))


def bd867n(spec, wl):
    """0.867 um crystalline hematite band depth, narrow continuum 805/908 nm."""
    return band_depth(spec, wl, 805.0, 867.0, 908.0, (5, 5, 5))


# -------------------------------------------------- IR (L detector), Table 2
# Non-retired Table 2 parameters between 1.0 and 3.9 um. Where the printed
# table and the CRISM Analysis Toolkit (CAT v7.4) differ, the choice is noted
# on the function.


def r1080(spec, wl):
    """1.08 um reflectance (FAL browse component)."""
    return reflectance(spec, wl, 1080.0, 5)


def r1330(spec, wl):
    """IR albedo at 1.33 um."""
    return reflectance(spec, wl, 1330.0, 11)


def r1506(spec, wl):
    """1.51 um reflectance (FAL browse component)."""
    return reflectance(spec, wl, 1506.0, 5)


def r2529(spec, wl):
    """2.53 um reflectance (FAL and TAN browse component)."""
    return reflectance(spec, wl, 2529.0, 5)


def r3920(spec, wl):
    """3.92 um reflectance (IC2 browse component)."""
    return reflectance(spec, wl, 3920.0, 5)


def bd1300(spec, wl):
    """1.3 um Fe2+ in plagioclase band. Also fires on HCP and Fe phyllosilicates."""
    return band_depth(spec, wl, 1080.0, 1320.0, 1750.0, (5, 15, 5))


def olindex3(spec, wl):
    """Olivine index: weighted relative band depths across the broad 1 um Fe2+ band."""
    c = continuum_line(spec, wl, 1750.0, 2400.0, (7, 7))
    terms = [(1080, .03), (1152, .03), (1210, .03), (1250, .03),
             (1263, .07), (1276, .07), (1330, .12), (1368, .12),
             (1395, .14), (1427, .18), (1470, .18)]
    return sum(wgt * relative_band_depth(spec, wl, c, w, 7) for w, wgt in terms)


def lcpindex2(spec, wl):
    """Low-Ca pyroxene index: broad absorption centred near 1.81 um."""
    c = continuum_line(spec, wl, 1560.0, 2450.0, (7, 7))
    terms = [(1690, .20), (1750, .20), (1810, .30), (1870, .30)]
    return sum(wgt * relative_band_depth(spec, wl, c, w, 7) for w, wgt in terms)


def hcpindex2(spec, wl):
    """High-Ca pyroxene index: broad absorption centred near 2.12 um."""
    c = continuum_line(spec, wl, 1810.0, 2530.0, (7, 7))
    terms = [(2120, .10), (2140, .10), (2230, .15), (2250, .30),
             (2430, .20), (2460, .15)]
    return sum(wgt * relative_band_depth(spec, wl, c, w, 7) for w, wgt in terms)


def var(spec, wl):
    """1.0-2.3 um spectral variance: summed squared residuals from a line through
    the 1014 and 2287 nm channels over the channels between 1021 and 2253 nm
    (CAT hyperspectral form). NaN anywhere in the range gives NaN."""
    spec, wl = as_bands_last(spec, wl)
    i1, i2 = nearest_channel(wl, 1014.0, np.inf), nearest_channel(wl, 2287.0, np.inf)
    w1, w2 = wl[i1], wl[i2]
    lo, hi = sorted((nearest_channel(wl, 1021.0, np.inf), nearest_channel(wl, 2253.0, np.inf)))
    wls = wl[lo:hi + 1]
    flat, lead = _flat(spec)
    out = np.empty(flat.shape[0], dtype=np.float32)
    for s in range(0, flat.shape[0], _CHUNK):
        blk = flat[s:s + _CHUNK].astype(np.float64)
        r1, r2 = blk[:, i1], blk[:, i2]
        with np.errstate(invalid="ignore", divide="ignore"):
            slope = (r2 - r1) / (w2 - w1)
            pred = r1[:, None] + slope[:, None] * (wls[None, :] - w1)
            out[s:s + _CHUNK] = np.sum((blk[:, lo:hi + 1] - pred) ** 2, axis=-1)
    return out.reshape(lead)


def islope1(spec, wl):
    """Spectral slope 1815 -> 2530 nm in reflectance per um (CAT x1000 convention);
    ferric coatings on dark rock are positive. Uses the nominal 715 nm span."""
    r1815 = reflectance(spec, wl, 1815.0, 5)
    r2530 = reflectance(spec, wl, 2530.0, 5)
    return 1000.0 * (r1815 - r2530) / (2530.0 - 1815.0)


def bd1400(spec, wl):
    """1.4 um H2O / structural OH overtone."""
    return band_depth(spec, wl, 1330.0, 1395.0, 1467.0, (5, 3, 5))


def bd1435(spec, wl):
    """1.435 um band: CO2 ice, some hydrated minerals."""
    return band_depth(spec, wl, 1370.0, 1435.0, 1470.0, (3, 1, 3))


def bd1500_2(spec, wl):
    """1.5 um H2O ice band (CAT hyperspectral form: one band depth at 1525 nm)."""
    return band_depth(spec, wl, 1367.0, 1525.0, 1808.0, (5, 11, 5))


def icer1_2(spec, wl):
    """CO2/H2O ice band-depth ratio: 1 - (1 - BD1435)/(1 - BD1500_2)."""
    d1 = bd1435(spec, wl)
    d2 = bd1500_2(spec, wl)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - (1.0 - d1) / (1.0 - d2)


def bd1750_2(spec, wl):
    """1.75 um band: alunite, gypsum."""
    return band_depth(spec, wl, 1690.0, 1750.0, 1815.0, (5, 3, 5))


def bd1900_2(spec, wl):
    """1.9 um bound-H2O band. CAT hyperspectral form with one centre at 1930 nm;
    the printed two-centre (1930, 1985) average is the multispectral form."""
    return band_depth(spec, wl, 1850.0, 1930.0, 2067.0, (5, 5, 5))


def bd1900r2(spec, wl):
    """1.9 um bound-H2O band, ratio form over six channels each side. No kernels
    (Table 2 prints '-')."""
    c = continuum_line(spec, wl, 1850.0, 2060.0, (1, 1))
    num = sum(continuum_removed(spec, wl, c, w, 1)
              for w in (1908, 1914, 1921, 1928, 1934, 1941))
    den = sum(continuum_removed(spec, wl, c, w, 1)
              for w in (1862, 1869, 1875, 2112, 2120, 2126))
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - num / den


def bdi1000ir(spec, wl):
    """1 um integrated band depth, IR side (1020-1255 nm); crystalline Fe2+ silicates."""
    return _integrated_band_depth(spec, wl, 1020.0, 1255.0)


def bdi2000(spec, wl):
    """2 um integrated band depth (1660-2390 nm); pyroxene."""
    return _integrated_band_depth(spec, wl, 1660.0, 2390.0)


def _bdi_continuum(blk, wl, lo1, hi1, lo2, hi2):
    """Continuum for the integrated band depths (CAT crism_sumutil_bdicont,
    hyperspectral form): the 75th percentile of reflectance over 1330-1875 nm
    at that interval's median wavelength, to the median over 2430-2600 nm at
    its median wavelength. Table 2 prints 'median' for the short anchor;
    CAT v7.4 ships the 75th percentile, which is followed here."""
    with np.errstate(invalid="ignore"):
        p75 = np.nanpercentile(blk[:, lo1:hi1 + 1], 75, axis=-1)
        med2 = np.nanmedian(blk[:, lo2:hi2 + 1], axis=-1)
    w1 = float(np.median(wl[lo1:hi1 + 1]))
    w2 = float(np.median(wl[lo2:hi2 + 1]))
    with np.errstate(invalid="ignore", divide="ignore"):
        slope = (med2 - p75) / (w2 - w1)
    return slope, med2, w2


def _integrated_band_depth(spec, wl, wlo, whi):
    """Trapezoidal integral of (1 - R/continuum) over [wlo, whi], in um."""
    spec, wl = as_bands_last(spec, wl)
    lo, hi = sorted((nearest_channel(wl, wlo, np.inf), nearest_channel(wl, whi, np.inf)))
    lo1, hi1 = sorted((nearest_channel(wl, 1330.0, np.inf), nearest_channel(wl, 1875.0, np.inf)))
    lo2, hi2 = sorted((nearest_channel(wl, 2430.0, np.inf), nearest_channel(wl, 2600.0, np.inf)))
    wls = wl[lo:hi + 1]
    um = wls / 1000.0
    flat, lead = _flat(spec)
    out = np.empty(flat.shape[0], dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for s in range(0, flat.shape[0], _CHUNK):
            blk = flat[s:s + _CHUNK].astype(np.float64)
            slope, med2, w2 = _bdi_continuum(blk, wl, lo1, hi1, lo2, hi2)
            with np.errstate(invalid="ignore", divide="ignore"):
                cont = med2[:, None] + slope[:, None] * (wls[None, :] - w2)
                out[s:s + _CHUNK] = _trapz(1.0 - blk[:, lo:hi + 1] / cont, um, axis=-1)
    return out.reshape(lead)


def bd2100_2(spec, wl):
    """2.1 um shifted H2O band; monohydrated sulfates."""
    return band_depth(spec, wl, 1930.0, 2132.0, 2250.0, (3, 5, 3))


def bd2165(spec, wl):
    """2.165 um Al-OH band; kaolinite group."""
    return band_depth(spec, wl, 2120.0, 2165.0, 2230.0, (5, 3, 3))


def bd2190(spec, wl):
    """2.19 um Al-OH band (centre 2185 nm); beidellite, allophane, imogolite."""
    return band_depth(spec, wl, 2120.0, 2185.0, 2250.0, (5, 3, 3))


def min2200(spec, wl):
    """Lesser of the 2.165 and 2.21 um band depths; kaolinite group."""
    d1 = band_depth(spec, wl, 2120.0, 2165.0, 2350.0, (5, 3, 3))
    d2 = band_depth(spec, wl, 2120.0, 2210.0, 2350.0, (5, 3, 3))
    return np.minimum(d1, d2)


def bd2210_2(spec, wl):
    """2.21 um Al-OH band. Table 2 prints the 2250 nm long shoulder (used here);
    CAT v7.4 ships 2290 nm."""
    return band_depth(spec, wl, 2165.0, 2210.0, 2250.0, (5, 5, 5))


def d2200(spec, wl):
    """2.2 um dropoff; Al-OH minerals."""
    c = continuum_line(spec, wl, 1815.0, 2430.0, (5, 5))
    num = continuum_removed(spec, wl, c, 2210, 3) + continuum_removed(spec, wl, c, 2230, 3)
    den = 2.0 * continuum_removed(spec, wl, c, 2165, 5)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - num / den


def bd2230(spec, wl):
    """2.23 um band (centre 2235 nm); hydroxylated ferric sulfate."""
    return band_depth(spec, wl, 2210.0, 2235.0, 2252.0, (3, 3, 3))


def bd2250(spec, wl):
    """2.25 um broad Al-OH / Si-OH band (centre 2245 nm); opal."""
    return band_depth(spec, wl, 2120.0, 2245.0, 2340.0, (5, 7, 3))


def min2250(spec, wl):
    """Lesser of the 2.21 and 2.265 um band depths; opal doublet."""
    d1 = band_depth(spec, wl, 2165.0, 2210.0, 2350.0, (5, 3, 5))
    d2 = band_depth(spec, wl, 2165.0, 2265.0, 2350.0, (5, 3, 5))
    return np.minimum(d1, d2)


def bd2265(spec, wl):
    """2.265 um band; jarosite, gibbsite, acid-leached nontronite."""
    return band_depth(spec, wl, 2210.0, 2265.0, 2340.0, (5, 3, 5))


def bd2290(spec, wl):
    """2.29 um Mg,Fe-OH band (also 2.292 um CO2 ice)."""
    return band_depth(spec, wl, 2250.0, 2290.0, 2350.0, (5, 5, 5))


def d2300(spec, wl):
    """2.3 um dropoff; hydroxylated Fe/Mg silicates (Mg carbonate also fires).

    The third numerator band is 2330 nm as printed in Table 2.
    """
    c = continuum_line(spec, wl, 1815.0, 2530.0, (5, 5))
    num = (continuum_removed(spec, wl, c, 2290, 3) + continuum_removed(spec, wl, c, 2320, 3)
           + continuum_removed(spec, wl, c, 2330, 3))
    den = (continuum_removed(spec, wl, c, 2120, 5) + continuum_removed(spec, wl, c, 2170, 5)
           + continuum_removed(spec, wl, c, 2210, 5))
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - num / den


def bd2355(spec, wl):
    """2.35 um metal-OH band; chlorite, prehnite, pumpellyite, also carbonate."""
    return band_depth(spec, wl, 2300.0, 2355.0, 2450.0, (5, 5, 5))


def sindex2(spec, wl):
    """Convexity at 2.29 um between the 2.1 and 2.4 um bands (hydrated sulfates).
    A band height: the modelled continuum sits in the numerator."""
    return shoulder_height(spec, wl, 2120.0, 2290.0, 2400.0, (5, 7, 3))


def icer2_2(spec, wl):
    """2.6 um relative band depth, continuum 2456-2530 nm; CO2 ice positive, H2O ice negative."""
    c = continuum_line(spec, wl, 2456.0, 2530.0, (5, 5))
    return relative_band_depth(spec, wl, c, 2600, 5)


def min2295_2480(spec, wl):
    """Mg carbonate: lesser of the 2.295 and 2.48 um overtone depths."""
    d1 = band_depth(spec, wl, 2165.0, 2295.0, 2364.0, (5, 5, 5))
    d2 = band_depth(spec, wl, 2364.0, 2480.0, 2570.0, (5, 5, 5))
    return np.minimum(d1, d2)


def min2345_2537(spec, wl):
    """Ca/Fe carbonate: lesser of the 2.345 and 2.537 um overtone depths."""
    d1 = band_depth(spec, wl, 2250.0, 2345.0, 2430.0, (5, 5, 5))
    d2 = band_depth(spec, wl, 2430.0, 2537.0, 2602.0, (5, 5, 5))
    return np.minimum(d1, d2)


def bd2500_2(spec, wl):
    """Mg carbonate overtone near 2.5 um (centre 2480 nm); pairs with D2300."""
    return band_depth(spec, wl, 2364.0, 2480.0, 2570.0, (5, 5, 5))


def bd2600(spec, wl):
    """2.6 um H2O vapour band."""
    return band_depth(spec, wl, 2530.0, 2600.0, 2630.0, (5, 5, 5))


def bd3000(spec, wl):
    """3 um bound-H2O band; continuum extrapolated as R2530 * (R2530 / R2210)."""
    r2210 = reflectance(spec, wl, 2210.0, 5)
    r2530 = reflectance(spec, wl, 2530.0, 5)
    r3000 = reflectance(spec, wl, 3000.0, 5)
    return 1.0 - _safe_div(r3000, r2530 * _safe_div(r2530, r2210))


def bd3100(spec, wl):
    """3.1 um H2O ice band (centre 3120 nm)."""
    return band_depth(spec, wl, 3000.0, 3120.0, 3250.0, (5, 5, 5))


def bd3200(spec, wl):
    """3.2 um CO2 ice band (centre 3320 nm)."""
    return band_depth(spec, wl, 3250.0, 3320.0, 3390.0, (5, 5, 5))


def bd3400_2(spec, wl):
    """3.4 um carbonate band (centre 3420 nm, between the Ca and Mg doublet)."""
    return band_depth(spec, wl, 3250.0, 3420.0, 3630.0, (10, 15, 10))


def cindex2(spec, wl):
    """Convexity at 3.6 um between the 3.4 and 3.9 um carbonate bands. A band
    height like SINDEX2: the modelled continuum sits in the numerator."""
    return shoulder_height(spec, wl, 3450.0, 3610.0, 3875.0, (9, 11, 7))


def irr2(spec, wl):
    """IR ratio 2, R2530/R2210; aphelion ice clouds versus dust."""
    return spectral_ratio(spec, wl, 2530.0, 2210.0)


def irr3(spec, wl):
    """IR ratio 3, R3500/R3390; aphelion ice clouds versus dust."""
    return spectral_ratio(spec, wl, 3500.0, 3390.0, (7, 7))


# ---------------------------------------------------------------- registry

PARAMS: dict[str, Callable] = {
    # VNIR, Table 2
    "R440": r440, "R530": r530, "R600": r600, "R770": r770,
    "RBR": rbr, "IRR1": irr1,
    "BD530_2": bd530_2, "SH600_2": sh600_2, "SH770": sh770, "BD640_2": bd640_2,
    "BD860_2": bd860_2, "BD920_2": bd920_2, "RPEAK1": rpeak1, "BDI1000VIS": bdi1000vis,
    # VNIR, Curiosity conventions
    "VIS_SLOPE": vis_slope, "S600_840": s600_840, "S750_840": s750_840,
    "RED_BLUE": red_blue, "RATIO_600_840": ratio_600_840, "NIR_RATIO": nir_ratio,
    "BD535": bd535, "BD867W": bd867w, "BD867N": bd867n,
    # IR, Table 2 order
    "BDI1000IR": bdi1000ir, "R1330": r1330, "BD1300": bd1300,
    "OLINDEX3": olindex3, "LCPINDEX2": lcpindex2, "HCPINDEX2": hcpindex2,
    "VAR": var, "ISLOPE1": islope1, "BD1400": bd1400, "BD1435": bd1435,
    "BD1500_2": bd1500_2, "ICER1_2": icer1_2, "BD1750_2": bd1750_2,
    "BD1900_2": bd1900_2, "BD1900r2": bd1900r2, "BDI2000": bdi2000,
    "BD2100_2": bd2100_2, "BD2165": bd2165, "BD2190": bd2190, "MIN2200": min2200,
    "BD2210_2": bd2210_2, "D2200": d2200, "BD2230": bd2230, "BD2250": bd2250,
    "MIN2250": min2250, "BD2265": bd2265, "BD2290": bd2290, "D2300": d2300,
    "BD2355": bd2355, "SINDEX2": sindex2, "ICER2_2": icer2_2,
    "MIN2295_2480": min2295_2480, "MIN2345_2537": min2345_2537, "BD2500_2": bd2500_2,
    "BD3000": bd3000, "BD3100": bd3100, "BD3200": bd3200, "BD3400_2": bd3400_2,
    "CINDEX2": cindex2, "R1080": r1080, "R1506": r1506, "R2529": r2529,
    "BD2600": bd2600, "IRR2": irr2, "IRR3": irr3, "R3920": r3920,
}

# Verification status keys used in PARAM_INFO["status"]:
#   "table2"     wavelengths, weights and kernels as printed in Table 2
#   "table2-cat" Table 2 parameter in the CAT v7.4 hyperspectral form (noted)
#   "method"     follows the method description; channel lists not re-checked
#   "adapted"    not in Viviano-Beck et al. (2014); see the docstring
_T2, _CAT, _METHOD, _ADAPTED = "table2", "table2-cat", "method", "adapted"

_INFO = [
    # name, detector, kind, bands [kernels] in nm, description, status
    ("R440", "S", "reflectance", "440 [5]", "0.44 um reflectance; TRU blue, clouds and haze", _T2),
    ("R530", "S", "reflectance", "530 [5]", "0.53 um reflectance; TRU green", _T2),
    ("R600", "S", "reflectance", "600 [5]", "0.60 um reflectance; TRU red", _T2),
    ("R770", "S", "reflectance", "770 [5]", "0.77 um reflectance; VNIR albedo", _T2),
    ("RBR", "S", "ratio", "770 / 440 [5, 5]", "red/blue ratio; nanophase ferric oxide", _T2),
    ("IRR1", "S", "ratio", "800 / 997 [5, 5]", "ice clouds versus dust", _T2),
    ("BD530_2", "S", "band depth", "440 / 530 / 614 [5, 5, 5]", "crystalline hematite", _T2),
    ("SH600_2", "S", "shoulder height", "533 / 600 / 716 [5, 5, 3]", "ferric minerals, textured surfaces", _T2),
    ("SH770", "S", "shoulder height", "716 / 775 / 860 [3, 5, 5]", "ferric minerals", _T2),
    ("BD640_2", "S", "band depth", "600 / 624 / 760 [5, 3, 5]", "maghemite", _T2),
    ("BD860_2", "S", "band depth", "755 / 860 / 977 [5, 5, 5]", "crystalline hematite", _T2),
    ("BD920_2", "S", "band depth", "807 / 920 / 984 [5, 5, 5]", "crystalline ferric minerals, low-Ca pyroxene", _T2),
    ("RPEAK1", "S", "polynomial peak", "fit 440-1000", "wavelength of the VNIR reflectance peak", _METHOD),
    ("BDI1000VIS", "S", "integrated depth", "833-1023 against the fitted peak", "1 um band, VNIR side; olivine, pyroxene, glass", _METHOD),
    ("VIS_SLOPE", "S", "slope", "442 -> 639 [5, 5]", "visible slope; oxidation, dust", _ADAPTED),
    ("S600_840", "S", "slope", "600 -> 840 [5, 5]", "red-to-NIR slope; ferric trend", _ADAPTED),
    ("S750_840", "S", "slope", "750 -> 840 [5, 5]", "NIR downturn; red crystalline hematite negative", _ADAPTED),
    ("RED_BLUE", "S", "ratio", "639 / 527 [5, 5]", "red/blue oxidation ratio", _ADAPTED),
    ("RATIO_600_840", "S", "ratio", "600 / 840 [5, 5]", "oxidation ratio", _ADAPTED),
    ("NIR_RATIO", "S", "ratio", "1013 / 751 [5, 5]", "ferric (above 1) versus ferrous (below 1)", _ADAPTED),
    ("BD535", "S", "band depth", "500 / 535 / 600 [5, 5, 5]", "ferric band", _ADAPTED),
    ("BD867W", "S", "band depth", "751 / 867 / 1012 [5, 5, 5]", "crystalline hematite, wide continuum", _ADAPTED),
    ("BD867N", "S", "band depth", "805 / 867 / 908 [5, 5, 5]", "crystalline hematite, narrow continuum", _ADAPTED),
    ("BDI1000IR", "L", "integrated depth", "1020-1255, continuum 1330-1875 / 2430-2600", "1 um band, IR side; Fe2+ silicates", _CAT),
    ("R1330", "L", "reflectance", "1330 [11]", "IR albedo", _T2),
    ("BD1300", "L", "band depth", "1080 / 1320 / 1750 [5, 15, 5]", "Fe2+ plagioclase", _T2),
    ("OLINDEX3", "L", "weighted RB", "11 bands 1080-1470, continuum 1750-2400 [7]", "olivine", _T2),
    ("LCPINDEX2", "L", "weighted RB", "1690-1870, continuum 1560-2450 [7]", "low-Ca pyroxene", _T2),
    ("HCPINDEX2", "L", "weighted RB", "2120-2460, continuum 1810-2530 [7]", "high-Ca pyroxene", _T2),
    ("VAR", "L", "variance", "line 1014-2287, residuals 1021-2253", "spectral variance 1.0-2.3 um", _CAT),
    ("ISLOPE1", "L", "slope", "1815 -> 2530 [5, 5], per um", "ferric coatings on dark rock", _T2),
    ("BD1400", "L", "band depth", "1330 / 1395 / 1467 [5, 3, 5]", "H2O, structural OH", _T2),
    ("BD1435", "L", "band depth", "1370 / 1435 / 1470 [3, 1, 3]", "CO2 ice, some hydrated minerals", _T2),
    ("BD1500_2", "L", "band depth", "1367 / 1525 / 1808 [5, 11, 5]", "H2O ice", _CAT),
    ("ICER1_2", "L", "ratio of depths", "BD1435, BD1500_2", "CO2 versus H2O ice", _T2),
    ("BD1750_2", "L", "band depth", "1690 / 1750 / 1815 [5, 3, 5]", "alunite, gypsum", _T2),
    ("BD1900_2", "L", "band depth", "1850 / 1930 / 2067 [5, 5, 5]", "bound H2O", _CAT),
    ("BD1900r2", "L", "ratio form", "1908-1941 over 1862-1875 and 2112-2126, line 1850-2060", "bound H2O", _T2),
    ("BDI2000", "L", "integrated depth", "1660-2390, continuum 1330-1875 / 2430-2600", "2 um band; pyroxene", _CAT),
    ("BD2100_2", "L", "band depth", "1930 / 2132 / 2250 [3, 5, 3]", "monohydrated sulfates", _T2),
    ("BD2165", "L", "band depth", "2120 / 2165 / 2230 [5, 3, 3]", "kaolinite group", _T2),
    ("BD2190", "L", "band depth", "2120 / 2185 / 2250 [5, 3, 3]", "beidellite, allophane, imogolite", _T2),
    ("MIN2200", "L", "min of depths", "2165 and 2210, shoulders 2120 / 2350 [5, 3, 3]", "kaolinite group", _T2),
    ("BD2210_2", "L", "band depth", "2165 / 2210 / 2250 [5, 5, 5]", "Al-OH", _T2),
    ("D2200", "L", "dropoff", "2210, 2230 over 2165, line 1815-2430", "Al-OH", _T2),
    ("BD2230", "L", "band depth", "2210 / 2235 / 2252 [3, 3, 3]", "hydroxylated ferric sulfate", _T2),
    ("BD2250", "L", "band depth", "2120 / 2245 / 2340 [5, 7, 3]", "opal", _T2),
    ("MIN2250", "L", "min of depths", "2210 and 2265, shoulders 2165 / 2350 [5, 3, 5]", "opal doublet", _T2),
    ("BD2265", "L", "band depth", "2210 / 2265 / 2340 [5, 3, 5]", "jarosite, gibbsite", _T2),
    ("BD2290", "L", "band depth", "2250 / 2290 / 2350 [5, 5, 5]", "Mg,Fe-OH; CO2 ice", _T2),
    ("D2300", "L", "dropoff", "2290, 2320, 2330 over 2120, 2170, 2210, line 1815-2530", "Fe/Mg phyllosilicates", _T2),
    ("BD2355", "L", "band depth", "2300 / 2355 / 2450 [5, 5, 5]", "chlorite, prehnite, carbonate", _T2),
    ("SINDEX2", "L", "shoulder height", "2120 / 2290 / 2400 [5, 7, 3]", "hydrated sulfates", _T2),
    ("ICER2_2", "L", "RB", "2600, line 2456-2530 [5]", "CO2 ice positive, H2O ice negative", _T2),
    ("MIN2295_2480", "L", "min of depths", "2165 / 2295 / 2364 and 2364 / 2480 / 2570 [5, 5, 5]", "Mg carbonate", _T2),
    ("MIN2345_2537", "L", "min of depths", "2250 / 2345 / 2430 and 2430 / 2537 / 2602 [5, 5, 5]", "Ca/Fe carbonate", _T2),
    ("BD2500_2", "L", "band depth", "2364 / 2480 / 2570 [5, 5, 5]", "Mg carbonate", _T2),
    ("BD3000", "L", "band depth", "3000 against R2530 * R2530 / R2210 [5]", "bound H2O", _T2),
    ("BD3100", "L", "band depth", "3000 / 3120 / 3250 [5, 5, 5]", "H2O ice", _T2),
    ("BD3200", "L", "band depth", "3250 / 3320 / 3390 [5, 5, 5]", "CO2 ice", _T2),
    ("BD3400_2", "L", "band depth", "3250 / 3420 / 3630 [10, 15, 10]", "carbonate", _T2),
    ("CINDEX2", "L", "shoulder height", "3450 / 3610 / 3875 [9, 11, 7]", "carbonate", _T2),
    ("R1080", "L", "reflectance", "1080 [5]", "FAL blue", _T2),
    ("R1506", "L", "reflectance", "1506 [5]", "FAL green", _T2),
    ("R2529", "L", "reflectance", "2529 [5]", "FAL red", _T2),
    ("BD2600", "L", "band depth", "2530 / 2600 / 2630 [5, 5, 5]", "H2O vapour", _T2),
    ("IRR2", "L", "ratio", "2530 / 2210 [5, 5]", "ice clouds versus dust", _T2),
    ("IRR3", "L", "ratio", "3500 / 3390 [7, 7]", "ice clouds versus dust", _T2),
    ("R3920", "L", "reflectance", "3920 [5]", "IC2 red", _T2),
]

PARAM_INFO: dict[str, dict[str, str]] = {
    name: {"detector": det, "kind": kind, "bands": bands, "description": desc, "status": status}
    for name, det, kind, bands, desc, status in _INFO
}
assert set(PARAM_INFO) == set(PARAMS), "registry and metadata disagree"

VNIR_PARAMS: list[str] = [n for n, i in PARAM_INFO.items() if i["detector"] == "S" and i["status"] != _ADAPTED]
MSL_PARAMS: list[str] = [n for n, i in PARAM_INFO.items() if i["status"] == _ADAPTED]
IR_PARAMS: list[str] = [n for n, i in PARAM_INFO.items() if i["detector"] == "L"]


def default_names(wl) -> list[str]:
    """Parameters whose detector range the wavelength axis covers.

    VNIR (and the adapted slopes) if the axis reaches below 1000 nm, IR if it
    reaches above 1050 nm; a joined S+L axis gets both.
    """
    wl = to_nm(wl)
    finite = wl[np.isfinite(wl)]
    names: list[str] = []
    if finite.size and finite.min() < 1000.0:
        names += VNIR_PARAMS + MSL_PARAMS
    if finite.size and finite.max() > 1050.0:
        names += IR_PARAMS
    return names


def compute(spec, wl, names=None) -> dict[str, np.ndarray]:
    """Compute parameters -> {name: array with the leading shape of ``spec``}.

    ``names`` defaults to every parameter the wavelength axis covers.
    RPEAK1 and BDI1000VIS share one polynomial fit, computed once.
    """
    spec, wl = as_bands_last(spec, wl)
    names = list(names) if names is not None else default_names(wl)
    unknown = [n for n in names if n not in PARAMS]
    if unknown:
        raise KeyError(f"unknown parameter(s) {unknown}; available: {list(PARAMS)}")
    out: dict[str, np.ndarray] = {}
    poly = None
    for name in names:
        if name in ("RPEAK1", "BDI1000VIS"):
            if poly is None:
                poly = vnir_polynomial_products(spec, wl)
            out[name] = poly[name]
        else:
            out[name] = np.asarray(PARAMS[name](spec, wl))
    return out


def stack_maps(maps: dict[str, np.ndarray], names=None) -> tuple[np.ndarray, list[str]]:
    """Stack parameter maps into one (..., n_params) float32 array, bands last."""
    names = list(maps) if names is None else list(names)
    stack = np.stack([np.asarray(maps[n], dtype=np.float32) for n in names], axis=-1)
    return stack, names
