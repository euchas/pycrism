"""Synthetic CRISM-like cubes and raw-file writers for the tests.

Nothing here downloads data. Spectra are a continuum multiplied by
(1 - sum of Gaussian absorptions), with a per-pixel depth factor ``k`` so
that every pixel is different and expected values can be bracketed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

STEP = 6.55  # nominal CRISM sampling in nm


def vnir_grid() -> np.ndarray:
    """S-detector-like axis 371-1013 nm with the two VRDR dropped-band gaps."""
    wl = 371.12 + STEP * np.arange(99)
    gap = ((wl > 390) & (wl < 404)) | ((wl > 644) & (wl < 684))
    return wl[~gap]


def ir_grid() -> np.ndarray:
    """L-detector-like axis 1001-3916 nm, regular."""
    return 1001.35 + STEP * np.arange(446)


def gaussian(wl, centre, sigma):
    return np.exp(-0.5 * ((wl - centre) / sigma) ** 2)


def make_cube(wl, features, lines=4, samples=3, base=0.25, slope=0.0, dtype=np.float32):
    """(lines, samples, bands) cube and the (lines, samples) depth factor k.

    ``features`` is a list of (centre_nm, sigma_nm, depth). Pixel (i, j) has
    reflectance cont(wl) * (1 - k_ij * sum(depth * gaussian)).
    """
    wl = np.asarray(wl, dtype=np.float64)
    k = np.linspace(0.5, 1.5, lines * samples).reshape(lines, samples)
    cont = base + slope * (wl - wl[0])
    absorb = np.zeros_like(wl)
    for centre, sigma, depth in features:
        absorb += depth * gaussian(wl, centre, sigma)
    cube = cont[None, None, :] * (1.0 - k[:, :, None] * absorb[None, None, :])
    return cube.astype(dtype), k


def write_raw_envi(path: Path, cube, wl, interleave="bil", bbl=None,
                   units="Nanometers", byte_order=0, ignore=65535.0):
    """Write a (lines, samples, bands) cube as a raw ENVI file, independently
    of the package's own writer."""
    cube = np.asarray(cube, dtype=np.float32)
    lines, samples, bands = cube.shape
    dt = "<f4" if byte_order == 0 else ">f4"
    if interleave == "bil":
        raw = cube.transpose(0, 2, 1)
    elif interleave == "bsq":
        raw = cube.transpose(2, 0, 1)
    else:
        raw = cube
    np.ascontiguousarray(raw).astype(dt).tofile(path)
    hdr = ["ENVI", "description = {synthetic test cube}",
           f"samples = {samples}", f"lines = {lines}", f"bands = {bands}",
           "header offset = 0", "file type = ENVI Standard", "data type = 4",
           f"interleave = {interleave}", f"byte order = {byte_order}",
           f"data ignore value = {ignore:g}",
           f"wavelength units = {units}",
           "wavelength = {\n " + ",\n ".join(f"{w:.6f}" for w in wl) + "}"]
    if bbl is not None:
        hdr.append("bbl = {" + ", ".join(str(int(b)) for b in bbl) + "}")
    Path(path).with_suffix(".hdr").write_text("\n".join(hdr) + "\n")
    return Path(path)


def write_raw_pds(path: Path, cube, wa_name="wa_test_l.img", extra=""):
    """Write a (lines, samples, bands) cube as a PDS3 LINE_INTERLEAVED PC_REAL
    image with a detached label."""
    cube = np.asarray(cube, dtype="<f4")
    lines, samples, bands = cube.shape
    np.ascontiguousarray(cube.transpose(0, 2, 1)).tofile(path)
    label = f"""PDS_VERSION_ID = PDS3
RECORD_TYPE = FIXED_LENGTH
PRODUCT_ID = "SYNTHETIC_TEST"
OBJECT = FILE
  OBJECT = IMAGE
    LINES = {lines}
    LINE_SAMPLES = {samples}
    BANDS = {bands}
    SAMPLE_TYPE = PC_REAL
    SAMPLE_BITS = 32
    BAND_STORAGE_TYPE = LINE_INTERLEAVED
    MRO:WAVELENGTH_FILE_NAME = "{wa_name.upper()}"
{extra}  END_OBJECT = IMAGE
END_OBJECT = FILE
END
"""
    Path(path).with_suffix(".lbl").write_text(label)
    return Path(path)


@pytest.fixture
def vnir_wl():
    return vnir_grid()


@pytest.fixture
def ir_wl():
    return ir_grid()


@pytest.fixture
def vnir_cube(vnir_wl):
    """Hematite-like VNIR cube: bands at 530 and 860 nm on a reddish slope."""
    cube, k = make_cube(vnir_wl, [(530.0, 30.0, 0.10), (860.0, 35.0, 0.08)],
                        base=0.15, slope=2.0e-4)
    return cube, k


@pytest.fixture
def ir_cube(ir_wl):
    """IR cube with a broad 1 um band, a 1.9 um band, a 2.3 um drop and a 2.5 um band."""
    cube, k = make_cube(ir_wl, [(1050.0, 150.0, 0.15), (1930.0, 30.0, 0.10),
                                (2300.0, 20.0, 0.05), (2500.0, 20.0, 0.05)],
                        base=0.20, slope=-1.0e-5)
    return cube, k
