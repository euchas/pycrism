"""Readers and a writer for CRISM cubes: ENVI (.img + .hdr) and PDS3 (.img + .lbl).

Everything is read with NumPy straight from the raw files, so there is no
dependency on ENVI, the CRISM Analysis Toolkit (CAT) or GDAL.

Conventions
-----------
* cubes are returned as (lines, samples, bands) float32, C-contiguous
* the CRISM null value 65535 and any non-finite value become NaN
* wavelengths are returned in nanometres, ascending; a cube stored in
  descending wavelength order (raw TRR3 L-detector cubes) is flipped
* ENVI bad-band lists (``bbl``) are applied by setting those bands to NaN

Memory
------
The raw file is memory-mapped and converted once into the float32 output,
so the peak footprint is about one float32 copy of the cube.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

import numpy as np

NULL = 65535.0

_ENVI_DTYPES = {1: "u1", 2: "<i2", 3: "<i4", 4: "<f4", 5: "<f8",
                12: "<u2", 13: "<u4", 14: "<i8", 15: "<u8"}


# ---------------------------------------------------------------- ENVI header


def parse_envi_header(hdr_path: str | os.PathLike) -> dict[str, Any]:
    """Parse an ENVI header into a dict with lower-case keys.

    Brace lists (``wavelength = { ... }``, possibly spanning several lines)
    become Python lists of strings; everything else stays a string.
    """
    text = Path(hdr_path).read_text(errors="ignore")
    meta: dict[str, Any] = {}
    for m in re.finditer(r"^\s*([^=\n{}]+?)\s*=\s*(\{.*?\}|[^\n]*)", text, re.M | re.S):
        key = m.group(1).strip().lower()
        val = m.group(2).strip()
        if val.startswith("{"):
            inner = val.strip("{}").replace("\n", " ")
            meta[key] = [v.strip() for v in inner.split(",") if v.strip()]
        else:
            meta[key] = val
    return meta


def envi_wavelengths(meta: dict[str, Any]) -> np.ndarray | None:
    """Band centres in nm from a parsed header (None if the header has none)."""
    if not meta.get("wavelength"):
        return None
    wl = np.array([float(v) for v in meta["wavelength"]], dtype=np.float64)
    units = str(meta.get("wavelength units", "")).lower()
    finite = wl[np.isfinite(wl) & (wl < NULL - 1.0)]
    if units.startswith("micro") or (not units and finite.size and finite.max() < 100.0):
        wl = wl * 1000.0
    wl[wl >= NULL - 1.0] = np.nan
    return wl


def _sibling(img_path: str | os.PathLike, ext: str) -> Path:
    img_path = Path(img_path)
    for cand in (img_path.with_suffix(ext), img_path.with_suffix(ext.upper()),
                 Path(str(img_path) + ext)):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"no {ext} file next to {img_path}")


def read_envi(img_path: str | os.PathLike, hdr_path: str | os.PathLike | None = None):
    """Read an ENVI cube -> (cube (lines, samples, bands) float32, wl nm or None, header dict).

    Honours interleave, byte order, data type, header offset, data ignore
    value and the bad-band list. Works for the ENVI-labelled CRISM TER3,
    MTR3 and map-projected VNIR (VRDR) products.
    """
    hdr_path = Path(hdr_path) if hdr_path else _sibling(img_path, ".hdr")
    meta = parse_envi_header(hdr_path)
    lines, samples, bands = int(meta["lines"]), int(meta["samples"]), int(meta["bands"])
    dt = _ENVI_DTYPES[int(meta.get("data type", 4))]
    if int(meta.get("byte order", 0)) == 1 and dt[0] == "<":
        dt = ">" + dt[1:]
    offset = int(meta.get("header offset", 0))
    interleave = str(meta.get("interleave", "bsq")).lower()
    shapes = {"bil": (lines, bands, samples), "bsq": (bands, lines, samples),
              "bip": (lines, samples, bands)}
    if interleave not in shapes:
        raise ValueError(f"unknown interleave {interleave!r} in {hdr_path}")
    mm = np.memmap(img_path, dtype=np.dtype(dt), mode="r", offset=offset,
                   shape=shapes[interleave])
    if interleave == "bil":
        view = mm.transpose(0, 2, 1)
    elif interleave == "bsq":
        view = mm.transpose(1, 2, 0)
    else:
        view = mm
    cube = np.array(view, dtype=np.float32, order="C")
    del view, mm
    ignore = float(meta.get("data ignore value", NULL))
    cube[(cube == ignore) | (cube >= NULL - 1.0) | ~np.isfinite(cube)] = np.nan
    if meta.get("bbl"):
        bbl = np.array([int(float(x)) for x in meta["bbl"]])
        if bbl.size == bands:
            cube[..., bbl == 0] = np.nan
    wl = envi_wavelengths(meta)
    if wl is not None and wl.size != bands:
        raise ValueError(f"{hdr_path}: {wl.size} wavelengths for {bands} bands")
    cube, wl = _ascending(cube, wl)
    return cube, wl, meta


def write_envi(img_path: str | os.PathLike, data: np.ndarray, wl=None,
               band_names=None, description: str = "", ignore_value: float = NULL):
    """Write a (lines, samples, bands) array as a float32 BSQ ENVI file with a .hdr.

    NaN is written as `ignore_value` and declared as ``data ignore value``.
    """
    img_path = Path(img_path)
    data = np.asarray(data, dtype=np.float32)
    if data.ndim == 2:
        data = data[..., None]
    lines, samples, bands = data.shape
    out = np.where(np.isfinite(data), data, np.float32(ignore_value))
    np.ascontiguousarray(out.transpose(2, 0, 1)).tofile(img_path)
    hdr = ["ENVI", f"description = {{{description or img_path.name}}}",
           f"samples = {samples}", f"lines = {lines}", f"bands = {bands}",
           "header offset = 0", "file type = ENVI Standard", "data type = 4",
           "interleave = bsq", "byte order = 0",
           f"data ignore value = {ignore_value:g}"]
    if band_names is not None:
        hdr.append("band names = {" + ", ".join(str(b) for b in band_names) + "}")
    if wl is not None:
        hdr.append("wavelength units = Nanometers")
        hdr.append("wavelength = {" + ", ".join(f"{float(w):.6g}" for w in wl) + "}")
    img_path.with_suffix(".hdr").write_text("\n".join(hdr) + "\n")
    return img_path


# ---------------------------------------------------------------- PDS3 labels


_LABEL_KEYS = ("LINES", "LINE_SAMPLES", "BANDS", "SAMPLE_BITS", "SAMPLE_TYPE",
               "BAND_STORAGE_TYPE", "PRODUCT_ID", "START_TIME", "MRO:SENSOR_ID",
               "MRO:WAVELENGTH_FILE_NAME", "SOLAR_LONGITUDE", "OBSERVATION_TYPE",
               "MRO:DETECTOR_TEMPERATURE")


def read_label(path: str | os.PathLike) -> dict[str, Any]:
    """Parse the flat KEYWORD = VALUE pairs needed from a PDS3 label."""
    txt = Path(path).read_text(errors="ignore")
    out: dict[str, Any] = {}
    for key in _LABEL_KEYS:
        m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(.+?)\s*$", txt, re.M)
        if m:
            out[key] = m.group(1).strip().strip('"')
    for k in ("LINES", "LINE_SAMPLES", "BANDS", "SAMPLE_BITS"):
        if k in out:
            out[k] = int(out[k])
    return out


def _pds_dtype(lbl: dict[str, Any]) -> np.dtype:
    st = lbl.get("SAMPLE_TYPE", "PC_REAL")
    bits = int(lbl.get("SAMPLE_BITS", 32))
    if "REAL" not in st:
        raise ValueError(f"unsupported SAMPLE_TYPE {st}")
    endian = "<" if st.startswith("PC") else ">"
    return np.dtype(f"{endian}f{bits // 8}")


def read_pds_cube(img_path: str | os.PathLike, lbl_path: str | os.PathLike | None = None):
    """Read a PDS3 CRISM cube (detached label) -> (cube (l, s, b) float32, label dict).

    The image is taken as the last LINES x LINE_SAMPLES x BANDS samples of the
    file, which is the whole file for CRISM products with detached labels.
    """
    lbl_path = Path(lbl_path) if lbl_path else _sibling(img_path, ".lbl")
    lbl = read_label(lbl_path)
    lines, samples, bands = lbl["LINES"], lbl["LINE_SAMPLES"], lbl["BANDS"]
    dt = _pds_dtype(lbl)
    storage = lbl.get("BAND_STORAGE_TYPE", "LINE_INTERLEAVED")
    shapes = {"LINE_INTERLEAVED": (lines, bands, samples),
              "BAND_SEQUENTIAL": (bands, lines, samples),
              "SAMPLE_INTERLEAVED": (lines, samples, bands)}
    if storage not in shapes:
        raise ValueError(f"unsupported BAND_STORAGE_TYPE {storage}")
    expected = lines * samples * bands * dt.itemsize
    size = os.path.getsize(img_path)
    if size < expected:
        raise ValueError(f"{img_path}: {size} bytes, expected at least {expected}")
    mm = np.memmap(img_path, dtype=dt, mode="r", offset=size - expected, shape=shapes[storage])
    if storage == "LINE_INTERLEAVED":
        view = mm.transpose(0, 2, 1)
    elif storage == "BAND_SEQUENTIAL":
        view = mm.transpose(1, 2, 0)
    else:
        view = mm
    cube = np.array(view, dtype=np.float32, order="C")
    del view, mm
    cube[(cube >= NULL - 1.0) | ~np.isfinite(cube)] = np.nan
    return cube, lbl


def read_wa_cdr(wa_img: str | os.PathLike) -> np.ndarray:
    """CRISM WA calibration cube -> per-column wavelengths (nm), shape (samples, bands).

    The WA CDR is a one-line cube giving the centre wavelength of every
    (band, detector column) pair, i.e. it carries the spectral smile.
    """
    cube, _ = read_pds_cube(wa_img)
    wl = cube[0].astype(np.float64)
    wl[(wl <= 0) | ~np.isfinite(wl)] = np.nan
    return wl


DDR_BANDS = {  # CRISM DDR band order (1-indexed) from the BAND_NAME list of the DDR label
    "incidence_areoid": 1, "emission_areoid": 2, "phase": 3,
    "lat": 4, "lon": 5, "incidence": 6, "emission": 7,
    "slope": 8, "slope_azimuth": 9, "elevation": 10,
    "thermal_inertia": 11, "albedo": 12, "local_solar_time": 13, "altitude": 14,
}


def read_ddr(ddr_img: str | os.PathLike) -> dict[str, Any]:
    """Read a CRISM DDR -> dict of named geometry planes, each (lines, samples)."""
    cube, lbl = read_pds_cube(ddr_img)
    out: dict[str, Any] = {"_bands": cube.shape[2], "_label": lbl}
    for k, b in DDR_BANDS.items():
        if b - 1 < cube.shape[2]:
            out[k] = cube[:, :, b - 1]
    return out


# ------------------------------------------------------------------ dispatch


def _descending(wl: np.ndarray | None) -> bool:
    if wl is None:
        return False
    good = np.isfinite(wl)
    return bool(good.sum() > 2 and np.nanmedian(np.diff(wl[good])) < 0)


def _ascending(cube: np.ndarray, wl: np.ndarray | None):
    """Flip the band axis if the wavelength axis runs downwards."""
    if _descending(wl):
        return np.ascontiguousarray(cube[..., ::-1]), wl[::-1].copy()
    return cube, wl


def read_cube(path: str | os.PathLike, wa: str | os.PathLike | None = None):
    """Read a CRISM cube by file type -> (cube (l, s, b) float32, wl (b,) nm, meta dict).

    * ``.img`` with a sibling ``.hdr``: ENVI (TER3, MTR3, VRDR); wavelengths
      from the header.
    * ``.img`` with a sibling ``.lbl``: PDS3 (TRR3 and other raw cubes);
      wavelengths come from the WA CDR named in the label. Pass its path as
      ``wa``, or put it next to the cube (or in a ``CDR`` folder beside it).
      The per-column wavelengths are averaged to one axis; the full
      (samples, bands) array is returned as ``meta["wl_columns"]``.
    """
    path = Path(path)
    try:
        hdr = _sibling(path, ".hdr")
    except FileNotFoundError:
        hdr = None
    if hdr is not None:
        return read_envi(path, hdr)
    cube, lbl = read_pds_cube(path)
    meta: dict[str, Any] = dict(lbl)
    wl = None
    if wa is None and lbl.get("MRO:WAVELENGTH_FILE_NAME"):
        name = lbl["MRO:WAVELENGTH_FILE_NAME"].lower().replace(".img", "")
        for cand in (path.parent / (name + ".img"), path.parent / (name + ".IMG"),
                     path.parent / "CDR" / (name + ".img")):
            if cand.exists():
                wa = cand
                break
    if wa is not None:
        cols = read_wa_cdr(wa)
        if cols.shape[0] != cube.shape[1]:
            raise ValueError(f"{path.name}: WA CDR has {cols.shape[0]} columns, "
                             f"cube has {cube.shape[1]}")
        with np.errstate(invalid="ignore"):
            wl = np.nanmean(cols, axis=0)
        if _descending(wl):
            cols = cols[:, ::-1].copy()
        meta["wl_columns"] = cols
    cube, wl = _ascending(cube, wl)
    return cube, wl, meta
