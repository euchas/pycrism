import numpy as np
import pytest

from conftest import make_cube, write_raw_envi, write_raw_pds
from pycrism import readers as R


def test_parse_envi_header_multiline_lists(tmp_path):
    hdr = tmp_path / "x.hdr"
    hdr.write_text("ENVI\nsamples = 3\nlines = 2\nbands = 2\n"
                   "wavelength = {\n 500.0,\n 600.0 }\n"
                   "band names = {a, b}\ndata type = 4\n")
    meta = R.parse_envi_header(hdr)
    assert meta["samples"] == "3"
    assert meta["wavelength"] == ["500.0", "600.0"]
    assert meta["band names"] == ["a", "b"]
    assert R.envi_wavelengths(meta).tolist() == [500.0, 600.0]


@pytest.mark.parametrize("interleave", ["bil", "bsq", "bip"])
def test_read_envi_interleaves_nulls_and_bad_bands(tmp_path, vnir_wl, interleave):
    cube, _ = make_cube(vnir_wl, [(530.0, 30.0, 0.1)])
    cube[2, 1, :] = 65535.0                      # null pixel
    bbl = np.ones(vnir_wl.size, dtype=int)
    bbl[5] = 0                                   # one bad band
    path = write_raw_envi(tmp_path / f"cube_{interleave}.img", cube, vnir_wl, interleave, bbl=bbl)
    got, wl, meta = R.read_envi(path)
    assert got.shape == cube.shape and got.dtype == np.float32
    assert np.allclose(wl, vnir_wl)
    assert np.isnan(got[2, 1, :]).all()
    assert np.isnan(got[:, :, 5]).all()
    keep = np.ones(cube.shape, dtype=bool)
    keep[2, 1, :] = False
    keep[:, :, 5] = False
    assert np.array_equal(got[keep], cube[keep])


def test_read_envi_micrometre_axis_and_big_endian(tmp_path, vnir_wl):
    cube, _ = make_cube(vnir_wl, [])
    path = write_raw_envi(tmp_path / "um.img", cube, vnir_wl / 1000.0, "bil",
                          units="Micrometers", byte_order=1)
    got, wl, _ = R.read_envi(path)
    assert np.allclose(wl, vnir_wl)
    assert np.allclose(got, cube)


def test_write_envi_round_trip(tmp_path, vnir_wl):
    maps = np.random.default_rng(0).normal(size=(3, 4, 2)).astype(np.float32)
    maps[0, 0, 1] = np.nan
    path = R.write_envi(tmp_path / "out.img", maps, band_names=["BD530_2", "SH600_2"])
    got, wl, meta = R.read_envi(path)
    assert wl is None
    assert meta["band names"] == ["BD530_2", "SH600_2"]
    assert np.isnan(got[0, 0, 1]) and np.allclose(got[np.isfinite(got)], maps[np.isfinite(maps)])
    got2, wl2, _ = R.read_cube(path)
    assert np.array_equal(np.isfinite(got2), np.isfinite(got))


def test_read_pds_cube_with_wa_cdr_flips_descending_axis(tmp_path):
    wl_desc = np.linspace(2000.0, 1000.0, 10)             # stored high to low
    cube, _ = make_cube(wl_desc, [(1500.0, 60.0, 0.1)], lines=4, samples=3)
    img = write_raw_pds(tmp_path / "frt_test_l_trr3.img", cube)
    cols = np.tile(wl_desc, (3, 1)) + np.array([[-0.5], [0.0], [0.5]])   # (samples, bands) smile
    write_raw_pds(tmp_path / "wa_test_l.img", cols[None, :, :])
    got, wl, meta = R.read_cube(img)
    assert got.shape == cube.shape
    assert np.all(np.diff(wl) > 0)
    assert np.allclose(wl, wl_desc[::-1])
    assert np.allclose(got, cube[..., ::-1])
    assert meta["wl_columns"].shape == (3, 10)
    assert np.all(np.diff(meta["wl_columns"][0]) > 0)
    assert meta["PRODUCT_ID"] == "SYNTHETIC_TEST"
