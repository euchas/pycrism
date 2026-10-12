import numpy as np
import pytest

from conftest import gaussian, make_cube
from pycrism import parameters as P


def test_reflectance_window_is_channel_count_on_regular_grid(ir_wl, ir_cube):
    cube, _ = ir_cube
    i = int(np.argmin(np.abs(ir_wl - 2210.0)))
    expected = np.median(cube[..., i - 2:i + 3], axis=-1)
    got = P.reflectance(cube, ir_wl, 2210.0, 5)
    assert got.shape == cube.shape[:-1]
    assert np.allclose(got, expected)
    assert np.array_equal(P.reflectance(cube, ir_wl, 2210.0, 1), cube[..., i])


def test_reflectance_gap_and_out_of_range_are_nan(vnir_wl, vnir_cube):
    cube, _ = vnir_cube
    assert np.isnan(P.reflectance(cube, vnir_wl, 660.0, 5)).all()      # in the 644-684 gap
    assert np.isnan(P.reflectance(cube, vnir_wl, 1300.0, 5)).all()     # beyond the detector
    assert np.isfinite(P.reflectance(cube, vnir_wl, 637.0, 5)).all()   # last clean channel


def test_window_stops_at_gap(vnir_wl):
    # spectrum equal to its own wavelength; kernel 3 at the last channel before
    # the 644 nm gap must average only the two channels below the gap
    spec = np.tile(vnir_wl.astype(np.float32), (2, 2, 1))
    i = int(np.argmin(np.abs(vnir_wl - 639.0)))
    assert vnir_wl[i] < 644 and vnir_wl[i + 1] > 684
    got = P.reflectance(spec, vnir_wl, 639.0, 3)
    assert np.allclose(got, (vnir_wl[i - 1] + vnir_wl[i]) / 2.0, atol=1e-3)


def test_band_depth_is_zero_on_flat_and_small_on_linear(vnir_wl):
    flat = np.full((2, 2, vnir_wl.size), 0.3, dtype=np.float32)
    assert np.abs(P.bd530_2(flat, vnir_wl)).max() < 1e-6
    assert np.abs(P.sh600_2(flat, vnir_wl)).max() < 1e-6
    lin, _ = make_cube(vnir_wl, [], base=0.2, slope=2e-4)
    assert np.abs(P.bd530_2(lin, vnir_wl)).max() < 5e-3


def test_bd530_2_brackets_the_synthetic_absorption(vnir_wl, vnir_cube):
    cube, k = vnir_cube
    bd = P.bd530_2(cube, vnir_wl)
    expected = 0.10 * k
    assert np.all(bd > 0.85 * expected)
    assert np.all(bd < 1.10 * expected)        # nominal-wavelength continuum on a sloped spectrum
    assert np.all(np.diff(bd.ravel()) > 0)          # deeper with k
    bd860 = P.bd860_2(cube, vnir_wl)
    assert np.all((bd860 > 0.8 * 0.08 * k) & (bd860 < 0.08 * k * 1.10))


def test_shoulder_height_positive_for_a_peak(vnir_wl):
    spec = (0.3 * (1.0 + 0.1 * gaussian(vnir_wl, 600.0, 30.0))).astype(np.float32)
    assert P.sh600_2(spec, vnir_wl) > 0.05
    assert P.band_depth(spec, vnir_wl, 533.0, 600.0, 716.0) < -0.05


def test_slope_and_ratio_on_linear_spectrum(vnir_wl):
    lin, _ = make_cube(vnir_wl, [], base=0.2, slope=1e-4)
    assert np.allclose(P.s600_840(lin, vnir_wl), 1e-4, rtol=0.02)
    # the 639 nm window is cut by the 644-684 nm gap, so this one is less exact
    assert np.allclose(P.vis_slope(lin, vnir_wl), 1e-4, rtol=0.05)
    assert np.all(P.red_blue(lin, vnir_wl) > 1.0)
    assert np.all(P.rbr(lin, vnir_wl) > 1.0)


def test_bands_first_layout_gives_same_result(vnir_wl, vnir_cube):
    cube, _ = vnir_cube
    a = P.compute(cube, vnir_wl, ["BD530_2", "SH770", "RBR"])
    b = P.compute(np.moveaxis(cube, -1, 0), vnir_wl, ["BD530_2", "SH770", "RBR"])
    for n in a:
        assert a[n].shape == cube.shape[:-1]
        assert np.allclose(a[n], b[n], equal_nan=True)


def test_rpeak1_and_bdi1000vis(vnir_wl):
    hump = 0.30 - 4e-7 * (vnir_wl - 750.0) ** 2             # smooth peak at 750 nm
    plain = np.tile(hump.astype(np.float32), (1, 1, 1))
    absorbed = (hump * (1.0 - 0.2 * gaussian(vnir_wl, 1000.0, 80.0))).astype(np.float32)[None, None, :]
    out = P.vnir_polynomial_products(np.concatenate([plain, absorbed], axis=1), vnir_wl)
    rpeak, bdi = out["RPEAK1"][0], out["BDI1000VIS"][0]
    assert abs(rpeak[0] - 750.0) < 15.0
    assert bdi[1] > bdi[0] + 0.005
    assert bdi[1] > 0.0


def test_ir_indices_respond_to_the_right_features(ir_wl, ir_cube):
    cube, k = ir_cube
    maps = P.compute(cube, ir_wl, ["OLINDEX3", "BD1900r2", "BD1900_2", "D2300", "BD2500_2",
                                   "MIN2295_2480", "BDI1000IR", "LCPINDEX2"])
    for n in ("OLINDEX3", "BD1900r2", "BD1900_2", "D2300", "BD2500_2", "MIN2295_2480", "BDI1000IR"):
        assert maps[n].shape == cube.shape[:-1]
        assert np.all(maps[n] > 0.0), n
        assert np.all(np.diff(maps[n].ravel()) > 0), n     # grows with k
    flat = np.full((2, 2, ir_wl.size), 0.3, dtype=np.float32)
    for n, arr in P.compute(flat, ir_wl, ["OLINDEX3", "D2300", "BD1900r2", "BD2500_2", "VAR"]).items():
        assert np.abs(arr).max() < 1e-5, n


def test_compute_defaults_unknown_and_nan_pixels(vnir_wl, vnir_cube, ir_wl):
    cube, _ = vnir_cube
    names = P.default_names(vnir_wl)
    assert "BD530_2" in names and "OLINDEX3" not in names
    names_ir = P.default_names(ir_wl)
    assert "OLINDEX3" in names_ir and "BD530_2" not in names_ir
    with pytest.raises(KeyError):
        P.compute(cube, vnir_wl, ["NOT_A_PARAMETER"])
    cube = cube.copy()
    cube[1, 2, :] = np.nan
    maps = P.compute(cube, vnir_wl)
    assert set(maps) == set(names)
    for n, arr in maps.items():
        assert arr.shape == cube.shape[:-1], n
        assert np.isnan(arr[1, 2]), n
        assert np.isfinite(arr[0, 0]), n


def test_float32_input_is_neither_upcast_nor_modified(vnir_wl, vnir_cube):
    cube32, _ = vnir_cube
    before = cube32.copy()
    bd32 = P.bd530_2(cube32, vnir_wl)
    assert bd32.dtype == np.float32
    assert np.array_equal(cube32, before)
    bd64 = P.bd530_2(cube32.astype(np.float64), vnir_wl)
    assert np.allclose(bd32, bd64, atol=1e-5)


def test_registry_metadata_is_complete():
    assert set(P.PARAM_INFO) == set(P.PARAMS)
    assert set(P.VNIR_PARAMS) | set(P.MSL_PARAMS) | set(P.IR_PARAMS) == set(P.PARAMS)
    assert len(P.VNIR_PARAMS) == 14 and len(P.MSL_PARAMS) == 9 and len(P.IR_PARAMS) == 46
