import numpy as np

from conftest import make_cube, write_raw_envi
from pycrism import readers as R
from pycrism.cli import main


def test_cli_list(capsys):
    assert main(["list"]) == 0
    out = capsys.readouterr().out
    assert "BD530_2" in out and "OLINDEX3" in out and "FM2" in out
    assert main(["list", "--detector", "S"]) == 0
    out = capsys.readouterr().out
    assert "BD530_2" in out and "OLINDEX3" not in out


def test_cli_params_writes_envi(tmp_path, vnir_wl, capsys):
    cube, _ = make_cube(vnir_wl, [(530.0, 30.0, 0.1)])
    img = write_raw_envi(tmp_path / "tile.img", cube, vnir_wl, "bsq")
    out = tmp_path / "params.img"
    assert main(["params", str(img), "-p", "BD530_2", "RPEAK1", "R770", "-o", str(out)]) == 0
    text = capsys.readouterr().out
    assert "BD530_2" in text and "valid 100.0 %" in text
    got, wl, meta = R.read_envi(out)
    assert got.shape == (4, 3, 3)
    assert meta["band names"] == ["BD530_2", "RPEAK1", "R770"]
    assert np.isfinite(got).all()


def test_cli_browse_envi_output(tmp_path, vnir_wl):
    cube, _ = make_cube(vnir_wl, [(530.0, 30.0, 0.1), (920.0, 40.0, 0.05)])
    img = write_raw_envi(tmp_path / "tile.img", cube, vnir_wl, "bil")
    out = tmp_path / "tile_FM2.img"
    assert main(["browse", str(img), "fm2", "-o", str(out)]) == 0
    rgb, _, meta = R.read_envi(out)
    assert rgb.shape == (4, 3, 3)
    assert meta["band names"] == ["BD530_2", "BD920_2", "BDI1000VIS"]
    assert np.nanmin(rgb) >= 0.0 and np.nanmax(rgb) <= 1.0
