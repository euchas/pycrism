"""Browse-product recipes (RGB triplets of summary parameters) and stretching.

The recipes are the CRISM standard browse products of Viviano-Beck et al.
(2014), Table 4, as delivered with the map-projected MTRDR and VRDR
products: each channel of the RGB image is one summary parameter.
"""
from __future__ import annotations

import numpy as np

from .parameters import compute

RECIPES: dict[str, tuple[str, str, str]] = {
    # VNIR (S detector)
    "TRU": ("R600", "R530", "R440"),
    "VNA": ("R770", "R770", "R770"),
    "FEM": ("BD530_2", "SH600_2", "BDI1000VIS"),
    "FM2": ("BD530_2", "BD920_2", "BDI1000VIS"),
    # IR (L detector)
    "FAL": ("R2529", "R1506", "R1080"),
    "MAF": ("OLINDEX3", "LCPINDEX2", "HCPINDEX2"),
    "HYD": ("SINDEX2", "BD2100_2", "BD1900_2"),
    "PHY": ("D2300", "D2200", "BD1900r2"),
    "PFM": ("BD2355", "D2300", "BD2290"),
    "PAL": ("BD2210_2", "BD2190", "BD2165"),
    "HYS": ("MIN2250", "BD2250", "BD1900r2"),
    "ICE": ("BD1900_2", "BD1500_2", "BD1435"),
    "IC2": ("R3920", "BD1500_2", "BD1435"),
    "CHL": ("ISLOPE1", "BD3000", "IRR2"),
    "CAR": ("D2300", "BD2500_2", "BD1900_2"),
    "CR2": ("MIN2295_2480", "MIN2345_2537", "CINDEX2"),
}

RECIPE_TITLES: dict[str, str] = {
    "TRU": "true colour", "VNA": "VNIR albedo", "FEM": "ferric minerals",
    "FM2": "ferric minerals 2", "FAL": "false colour", "MAF": "mafic minerals",
    "HYD": "hydrated minerals", "PHY": "phyllosilicates",
    "PFM": "Fe/Mg phyllosilicates", "PAL": "Al phyllosilicates",
    "HYS": "hydrated silica", "ICE": "ices", "IC2": "ices 2",
    "CHL": "chlorides", "CAR": "carbonates", "CR2": "carbonates 2",
}


def stretch(band: np.ndarray, percentiles: tuple[float, float] = (2.0, 98.0)) -> np.ndarray:
    """Linear percentile stretch to 0..1; NaN becomes 0."""
    band = np.asarray(band, dtype=np.float32)
    out = np.zeros(band.shape, dtype=np.float32)
    finite = np.isfinite(band)
    if finite.any():
        lo, hi = np.nanpercentile(band[finite], percentiles)
        if hi > lo:
            out[finite] = np.clip((band[finite] - lo) / (hi - lo), 0.0, 1.0)
    return out


def browse_rgb(spec, wl, recipe: str = "FM2",
               percentiles: tuple[float, float] = (2.0, 98.0),
               maps: dict[str, np.ndarray] | None = None) -> np.ndarray:
    """Percentile-stretched (lines, samples, 3) RGB image for a browse recipe.

    Pass already computed parameter maps as `maps` to avoid recomputing them.
    """
    recipe = recipe.upper()
    if recipe not in RECIPES:
        raise KeyError(f"unknown recipe {recipe!r}; choose from {list(RECIPES)}")
    names = RECIPES[recipe]
    maps = dict(maps or {})
    missing = [n for n in dict.fromkeys(names) if n not in maps]
    if missing:
        maps.update(compute(spec, wl, missing))
    return np.stack([stretch(maps[n], percentiles) for n in names], axis=-1)
