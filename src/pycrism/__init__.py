"""pycrism: CRISM summary parameters, readers and browse products in pure NumPy."""
from ._version import __version__
from .browse import RECIPES, browse_rgb, stretch
from .parameters import (
    IR_PARAMS,
    MSL_PARAMS,
    PARAM_INFO,
    PARAMS,
    VNIR_PARAMS,
    as_bands_last,
    band_depth,
    compute,
    continuum_line,
    continuum_removed,
    default_names,
    reflectance,
    relative_band_depth,
    shoulder_height,
    spectral_ratio,
    spectral_slope,
    stack_maps,
)
from .readers import (
    parse_envi_header,
    read_cube,
    read_ddr,
    read_envi,
    read_label,
    read_pds_cube,
    read_wa_cdr,
    write_envi,
)

__all__ = [
    "__version__", "RECIPES", "browse_rgb", "stretch",
    "IR_PARAMS", "MSL_PARAMS", "PARAM_INFO", "PARAMS", "VNIR_PARAMS",
    "as_bands_last", "band_depth", "compute", "continuum_line", "continuum_removed",
    "default_names", "reflectance", "relative_band_depth", "shoulder_height",
    "spectral_ratio", "spectral_slope", "stack_maps",
    "parse_envi_header", "read_cube", "read_ddr", "read_envi", "read_label",
    "read_pds_cube", "read_wa_cdr", "write_envi",
]
