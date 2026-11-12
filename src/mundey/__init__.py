"""Recalibration of TESS target pixel files for saturated bright stars."""

from __future__ import annotations

from . import archive, calibration, mundey
from ._contract import (
    HALOCVER,
    PIXEL_QUALITY_EXTNAME,
    SMEAR_ERR_COLNAME,
    PixelQuality,
    append_smear_err_column,
    atomic_publish,
    atomic_temp_path,
    make_atomic_tmp_path,
    read_mundey_provenance,
    read_pixel_quality,
    read_smear_err_column,
    resolve_cache_dir,
    write_mundey_provenance,
    write_pixel_quality,
)
from ._precision import cast_tree, in_float64, run_in
from .calibration import desmear_calibrated_cutout, recalibrate_tpf
from .mundey import mundey_tpf

name = "mundey"
__version__ = "2.0.0a0"

__all__ = [
    "HALOCVER",
    "PIXEL_QUALITY_EXTNAME",
    "SMEAR_ERR_COLNAME",
    "PixelQuality",
    "__version__",
    "append_smear_err_column",
    "archive",
    "atomic_publish",
    "atomic_temp_path",
    "calibration",
    "cast_tree",
    "desmear_calibrated_cutout",
    "in_float64",
    "make_atomic_tmp_path",
    "mundey",
    "mundey_tpf",
    "name",
    "read_mundey_provenance",
    "read_pixel_quality",
    "read_smear_err_column",
    "recalibrate_tpf",
    "resolve_cache_dir",
    "run_in",
    "write_mundey_provenance",
    "write_pixel_quality",
]
