"""Recalibration of TESS target pixel files for saturated bright stars."""

from __future__ import annotations

from ._contract import PixelQuality, read_pixel_quality, read_smear_err
from .calibration import CalibrationConfig, desmear_calibrated_cutout, recalibrate_tpf
from .mundey import mundey_tpf
from .plotting import plot_diagnostics

__version__ = "2.0.1"

__all__ = [
    "CalibrationConfig",
    "PixelQuality",
    "__version__",
    "desmear_calibrated_cutout",
    "mundey_tpf",
    "plot_diagnostics",
    "read_pixel_quality",
    "read_smear_err",
    "recalibrate_tpf",
]
