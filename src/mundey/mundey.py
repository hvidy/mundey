"""Compatibility ``mundey_tpf`` shim wired onto the ``mundey`` v2.0 JAX calibration engine."""

from __future__ import annotations

import functools
import importlib.resources
import os
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Callable
from pathlib import Path
from typing import Any, TypeVar
from xml.etree.ElementTree import Element

import astropy.coordinates as coord
import astropy.units as u
import lightkurve
import numpy as np
from astropy.io import fits

from .calibration import (
    CadenceMatchMode,
    MissingCadenceMode,
    SmearMode,
    _recalibrate_hdulist_inplace,
)

F = TypeVar("F", bound=Callable[..., Any])

_RANK_WARNING: type[Warning] = getattr(
    getattr(np, "exceptions", np), "RankWarning", UserWarning
)


def _deprecated(message: str) -> Callable[[F], F]:
    """PEP 702 ``@deprecated`` decorator compatible with Python 3.11+."""
    stdlib_dep = getattr(warnings, "deprecated", None)
    if stdlib_dep is not None:
        return stdlib_dep(message, category=DeprecationWarning, stacklevel=2)  # type: ignore[no-any-return]

    def decorator(fn: F) -> F:
        @functools.wraps(fn)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return fn(*args, **kwargs)

        return wrapper  # type: ignore[return-value]

    return decorator


def load_undershoot_xml() -> Element:
    """Load the bundled TESS LDE undershoot XML table via ``importlib.resources``."""
    resource = importlib.resources.files("mundey").joinpath("data/tess-undershoot.xml")
    with resource.open("rb") as fh:
        return ET.parse(fh).getroot()


def _maybe_display_union_link() -> None:
    """Display the union link only when running in an interactive IPython shell."""
    try:
        from IPython import get_ipython
        from IPython.core.display import HTML, display
    except ImportError:
        return
    if get_ipython() is not None:
        display(
            HTML(
                '<a href="https://www.australianunions.org.au/join" '
                'target="_blank" rel="noreferrer noopener">'
                "Join an Australian Union</a> or "
                '<a href="https://www.ituc-csi.org/?page=abook" '
                'target="_blank" rel="noreferrer noopener">'
                "find an International Trade Union Confederation affiliate in your country</a>."
            )
        )


def _query_simbad_target(
    ra: float,
    dec: float,
    *,
    simbad_cls: Any = None,
) -> str | None:
    """Query SIMBAD across legacy (``FLUX_V``/``MAIN_ID``) and modern (``V``/``main_id``) TAP."""
    if ra is None or dec is None:
        return None
    try:
        if simbad_cls is None:
            from astroquery.simbad import Simbad as cls
        else:
            cls = simbad_cls
        custom_simbad = cls()
        try:
            custom_simbad.add_votable_fields("V")
        except Exception:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                custom_simbad.add_votable_fields("flux(V)")

        sky = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
        result_table = custom_simbad.query_region(sky)
        if result_table is None or len(result_table) == 0:
            return None

        colnames = list(result_table.colnames)
        for v_col in ("V", "FLUX_V", "v", "flux_v"):
            if v_col in colnames:
                result_table.sort([v_col])
                break

        for id_col in ("main_id", "MAIN_ID"):
            if id_col in colnames:
                raw_id = result_table[0][id_col]
                if isinstance(raw_id, bytes):
                    raw_id = raw_id.decode("utf-8", errors="replace")
                resolved = str(raw_id).strip()
                if resolved:
                    return resolved
    except Exception:
        return None
    return None


def _resolve_target_name(
    hdu0_header: fits.Header,
    ra: float | None,
    dec: float | None,
    targetid: Any,
    *,
    target: str | None = None,
    use_simbad: bool | None = None,
    offline: bool = False,
    simbad_cls: Any = None,
) -> str:
    """Resolve a human-readable target label without requiring network access on HPC.

    Defaults ``use_simbad`` to ``False`` whenever ``OBJECT`` is already populated
    in ``hdu[0].header`` or ``offline=True``.
    """
    if target is not None:
        return str(target)

    header_object = str(hdu0_header.get("OBJECT", "")).strip()
    fallback = header_object if header_object else f"TIC {targetid}"

    if offline:
        should_query_simbad = False
    elif use_simbad is None:
        should_query_simbad = not bool(header_object)
    else:
        should_query_simbad = bool(use_simbad)

    if should_query_simbad and ra is not None and dec is not None:
        resolved = _query_simbad_target(ra, dec, simbad_cls=simbad_cls)
        if resolved is not None:
            return resolved

    return fallback


class mundey_tpf(lightkurve.TessTargetPixelFile):
    """Compatibility ``TessTargetPixelFile`` subclass delegating to ``mundey`` v2.0."""

    @_deprecated(
        "mundey_tpf.calibrate is deprecated in mundey v2.0; "
        "prefer functional mundey.recalibrate_tpf(tpf, ...)."
    )
    def calibrate(
        self,
        verbose: bool = True,
        smear: SmearMode = "automatic",
        ddir: str | Path | None = "./",
        cache_dir: str | Path | None = None,
        target: str | None = None,
        use_simbad: bool | None = None,
        cadence_match: CadenceMatchMode = "cadenceno",
        missing_cadence: MissingCadenceMode = "nan",
        offline: bool = False,
        strict_hpc: bool = True,
    ) -> None:
        """Recalibrate ``self.hdu`` in place via the ``mundey`` v2.0 JAX engine."""
        resolved_target = _resolve_target_name(
            self.hdu[0].header,
            self.ra,
            self.dec,
            self.targetid,
            target=target,
            use_simbad=use_simbad,
            offline=offline,
        )
        if ddir is not None:
            self.ddir = f"{Path(ddir)}/"

        if verbose:
            print(f"Beginning calibration of {resolved_target}")
            _maybe_display_union_link()

        raw_path = getattr(self, "path", None)
        src_path = (
            Path(raw_path).expanduser()
            if isinstance(raw_path, (str, os.PathLike)) and os.fspath(raw_path)
            else None
        )

        _recalibrate_hdulist_inplace(
            self.hdu,
            source_path=src_path,
            cache_dir=cache_dir,
            ddir=ddir,
            smear=smear,
            cadence_match=cadence_match,
            missing_cadence=missing_cadence,
            offline=offline,
            strict_hpc=strict_hpc,
            verbose=verbose,
            target_obj=self,
        )
