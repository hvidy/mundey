"""Tests for PEP 621 packaging, clean imports, XML resources, and SIMBAD offline safety."""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import jax
import numpy as np
from astropy.io import fits
from astropy.table import Table

import mundey
from mundey.mundey import (
    _RANK_WARNING,
    _maybe_display_union_link,
    _query_simbad_target,
    _resolve_target_name,
    load_undershoot_xml,
)


def test_clean_import_has_no_banned_modules_or_global_x64() -> None:
    """Importing mundey must not load IPython/pkg_resources/distutils or enable global x64."""
    assert "IPython" not in sys.modules
    assert "pkg_resources" not in sys.modules
    assert "numpy.distutils" not in sys.modules
    assert jax.config.jax_enable_x64 is False

    pkg_dir = Path(mundey.__file__).resolve().parent
    forbidden_top_level = {
        "IPython",
        "pkg_resources",
        "distutils",
        "fcntl",
        "tessbkgd",
        "halophot",
    }
    for py_file in sorted(pkg_dir.glob("*.py")):
        tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
        for stmt in tree.body:
            if isinstance(stmt, ast.Import):
                for alias in stmt.names:
                    root = alias.name.split(".")[0]
                    assert root not in forbidden_top_level, f"{py_file.name} imports {alias.name}"
                    assert not alias.name.startswith("numpy.distutils")
            elif isinstance(stmt, ast.ImportFrom) and stmt.module:
                root = stmt.module.split(".")[0]
                assert root not in forbidden_top_level, f"{py_file.name} imports {stmt.module}"
                assert not stmt.module.startswith("numpy.distutils")


def test_import_and_display_succeed_when_ipython_and_pkg_resources_blocked(monkeypatch) -> None:
    """Headless execution must succeed even when IPython and pkg_resources are blocked."""
    monkeypatch.setitem(sys.modules, "IPython", None)
    monkeypatch.setitem(sys.modules, "IPython.core.display", None)
    monkeypatch.setitem(sys.modules, "pkg_resources", None)

    _maybe_display_union_link()
    root = load_undershoot_xml()
    assert root.tag == "undershootData"


def test_undershoot_xml_loads_all_64_channels_via_importlib_resources() -> None:
    """Bundled tess-undershoot.xml must load cleanly and cover 4 cameras x 4 CCDs x 4 outputs."""
    root = load_undershoot_xml()
    channels = root.findall("channel")
    assert len(channels) == 64

    seen: set[tuple[int, int, str]] = set()
    for ch in channels:
        cam = int(ch.attrib["cameraNumber"])
        ccd = int(ch.attrib["ccdNumber"])
        out = ch.attrib["ccdOutput"]
        val = float(ch.attrib["undershoot"])
        assert 0.0 < val < 0.02
        seen.add((cam, ccd, out))

    assert len(seen) == 64
    assert hasattr(mundey, "mundey_tpf")


def test_numpy_rankwarning_compat() -> None:
    """RankWarning lookup must resolve to a Warning subclass on NumPy 1.x and 2.x."""
    assert isinstance(_RANK_WARNING, type)
    assert issubclass(_RANK_WARNING, Warning)
    expected = getattr(getattr(np, "exceptions", np), "RankWarning", UserWarning)
    assert _RANK_WARNING is expected


def test_resolve_target_defaults_offline_when_object_present_or_offline_flag() -> None:
    """SIMBAD must not be queried when OBJECT is in hdu[0].header or offline=True."""
    calls = 0

    class BoomSimbad:
        def __init__(self) -> None:
            nonlocal calls
            calls += 1
            raise AssertionError("Simbad should not be instantiated")

    hdr_with_obj = fits.Header({"OBJECT": "TIC 322899250"})
    assert (
        _resolve_target_name(
            hdr_with_obj,
            101.287,
            -16.716,
            322899250,
            simbad_cls=BoomSimbad,
        )
        == "TIC 322899250"
    )
    assert calls == 0

    hdr_empty = fits.Header()
    assert (
        _resolve_target_name(
            hdr_empty,
            101.287,
            -16.716,
            987654,
            use_simbad=True,
            offline=True,
            simbad_cls=BoomSimbad,
        )
        == "TIC 987654"
    )
    assert calls == 0


def test_simbad_modern_legacy_none_and_network_error_handling() -> None:
    """SIMBAD helper must support V/main_id, legacy FLUX_V/MAIN_ID, None, and network errors."""

    class ModernSimbad:
        def add_votable_fields(self, field: str) -> None:
            assert field == "V"

        def query_region(self, _sky) -> Table:
            return Table({"main_id": ["Faint Star", "* alf CMa"], "V": [8.2, -1.46]})

    class LegacySimbad:
        def add_votable_fields(self, field: str) -> None:
            if field == "V":
                raise ValueError("legacy astroquery expects flux(V)")
            assert field == "flux(V)"

        def query_region(self, _sky) -> Table:
            return Table({"MAIN_ID": [b"* alf CMa", b"Companion"], "FLUX_V": [-1.46, 8.5]})

    class NoneSimbad:
        def add_votable_fields(self, _field: str) -> None:
            pass

        def query_region(self, _sky) -> None:
            return None

    class OfflineErrorSimbad:
        def add_votable_fields(self, _field: str) -> None:
            pass

        def query_region(self, _sky) -> Table:
            raise ConnectionError("Compute node has no outbound network")

    assert _query_simbad_target(101.287, -16.716, simbad_cls=ModernSimbad) == "* alf CMa"
    assert _query_simbad_target(101.287, -16.716, simbad_cls=LegacySimbad) == "* alf CMa"
    assert _query_simbad_target(101.287, -16.716, simbad_cls=NoneSimbad) is None
    assert _query_simbad_target(101.287, -16.716, simbad_cls=OfflineErrorSimbad) is None

    hdr_empty = fits.Header()
    assert (
        _resolve_target_name(
            hdr_empty,
            101.287,
            -16.716,
            322899250,
            use_simbad=True,
            simbad_cls=OfflineErrorSimbad,
        )
        == "TIC 322899250"
    )
