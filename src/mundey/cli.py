"""Command-line interface for ``mundey`` (``mundey prefetch`` & ``mundey calibrate``)."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from astropy.io import fits

from .archive import DEFAULT_CHUNK_ROWS, CadenceKind, infer_cadence_kind, prefetch_collateral
from .calibration import (
    OUTPUT_ORDER,
    desmear_calibrated_cutout,
    output_for_rawx_col,
    recalibrate_tpf,
)


def _parse_outputs(raw: str) -> tuple[str, ...]:
    tokens = [tok.strip().upper() for tok in raw.replace(",", " ").split() if tok.strip()]
    if not tokens:
        raise ValueError("Expected at least one CCD output ('A', 'B', 'C', 'D').")
    for tok in tokens:
        if tok not in OUTPUT_ORDER:
            raise ValueError(f"Invalid CCD output {tok!r}; expected A, B, C, or D.")
    return tuple(out for out in OUTPUT_ORDER if out in tokens)


def _parse_cadences(raw: str) -> tuple[CadenceKind, ...]:
    val = raw.strip().lower()
    if val == "both":
        return ("2min", "20s")
    if val in ("2min", "20s"):
        return (val,)  # type: ignore[return-value]
    raise ValueError(f"Invalid --cadence {raw!r}; expected '2min', '20s', or 'both'.")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mundey",
        description="TESS bright-star recalibration and HPC collateral prefetching.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- mundey prefetch ---
    p_pre = sub.add_parser(
        "prefetch",
        help="Download static CCD models and build chunked -blf.fits collateral caches.",
    )
    p_pre.add_argument("--sector", type=int, default=None, help="TESS sector number.")
    p_pre.add_argument("--camera", type=int, default=None, help="TESS camera (1..4).")
    p_pre.add_argument("--ccd", type=int, default=None, help="TESS CCD (1..4).")
    p_pre.add_argument(
        "--task-id",
        type=int,
        default=None,
        help="0-based Slurm array index 0..15 mapping to camera=(id//4)+1, ccd=(id%%4)+1.",
    )
    p_pre.add_argument(
        "--tpf",
        nargs="+",
        default=None,
        help="Infer (sector, camera, ccd, outputs, cadence) from one or more TPF FITS files.",
    )
    p_pre.add_argument(
        "--outputs",
        type=str,
        default=None,
        help="Comma-separated CCD outputs to prefetch (default: A,B,C,D or inferred from --tpf).",
    )
    p_pre.add_argument(
        "--cadence",
        type=str,
        default=None,
        help="Collateral cadence: '2min', '20s', or 'both' (default: '2min' or inferred).",
    )
    p_pre.add_argument("--cache-dir", type=str, default=None, help="Cache directory path.")
    p_pre.add_argument("--ddir", type=str, default=None, help="Local reference directory.")
    p_pre.add_argument(
        "--stage-dir",
        type=str,
        default=None,
        help="Temporary staging directory for raw _col.fits (defaults to $JOBFS when set).",
    )
    p_pre.add_argument(
        "--chunk-rows",
        type=int,
        default=DEFAULT_CHUNK_ROWS,
        help=f"Cadence rows per streaming chunk (default: {DEFAULT_CHUNK_ROWS}).",
    )
    p_pre.add_argument(
        "--keep-raw",
        action="store_true",
        help="Retain staged raw _col.fits inside --stage-dir after publishing -blf.fits.",
    )
    p_pre.add_argument(
        "--offline",
        action="store_true",
        help="Forbid all outbound S3 network calls (verify local staging only).",
    )
    p_pre.add_argument(
        "--strict-hpc",
        action="store_true",
        help="Raise ValueError if --cache-dir resolves under /fred.",
    )
    p_pre.add_argument("-v", "--verbose", action="store_true", help="Print progress.")

    # --- mundey calibrate ---
    p_cal = sub.add_parser(
        "calibrate",
        help="Recalibrate a TESS TPF (or desmear an FFI cutout) to a HALOCVER=1 FITS file.",
    )
    p_cal.add_argument("tpf", type=str, help="Input TPF FITS path.")
    p_cal.add_argument("-o", "--output", type=str, required=True, help="Output FITS path.")
    p_cal.add_argument(
        "--smear",
        choices=("automatic", "standard", "alternate"),
        default="automatic",
        help="Smear estimation mode (default: automatic).",
    )
    p_cal.add_argument(
        "--cadence-match",
        choices=("cadenceno", "positional"),
        default="cadenceno",
        help="Collateral cadence pairing mode (default: cadenceno).",
    )
    p_cal.add_argument(
        "--missing-cadence",
        choices=("nan", "error"),
        default="nan",
        help="Policy when TPF cadences lack collateral CADENCENO (default: nan).",
    )
    p_cal.add_argument(
        "--desmear-cutout",
        action="store_true",
        help="Run collateral-free FFI cutout desmearing (TessCut mode).",
    )
    p_cal.add_argument("--cache-dir", type=str, default=None, help="Cache directory path.")
    p_cal.add_argument("--ddir", type=str, default=None, help="Local reference directory.")
    p_cal.add_argument("--offline", action="store_true", help="Forbid network access.")
    p_cal.add_argument(
        "--strict-hpc",
        action="store_true",
        help="Raise ValueError if --cache-dir resolves under /fred.",
    )
    p_cal.add_argument("--overwrite", action="store_true", help="Overwrite existing output.")
    p_cal.add_argument("-v", "--verbose", action="store_true", help="Print progress.")

    return parser


def _run_prefetch(args: argparse.Namespace) -> int:
    jobs: dict[tuple[int, int, int], tuple[set[str], set[CadenceKind]]] = {}

    if args.tpf:
        for tpf_str in args.tpf:
            tpf_path = Path(tpf_str).expanduser()
            with fits.open(tpf_path, mode="readonly", memmap=True) as hdul:
                hdr0 = hdul[0].header
                hdr1 = hdul[1].header
                hdr2 = hdul[2].header
                sec = int(hdr0["SECTOR"])
                cam = int(hdr0["CAMERA"])
                ccd = int(hdr0["CCD"])
                col0 = int(hdr2["CRVAL1P"])
                nx = int(hdr2["NAXIS1"])
                outs = {output_for_rawx_col(col0 + i) for i in range(nx)}
                cad = infer_cadence_kind(
                    tpf_path=tpf_path, hdu0_header=hdr0, hdu1_header=hdr1
                )
            key = (sec, cam, ccd)
            if key not in jobs:
                jobs[key] = (set(), set())
            jobs[key][0].update(outs)
            jobs[key][1].add(cad)
    else:
        if args.sector is None:
            raise ValueError("Provide --sector (with --camera/--ccd or --task-id) or --tpf.")
        if args.task_id is not None:
            if not (0 <= args.task_id <= 15):
                raise ValueError(f"--task-id must be in 0..15, got {args.task_id}.")
            cam = (args.task_id // 4) + 1
            ccd = (args.task_id % 4) + 1
        elif args.camera is not None and args.ccd is not None:
            cam = args.camera
            ccd = args.ccd
        else:
            raise ValueError("Provide either --camera and --ccd, or --task-id 0..15.")
        jobs[(args.sector, cam, ccd)] = (set(OUTPUT_ORDER), {"2min"})

    explicit_outs = _parse_outputs(args.outputs) if args.outputs is not None else None
    explicit_cads = _parse_cadences(args.cadence) if args.cadence is not None else None

    for (sec, cam, ccd), (outs_set, cads_set) in sorted(jobs.items()):
        outs_tuple = (
            explicit_outs
            if explicit_outs is not None
            else tuple(o for o in OUTPUT_ORDER if o in outs_set)
        )
        cads_tuple: tuple[CadenceKind, ...] = (
            explicit_cads
            if explicit_cads is not None
            else tuple(c for c in ("2min", "20s") if c in cads_set)  # type: ignore[assignment]
        )
        published = prefetch_collateral(
            sector=sec,
            camera=cam,
            ccd=ccd,
            outputs=outs_tuple,
            cadences=cads_tuple,
            cache_dir=args.cache_dir,
            ddir=args.ddir,
            stage_dir=args.stage_dir,
            chunk_rows=args.chunk_rows,
            offline=args.offline,
            strict_hpc=args.strict_hpc,
            cleanup_staged_raw=not args.keep_raw,
            verbose=args.verbose,
        )
        if args.verbose:
            print(
                f"Prefetched {len(published)} artefacts for s{sec:04d} cam={cam} ccd={ccd}."
            )
    return 0


def _run_calibrate(args: argparse.Namespace) -> int:
    if args.desmear_cutout:
        cal_tpf = desmear_calibrated_cutout(args.tpf, verbose=args.verbose)
    else:
        cal_tpf = recalibrate_tpf(
            args.tpf,
            cache_dir=args.cache_dir,
            ddir=args.ddir,
            smear=args.smear,
            cadence_match=args.cadence_match,
            missing_cadence=args.missing_cadence,
            offline=args.offline,
            strict_hpc=args.strict_hpc,
            verbose=args.verbose,
        )
    out_path = Path(args.output).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cal_tpf.to_fits(str(out_path), overwrite=args.overwrite)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``mundey`` CLI."""
    parser = _build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.command == "prefetch":
        return _run_prefetch(args)
    if args.command == "calibrate":
        return _run_calibrate(args)
    raise ValueError(f"Unhandled subcommand {args.command!r}.")


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
