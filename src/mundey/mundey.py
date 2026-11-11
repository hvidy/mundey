"""Legacy ``mundey_tpf`` calibration class modernized for Python 3.11+ and NumPy 2.x."""

from __future__ import annotations

import importlib.resources
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import Element

import astropy.coordinates as coord
import astropy.units as u
import lightkurve
import numpy as np
from astropy.io import fits

from ._contract import write_mundey_provenance

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable: Any, **_kwargs: Any) -> Any:
        return iterable


_RANK_WARNING: type[Warning] = getattr(
    getattr(np, "exceptions", np), "RankWarning", UserWarning
)


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


def _build_blf_time_cadence_columns(
    hdu1: fits.BinTableHDU,
) -> tuple[fits.Column, fits.Column]:
    """Construct ``TIME`` and ``CADENCENO`` (``'J'``) columns for ``-blf.fits`` caches.

    SPOC TPFs store ``TIMECORR`` (``'E'``, 32-bit float) in ``TTYPE2``/``TFORM2``
    and ``CADENCENO`` (``'J'``, 32-bit integer) in ``TTYPE3``/``TFORM3``.
    """
    hdr = hdu1.header
    col_time = fits.Column(
        name=hdr.get("TTYPE1", "TIME"),
        format=hdr.get("TFORM1", "D"),
        unit=hdr.get("TUNIT1", "BJD - 2457000, days"),
        disp=hdr.get("TDISP1", "D14.7"),
        array=hdu1.data["TIME"],
    )
    cad_name = hdr.get("TTYPE3", "CADENCENO")
    cad_form = hdr.get("TFORM3", "J")
    cad_disp = hdr.get("TDISP3", "I10")
    if cad_name != "CADENCENO" or not str(cad_form).strip().endswith(("J", "K", "I")):
        cad_name = "CADENCENO"
        cad_form = "J"
        cad_disp = "I10"
    col_cad = fits.Column(
        name=cad_name,
        format=cad_form,
        disp=cad_disp,
        array=np.asarray(hdu1.data["CADENCENO"], dtype=np.int32),
    )
    return col_time, col_cad


class mundey_tpf(lightkurve.TessTargetPixelFile):
    """TargetPixelFile subclass for TESS bright-star recalibration."""

    def calibrate(
        self,
        verbose: bool = True,
        smear: str = "automatic",
        ddir: str | Path = "./",
        target: str | None = None,
        use_simbad: bool | None = None,
        offline: bool = False,
    ) -> None:
        """Perform a recalibration of a TESS target pixel file from ``RAW_CNTS``.

        Parameters
        ----------
        verbose : bool, default True
            Print progress messages during calibration.
        smear : {'automatic', 'standard', 'alternate'}, default 'automatic'
            Smear estimation mode. ``'automatic'`` selects ``'alternate'`` when
            the cutout reaches row 2048 and ``'standard'`` otherwise.
        ddir : str or Path, default './'
            Directory containing collateral and static CCD reference files.
        target : str, optional
            Explicit target label override.
        use_simbad : bool, optional
            Whether to query SIMBAD for target identification. Defaults to
            ``False`` whenever ``OBJECT`` is present in ``hdu[0].header`` or
            ``offline=True``.
        offline : bool, default False
            Disable all network queries (e.g. on HPC compute nodes).
        """
        resolved_target = _resolve_target_name(
            self.hdu[0].header,
            self.ra,
            self.dec,
            self.targetid,
            target=target,
            use_simbad=use_simbad,
            offline=offline,
        )
        ddir_path = Path(ddir)
        self.ddir = f"{ddir_path}/"

        if verbose:
            print(f"Beginning calibration of {resolved_target}")

        flux = np.copy(self.hdu[1].data["RAW_CNTS"])
        cols = np.arange(self.column, self.column + self.shape[2])

        outputs: list[str] = []
        for col in cols:
            if 45 <= col <= 556:
                outputs.append("A")
            elif 557 <= col <= 1068:
                outputs.append("B")
            elif 1069 <= col <= 1580:
                outputs.append("C")
            elif 1581 <= col <= 2092:
                outputs.append("D")

        unique_outputs = [out for out in ("A", "B", "C", "D") if out in outputs]
        path_name = Path(self.path).name
        parts = path_name.split("-")

        x1 = self.hdu[2].header["CRVAL1P"] - 1
        dx = self.hdu[2].header["NAXIS1"]
        y1 = self.hdu[2].header["CRVAL2P"] - 1
        dy = self.hdu[2].header["NAXIS2"]

        if smear == "automatic":
            smear = "alternate" if (y1 + dy >= 2048) else "standard"

        if verbose:
            print(f"Sector: {self.sector}")
            print(f"Camera: {self.camera}")
            print(f"CCD: {self.ccd}")
            print(f"Outputs: {unique_outputs}")
            print("")
            print("Searching for collateral target pixel files")

        n_out = len(unique_outputs)
        sm_cal = np.zeros((n_out, flux.shape[0], 10, 512))
        tv_cal = np.zeros((n_out, flux.shape[0], 2078, 11))
        calcsmear = np.zeros(n_out, dtype=bool)
        calcblack = np.zeros(n_out, dtype=bool)

        for idx, output in enumerate(unique_outputs):
            stem_smrow = (
                f"{parts[0]}-{parts[1]}-smrow-{self.camera}-{self.ccd}-"
                f"{output.lower()}-{parts[3]}"
            )
            smearfile = ddir_path / f"{stem_smrow}-s_col-blf.fits"
            try:
                if verbose:
                    print(f"... loading smear file {smearfile}")
                with fits.open(smearfile) as smrow:
                    sm_cal[idx, :, 0, :] = smrow[1].data["SMEAR"]
                calcsmear[idx] = False
            except FileNotFoundError:
                smrowfile = ddir_path / f"{stem_smrow}-s_col.fits"
                if verbose:
                    print(f"... file not found: {smearfile}")
                    print("... smear correction will be calculated from raw counts")
                    print(f"... loading smear row file {smrowfile}")
                    print("")
                with fits.open(smrowfile) as smrow:
                    sm_cal[idx] = smrow[1].data["SMROW_RAW"]
                calcsmear[idx] = True

            stem_tvcol = (
                f"{parts[0]}-{parts[1]}-tvcol-{self.camera}-{self.ccd}-"
                f"{output.lower()}-{parts[3]}"
            )
            onedblackfile = ddir_path / f"{stem_tvcol}-s_col-blf.fits"
            try:
                if verbose:
                    print(f"... loading 1D black file {onedblackfile}")
                with fits.open(onedblackfile) as tvcol:
                    tv_cal[idx, :, :, 0] = tvcol[1].data["1DBLACK"]
                calcblack[idx] = False
            except FileNotFoundError:
                tvcolfile = ddir_path / f"{stem_tvcol}-s_col.fits"
                if verbose:
                    print(f"... file not found: {onedblackfile}")
                    print("... 1D black correction will be calculated from raw counts")
                    print(f"... loading trailing virtual column file {tvcolfile}")
                    print("")
                with fits.open(tvcolfile) as tvcol:
                    tv_cal[idx] = tvcol[1].data["TVCOL_RAW"]
                calcblack[idx] = True

        # Subtract fixed offset, add mean black
        for idx, output in enumerate(unique_outputs):
            indices = [i for i, x in enumerate(outputs) if x == output]
            flux[:, :, indices] = self.fixedoffset_meanblack(
                flux[:, :, indices], output, verbose=verbose
            )
            if calcsmear[idx]:
                sm_cal[idx] = self.fixedoffset_meanblack(sm_cal[idx], output, verbose=False)
            if calcblack[idx]:
                tv_cal[idx] = self.fixedoffset_meanblack(tv_cal[idx], output, verbose=False)

        # Remove 2D black
        black2dfile = ddir_path / f"tess2018323-{self.camera}-{self.ccd}-2dblack.fits"
        if verbose:
            print("")
            print("Loading 2D black model")
            print(f"... {black2dfile}")

        with fits.open(black2dfile) as black2dfits:
            black2dimg = np.array(black2dfits[1].data, copy=True)

        flux = self.twodblack(flux, black2dimg[y1 : y1 + dy, x1 : x1 + dx], verbose=verbose)

        output_bounds = {
            "A": (44, 556, 2092, 2103),
            "B": (556, 1068, 2103, 2114),
            "C": (1068, 1580, 2114, 2125),
            "D": (1580, 2092, 2125, 2136),
        }
        for idx, output in enumerate(unique_outputs):
            sx1, sx2, tx1, tx2 = output_bounds[output]
            if calcsmear[idx]:
                sm_cal[idx] = self.twodblack(
                    sm_cal[idx], black2dimg[2058:2068, sx1:sx2], verbose=False
                )
            if calcblack[idx]:
                tv_cal[idx] = self.twodblack(tv_cal[idx], black2dimg[:, tx1:tx2], verbose=False)

        # Calculate 1D black
        for idx, output in enumerate(unique_outputs):
            if calcblack[idx]:
                tv_cal[idx] = self.calconedblack(tv_cal[idx], output, verbose=verbose)

        # Remove 1D black
        tv_cal_1d = tv_cal[:, :, :, 0]
        if n_out == 1:
            flux = self.onedblack(flux, tv_cal_1d[0, :, y1 : y1 + dy], verbose=verbose)
            if calcsmear[0]:
                sm_cal[0] = self.onedblack(
                    sm_cal[0], tv_cal_1d[0, :, 2058:2068], verbose=verbose
                )
        else:
            flux = flux.astype(float, copy=False)
            for idx, output in enumerate(unique_outputs):
                indices = [i for i, x in enumerate(outputs) if x == output]
                flux[:, :, indices] = self.onedblack(
                    flux[:, :, indices], tv_cal_1d[idx, :, y1 : y1 + dy], verbose=verbose
                )
                if calcsmear[idx]:
                    sm_cal[idx] = self.onedblack(
                        sm_cal[idx], tv_cal_1d[idx, :, 2058:2068], verbose=verbose
                    )

        # Correct for non-linearity and gain
        for idx, output in enumerate(unique_outputs):
            indices = [i for i, x in enumerate(outputs) if x == output]
            flux[:, :, indices] = self.linearity_gain(
                flux[:, :, indices], output, verbose=verbose
            )
            if calcsmear[idx]:
                sm_cal[idx] = self.linearity_gain(sm_cal[idx], output, verbose=False)

        # Correct for LDE undershoot
        for idx, output in enumerate(unique_outputs):
            indices = [i for i, x in enumerate(outputs) if x == output]
            flux[:, :, indices] = self.undershoot(
                flux[:, :, indices], output, verbose=verbose
            )
            if calcsmear[idx]:
                sm_cal[idx] = self.undershoot(sm_cal[idx], output, verbose=False)

        # Calculate smear correction
        for idx, output in enumerate(unique_outputs):
            if calcsmear[idx]:
                sm_cal[idx] = self.calcsmear(
                    flux, sm_cal[idx], output, smear=smear, verbose=verbose
                )

        # Remove smear correction
        sm_cal_2d = sm_cal[:, :, 0, :]
        flux = self.smear(flux, sm_cal_2d, unique_outputs[0], verbose=verbose)

        # Flatfield correction
        flux = self.flatfield(flux, verbose=verbose)

        # Convert to electrons per second
        flux = self.to_electrons_per_second(flux, verbose=verbose)

        self.hdu[1].data["FLUX"] = flux
        col_names = self.hdu[1].data.names or ()
        if "FLUX_BKG" in col_names:
            self.hdu[1].data["FLUX_BKG"] = np.float32(0.0)
        if "FLUX_BKG_ERR" in col_names:
            self.hdu[1].data["FLUX_BKG_ERR"] = np.float32(0.0)
        write_mundey_provenance(self, bkgsub=False)

        if verbose:
            print("")
            print("Calibration Complete")

    def fixedoffset_meanblack(
        self,
        img: np.ndarray,
        output: str,
        verbose: bool = True,
    ) -> np.ndarray:
        """Subtract fixed electronic offset and add output mean black."""
        if verbose:
            print("")
            print(f"Adjusting for fixed offset and mean black - Output {output}")

        img = img - self.hdu[1].header["FXDOFF"]
        img = img + self.hdu[1].header[f"MEANBLC{output}"]
        return img

    def twodblack(
        self,
        img: np.ndarray,
        twodblackimg: np.ndarray,
        verbose: bool = True,
    ) -> np.ndarray:
        """Subtract static 2D black image."""
        if verbose:
            print("Performing 2-D black correction")
        return img - twodblackimg

    def calconedblack(
        self,
        tvcol: np.ndarray,
        output: str,
        verbose: bool = True,
    ) -> np.ndarray:
        """Fit per-cadence 1D black polynomial to virtual columns and cache ``-blf.fits``."""
        if verbose:
            print("")
            print("Calculating 1-D black correction")
            print("... this step may take a while;")
            print("... if you haven't already done so")
            print("... you might like to take this opportunity")
            print("... to join your union.")
            _maybe_display_union_link()

        nrows = tvcol.shape[1]
        x = np.arange(nrows) / nrows - 0.5
        ks = np.arange(20, 50)
        black1d = np.zeros((tvcol.shape[0], tvcol.shape[1]))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", _RANK_WARNING)
            for frame in tqdm(np.arange(tvcol.shape[0]), desc="Frame", disable=not verbose):
                aics = []
                y = np.nanmean(tvcol[frame, :, :-1], axis=1)

                for k in ks:
                    p = np.polynomial.Polynomial.fit(x, y, int(k))
                    aic = 2 * k + nrows * np.log(np.sum((y - p(x)) ** 2))
                    aics.append(aic)

                order = int(ks[np.argmin(np.asarray(aics))])
                p = np.polynomial.Polynomial.fit(x, y, order)
                black1d[frame, :] = p(x)

        col1, col2 = _build_blf_time_cadence_columns(self.hdu[1])
        col3 = fits.Column(
            name="1DBLACK",
            format=f"{black1d.shape[1]}D",
            unit="count",
            disp="F14.7",
            dim=f"({black1d.shape[1]})",
            array=black1d,
        )

        blackhdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1, col2, col3]))
        parts = Path(self.path).name.split("-")
        onedblackfile = Path(self.ddir) / (
            f"{parts[0]}-{parts[1]}-tvcol-{self.camera}-{self.ccd}-"
            f"{output.lower()}-{parts[3]}-s_col-blf.fits"
        )
        blackhdu.writeto(onedblackfile, overwrite=True)

        tvcol[:, :, 0] = black1d
        return tvcol

    def onedblack(
        self,
        img: np.ndarray,
        onedblackimg: np.ndarray,
        verbose: bool = True,
    ) -> np.ndarray:
        """Subtract 1D black profile across columns."""
        if verbose:
            print("Performing 1-D black correction")
        return img - np.expand_dims(onedblackimg, axis=2)

    def linearity_gain(
        self,
        img: np.ndarray,
        output: str,
        verbose: bool = True,
    ) -> np.ndarray:
        """Apply TESS nonlinearity polynomial and output gain."""
        if verbose:
            print("")
            print("Correcting for non-linearity and gain")

        factor = 1.0
        if bool(self.hdu[0].header["CRMITEN"]):
            factor = (self.hdu[0].header["CRBLKSZ"] - 2) / self.hdu[0].header["CRBLKSZ"]

        img = img / (factor * self.hdu[1].header["NUM_FRM"])

        lintree = ET.parse(Path(self.ddir) / "tess2018143203310-41006_100-linearity.xml")
        linroot = lintree.getroot()
        xpath = (
            f".//*[@cameraNumber='{self.camera}']"
            f"[@ccdNumber='{self.ccd}']"
            f"[@ccdOutput='{output}']/linearityPoly"
        )
        poly = linroot.find(xpath)
        if poly is None:
            raise KeyError(
                f"Linearity polynomial not found for camera={self.camera}, "
                f"ccd={self.ccd}, output={output}."
            )

        originx = float(poly.attrib["originx"])
        scalex = float(poly.attrib["scalex"])
        coeffs = [float(elem.attrib["value"]) for elem in poly.findall("coeffs/coeff")]
        p = np.polynomial.Polynomial(np.asarray(coeffs))

        return img * p((img - originx) * scalex) * self.hdu[1].header[f"GAIN{output}"]

    def undershoot(
        self,
        img: np.ndarray,
        output: str,
        verbose: bool = True,
    ) -> np.ndarray:
        """Correct for LDE serial readout undershoot using bundled ``tess-undershoot.xml``."""
        if verbose:
            print("")
            print("Correcting for LDE undershoot")

        underroot = load_undershoot_xml()
        xpath = (
            f".//*[@cameraNumber='{self.camera}']"
            f"[@ccdNumber='{self.ccd}']"
            f"[@ccdOutput='{output}']"
        )
        node = underroot.find(xpath)
        if node is None:
            raise KeyError(
                f"Undershoot factor not found for camera={self.camera}, "
                f"ccd={self.ccd}, output={output}."
            )
        undershoot_coeff = float(node.attrib["undershoot"])

        if output in ("A", "C"):
            img[:, :, 1:] = img[:, :, 1:] + undershoot_coeff * img[:, :, :-1]
            img[:, :, 0] = img[:, :, 0] * (1 + undershoot_coeff)
        elif output in ("B", "D"):
            img[:, :, :-1] = img[:, :, :-1] + undershoot_coeff * img[:, :, 1:]
            img[:, :, -1] = img[:, :, -1] * (1 + undershoot_coeff)

        return img

    def calcsmear(
        self,
        img: np.ndarray,
        smrow: np.ndarray,
        output: str,
        smear: str = "alternate",
        verbose: bool = True,
    ) -> np.ndarray:
        """Estimate collateral or TPF-anchored smear and cache ``-blf.fits``."""
        if verbose:
            print("")
            if smear == "alternate":
                print("Calculating photometric smear from TPF background")
            if smear == "standard":
                print("Calculating photometric smear from collateral data")

        sx1_map = {"A": 44, "B": 556, "C": 1068, "D": 1580}
        sx1 = sx1_map[output]

        x1 = self.hdu[2].header["CRVAL1P"] - 1 - sx1
        dx = self.hdu[2].header["NAXIS1"]
        dy = self.hdu[2].header["NAXIS2"]

        smcor = np.nanmedian(smrow, axis=1)

        if smear == "alternate":
            medimg = np.nanmedian(img, axis=0)
            maxrows = [
                int(np.argmax(medimg[:, col][medimg[:, col] < 1e5])) for col in np.arange(dx)
            ]
            maxrow = int(np.nanmin(maxrows)) - 10

            smmask = np.zeros((dy, dx), dtype=bool)
            for col in np.arange(dx):
                idx = medimg[:maxrow, col] < np.percentile(medimg[:maxrow, col], 20)
                smmask[:maxrow, col][idx] = True

            smbkgd = np.nanmean(
                img.T[smmask.T].reshape(dx, int(np.sum(smmask) / dx), img.shape[0]),
                axis=1,
            ).T

            msk = (np.arange(smrow.shape[2]) < x1) | (np.arange(smrow.shape[2]) > x1 + dx)
            bkgd = np.nanmin(smbkgd, axis=1) - np.nanmedian(smcor[:, msk], axis=1)

            smcor[:, np.max([x1, 0]) : x1 + dx] = smbkgd[
                :, -dx - np.min([x1, 0]) :
            ] - np.tile(bkgd[:, np.newaxis], dx + np.min([x1, 0]))

        col1, col2 = _build_blf_time_cadence_columns(self.hdu[1])
        col3 = fits.Column(
            name="SMEAR",
            format=f"{smcor.shape[1]}D",
            unit="count",
            disp="F14.7",
            dim=f"({smcor.shape[1]})",
            array=smcor,
        )

        smearhdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1, col2, col3]))
        parts = Path(self.path).name.split("-")
        smearfile = Path(self.ddir) / (
            f"{parts[0]}-{parts[1]}-smrow-{self.camera}-{self.ccd}-"
            f"{output.lower()}-{parts[3]}-s_col-blf.fits"
        )
        smearhdu.writeto(smearfile, overwrite=True)

        smrow[:, 0, :] = smcor
        return smrow

    def smear(
        self,
        img: np.ndarray,
        sm: np.ndarray,
        output: str,
        verbose: bool = True,
    ) -> np.ndarray:
        """Subtract 1D column smear profile across science rows."""
        if verbose:
            print("")
            print("Removing photometric smear")

        sx1_map = {"A": 44, "B": 556, "C": 1068, "D": 1580}
        sx1 = sx1_map[output]

        x1 = self.hdu[2].header["CRVAL1P"] - 1 - sx1
        dx = self.hdu[2].header["NAXIS1"]

        if sm.shape[0] == 2:
            newsm = np.zeros((sm.shape[1], 1024), dtype=sm.dtype)
            newsm[:, :512] = sm[0]
            newsm[:, 512:] = sm[1]
            sm_combined = newsm
        else:
            sm_combined = sm[0]

        return img - np.expand_dims(sm_combined[:, x1 : x1 + dx], axis=1)

    def flatfield(self, img: np.ndarray, verbose: bool = True) -> np.ndarray:
        """Divide by static TESS flat-field model."""
        flatfile = Path(self.ddir) / f"tess2018323-{self.camera}-{self.ccd}-flat.fits"
        if verbose:
            print("")
            print("Loading flat field")
            print(f"... {flatfile}")

        x1 = self.hdu[2].header["CRVAL1P"] - 1
        dx = self.hdu[2].header["NAXIS1"]
        y1 = self.hdu[2].header["CRVAL2P"] - 1
        dy = self.hdu[2].header["NAXIS2"]

        with fits.open(flatfile) as flatfits:
            if verbose:
                print("Performing flatfield correction")
            flat = np.array(flatfits[1].data[y1 : y1 + dy, x1 : x1 + dx], copy=True)

        return img / flat

    def to_electrons_per_second(self, img: np.ndarray, verbose: bool = True) -> np.ndarray:
        """Convert calibrated electrons per exposure to electrons per second."""
        if verbose:
            print("")
            print("Converting to electrons per second")
        return img / self.hdu[1].header["INT_TIME"]
