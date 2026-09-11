"""Aperture photometry functions, adapted from code by Zafar Rustamkulov."""

# Standard library imports
import pickle
from typing import Tuple

# Related third party imports
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.units import Unit
from astropy.wcs import WCS
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm

# Local
from spherex_source_discovery_tool.sdt_utils import get_exp_id


# Constants 
# Default pixel solid angle, in arcsec^2 per pixel, if no SA map exists.
DEFAULT_SA_ARCSEC2_PER_PIX = 6.15 * 6.15  # arcsec^2 / pix

# Convert arcsec^2 to steradians.
ARCSEC2_TO_SR = (np.pi / (180.0 * 3600.0))**2  # sr / arcsec^2

__all__ = ["interp2d_rg", "fetch_sapm", "sa_mjy_scale_map", "interp2d_like",
    "extract_photometry", "extract_dqf", "get_AB", "build_phot_table", "grab_star"]

# Functions
def interp2d_rg(x: np.ndarray, y: np.ndarray, z: np.ndarray, kind="linear", fill_value=np.nan) -> callable:
    """Creates a 2D interpolating function from a 2D array of values.

    Parameters
    ----------
    x : `numpy.ndarray`
        1D array of x pixel indices.
    y : `numpy.ndarray`
        1D array of y pixel indices.
    z : `numpy.ndarray`
        2D array of values at each (y, x) index.
    kind : `str`, optional
        Interpolation method, either "linear" or "nearest". Default is
        "linear".
    fill_value : `float`, optional
        Value to return for out-of-bounds points. Default is NaN.

    Returns
    -------
    f : `callable`
        f(x_new, y_new) that evaluates z[y, x] with the specified
        interpolation method.
    """
    # Initialize interpolator
    if np.any(np.diff(x) < 0):
        x = x[::-1]  # reverse x to be monotonically increasing
        z = z[:, ::-1]
    if np.any(np.diff(y) < 0):
        y = y[::-1]  # reverse y to be monotonically increasing
        z = z[::-1, :]
    interp = RegularGridInterpolator(
        (y, x), z, method=("linear" if kind == "linear" else "nearest"),
        bounds_error=False, fill_value=fill_value
    )

    # Make wrapper for interpolator to handle scalar inputs
    def f(x_new, y_new):
        Xn, Yn = np.meshgrid(np.atleast_1d(x_new), np.atleast_1d(y_new))
        pts = np.column_stack([Yn.ravel(), Xn.ravel()])
        return interp(pts).reshape(Yn.shape)
    return f

        
def fetch_sapm(sapm_s3: list[str]) -> Tuple[dict, dict]:
    """Fetch solid angle pixel maps and generate 2D interpolating functions.

    Parameters
    ----------
    sapm_s3 : `list` of `str`
        List of s3 paths for solid angle pixel maps, in order of
        detector.

    Returns
    -------
    SA_MAPS : `dict`
        Dictionary of integer keys (detector numbers) mapping to
        numpy array values (solid angle pixel maps).
    SA_MAPS_interp : `dict`
        Dictionary of integer keys (detector numbers) mapping to
        callable values (2D interpolating functions).
    """
    SA_MAPS = {}
    SA_MAPS_interp = {}
    for i, fpath in enumerate(sapm_s3):
        with fits.open(f"s3://{fpath}", fsspec_kwargs={"anon": True}) as hm:
            sapm = hm["IMAGE"].data
            SA_MAPS[i+1] = sapm
            ny, nx = sapm.shape
            x = np.arange(nx)
            y = np.arange(ny)
            interp = interp2d_rg(x, y, sapm, kind='linear')
            SA_MAPS_interp[i + 1] = interp

    return SA_MAPS, SA_MAPS_interp
    

def sa_mjy_scale_map(band: int, ny: int, nx: int, SA_MAPS: dict | None = None,
                     SA_MAPS_interp: dict | None = None,
                     default_sa_arcsec2_per_pix: float = DEFAULT_SA_ARCSEC2_PER_PIX) -> np.ndarray:
    """
    Build a per-pixel scale map converting IMAGE units [MJy/sr] to [mJy/pix].

    Parameters
    ----------
    band : `int`
        Detector number.
    ny, nx : `int`
        Number of pixels along Y/X axis.
    SA_MAPS : `dict`
        Dictionary of integer keys (detector numbers) mapping to
        numpy array values (solid angle pixel maps).
    SA_MAPS_interp : `dict`
        Dictionary of integer keys (detector numbers) mapping to
        callable values (2D interpolating functions).
    default_sa_arcsec2_per_pix : `float`
        Default pixel solid angle, in arcsec^2 per pixel, if no
        solid angle pixel map exists.

    Returns
    -------
    scale_map : `numpy.ndarray`
        Float array of shape (ny, nx). Multiply IMAGE [MJy/sr] by
        this to get [mJy/pix].

    Notes
    -----
    - If solid angle (SA) maps exist and contain per-pixel solid angle in 
      arcsec^2/pix, then:
          scale_map[pix] = SA_arcsec2_per_pix[pix] * ARCSEC2_TO_SR * 1e6 * 1e3
      The factors are:
          ARCSEC2_TO_SR : sr / arcsec^2
          1e6           : Jy / MJy
          1e3           : mJy / Jy
    - If SA maps are missing, this uses a constant solid angle:
          SA_arcsec2_per_pix = (6.15 arcsec)^2
    """
    use_fallback = True
    sa_arcsec2 = None

    if SA_MAPS is not None and (band in SA_MAPS):
        m = np.asarray(SA_MAPS[band], dtype=float)
        if m.shape == (ny, nx):
            sa_arcsec2 = m
            use_fallback = False
        elif SA_MAPS_interp is not None and (band in SA_MAPS_interp):
            # SA_MAPS_interp[band] should accept (x, y) arrays and return
            # (ny, nx).
            x = np.arange(nx, dtype=float)
            y = np.arange(ny, dtype=float)
            sa_arcsec2 = np.asarray(SA_MAPS_interp[band](x, y), dtype=float)
            if sa_arcsec2.shape == (ny, nx):
                use_fallback = False

    if use_fallback:
        sa_arcsec2 = np.full((ny, nx), float(default_sa_arcsec2_per_pix), dtype=float)

    scale_map = sa_arcsec2 * ARCSEC2_TO_SR * 1e6 * 1e3  # [mJy/pix] per [MJy/sr]
    return np.asarray(scale_map, dtype=float)


def interp2d_like(grid: np.ndarray) -> callable:
    """Creates a 2D interpolating function from a 2D grid of values.

    Parameters
    ----------
    grid : `numpy.ndarray`
        2D array of values.

    Returns
    -------
    f : `callable`
        f(x, y) that evaluates grid[y, x] with linear interpolation.
    """
    # Initialize interpolator
    ny, nx = grid.shape  # shape is (ny, nx)
    y_ax = np.arange(ny)
    x_ax = np.arange(nx)
    rgi = RegularGridInterpolator((y_ax, x_ax), grid,
                                  method="linear", bounds_error=False, fill_value=np.nan)

    # Make wrapper for interpolator to handle scalar inputs
    def f(x_new, y_new):
        x_arr = np.atleast_1d(x_new)
        y_arr = np.atleast_1d(y_new)
        pts = np.column_stack([y_arr, x_arr])  # order is (y, x)
        vals = rgi(pts)
        if np.isscalar(x_new) and np.isscalar(y_new):
            return float(vals[0])
        return vals

    return f


def extract_photometry(image: np.ndarray, positions: list[tuple], aperture_radii: list[float | int],
                       annulus_rad_in: float, annulus_rad_out: float, mask: np.ndarray = None) -> Table:
    """Extract aperture photometry for input positions and aperture radii
    within an image.

    Parameters
    ----------
    image : `numpy.ndarray`
        2D array of image data.
    positions : `list` of `tuple`
        List of (x, y) pixel positions for each source to photometer.
    aperture_radii : `list` of `float` or `int`
        List of aperture radii (in pixels).
    annulus_rad_in : `float`
        Inner radius of the background annulus (in pixels).
    annulus_rad_out : `float`
        Outer radius of the background annulus (in pixels).
    mask : `numpy.ndarray`, optional
        Boolean mask array where True indicates masked pixels, by default None.

    Returns
    -------
    out : `astropy.table.Table`
        Table containing aperture photometry results.

    Note
    ----
    This function returns sums in the SAME UNITS as `image`. If `image` is in
    surface brightness units (MJy/sr), the output sums are NOT fluxes. If
    `image` has been converted to (mJy/pix), the output sums are (mJy).
    """
    # Build bad data mask
    bad = ~np.isfinite(image)
    if mask is not None:
        bad |= mask

    # Initialize background annuli
    ann = CircularAnnulus(positions, r_in=annulus_rad_in, r_out=annulus_rad_out)
    ann_masks = ann.to_mask(method='center')

    # Calculate background stats for each aperture mask
    bkg_med = np.zeros(len(ann_masks), dtype=np.float32)
    bkg_mean = np.zeros(len(ann_masks), dtype=np.float32)
    bkg_std = np.zeros(len(ann_masks), dtype=np.float32)
    for i, m in enumerate(ann_masks):
        cut = m.multiply(image)
        vals = cut[m.data > 0]  # get data within annulus
        vals = vals[np.isfinite(vals)]  # only use finite values
        mean, med, std = sigma_clipped_stats(vals, sigma=3, maxiters=5)
        bkg_mean[i] = mean
        bkg_med[i] = med
        bkg_std[i] = std

    # Create output photometry table for this image
    out = Table()
    out['annulus_median'] = bkg_med
    out['annulus_mean'] = bkg_mean
    out['annulus_std'] = bkg_std

    # Perform aperture photometry for each desired radius
    for r in aperture_radii:
        ap = CircularAperture(positions, r=r)
        phot = aperture_photometry(image, ap, mask=bad, method='subpixel', subpixels=5)
        key = f'aperture_sum_r{r:.1f}'
        out[key] = phot['aperture_sum']
        out[f'aper_bkg_r{r:.1f}'] = bkg_med * ap.area
        out[f'aper_sum_bkgsub_r{r:.1f}'] = out[key] - out[f'aper_bkg_r{r:.1f}']
    return out


def extract_dqf(flag_mask: np.ndarray, positions: list[tuple], aperture_radii: list[float | int]) -> Table:
    """Extract data quality flag sums for input positions and aperture radii
    within an image.

    Parameters
    ----------
    flag_mask : `numpy.ndarray`
        2D boolean array where True indicates flagged pixels.
    positions : `list` of `tuple`
        List of (x, y) pixel positions for each source to photometer.
    aperture_radii : `list` of `float` or `int`
        List of aperture radii (in pixels).

    Returns
    -------
    out_table : `astropy.table.Table`
        Table containing data quality flag sums for each aperture radius.
    """
    out_table = Table()
    for r in aperture_radii:
        aperture = CircularAperture(positions, r=r)
        phot = aperture_photometry(flag_mask, aperture, method='center')
        colname = f'aperture_sum_r{r:.1f}'
        out_table[colname] = phot['aperture_sum']
    return out_table


def get_AB(f_jy: float | np.ndarray) -> float | np.ndarray:
    """Convert spectral flux density in Janskys to AB magnitude.

    Parameters
    ----------
    f_jy : `float` or `numpy.ndarray`
        Spectral flux density in Janskys.

    Returns
    -------
    ab : `float` or `numpy.ndarray`
        AB magnitude.
    """
    f = np.asarray(f_jy, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ab = -2.5 * np.log10(f / 3631.0)  # zero-point is 3631 Jy
    if f.ndim == 0:
        return float(ab) if np.isfinite(f) and f > 0 else np.nan  # handle scalars
    ab[(~np.isfinite(f)) | (f <= 0)] = np.nan
    return ab


def build_phot_table(results_table: Table, matched_stars: pd.DataFrame, aperture_radii: list[float],
                     phot_path: str, sapm_s3: list[str]):
    """Perform aperture photometry for all matched stars in all cutouts.

    Parameters
    ----------
    results_table : `astropy.table.Table`
        Table containing cutout results. Must contain 'hdus' and 'spectral_wcs'
        columns.
    matched_stars : `pandas.DataFrame`
        Dataframe of selected stars to photometer. Must contain 'ra' and 'dec'
        `float` columns (in degrees).
    aperture_radii : `list` of `float`
        List of aperture radii (in pixels) for photometry extraction.
    phot_path : `str`
        Path to save the photometry results (pickle file).
    sapm_s3 : `list` of `str`
        List of s3 paths for solid angle pixel maps, in order of
        detector.

    Note
    ----
    - If IMAGE is [MJy/sr], then we convert to [mJy/pix] before aperture sums:
          image_mjy_pix = image_MJySr * scale_map
      where scale_map converts [MJy/sr] to [mJy/pix].
    - Variance is propagated consistently:
          var_mjy2_pix2 = var_(MJy/sr)^2 * scale_map^2
      so summed variances are in [mJy^2].
    """
    # Initialize index map for HDU access
    # NOTE: We use an index map because the results table does not save
    #       HDULists, and the EXTNAMEs are unique to each cutout index.
    hdu_idx_map = {
        "IMAGE": 0,
        "FLAGS": 1,
        "VARIANCE": 2,
        "ZODI": 3,
        "PSF": 4,
        "WCS-WAVE": 5
    }

    # Fetch solid angle pixel maps and interpolating functions
    SA_MAPS, SA_MAPS_interp = fetch_sapm(sapm_s3)

    # Iterate over cutouts
    for i in tqdm(range(1, len(results_table) + 1),
                  total=len(results_table),
                  desc="Processing cutouts for photometry extraction"):
        # Fetch HDUs and spectral WCS
        row = results_table[i-1]
        hdus = row["hdus"]
        spectral_wcs = row["spectral_wcs"]

        # Don't attempt to perform photometry if HDUs or WCS are missing
        if hdus is None:
            print(f"No HDUs found for cutout {i} with URI: {row['uri']}. Skipping this cutout.")
            continue
        if spectral_wcs is None:
            print(f"No spectral WCS found for cutout {i} with URI: {row['uri']}. Skipping this cutout.")
            continue

        # Fetch image, variance, and flags data
        image_data = hdus[hdu_idx_map["IMAGE"]].data.copy()
        image_header = hdus[hdu_idx_map["IMAGE"]].header
        flux_unit = Unit(image_header["BUNIT"])
        var_img = hdus[hdu_idx_map["VARIANCE"]].data
        flags = hdus[hdu_idx_map["FLAGS"]].data
        flags_header = hdus[hdu_idx_map["FLAGS"]].header
        flag_dict = {key[3:]: flags_header[key] for key in flags_header if key.startswith('MP_')}
        flag_dict = dict(sorted(flag_dict.items(), key=lambda x: x[1]))
        flag_names = list(flag_dict.keys())
        flag_masks = {name: (flags & (1 << flag_dict[name])) > 0 for name in flag_names}

        # Fetch metadata and spatial WCS
        current_image_id = get_exp_id(hdus[hdu_idx_map["IMAGE"]])
        current_mjd = image_header['MJD-OBS']
        wcs = WCS(image_header)

        # Generate per-pixel wavelength and bandpass grids
        ny, nx = image_data.shape
        y_grid, x_grid = np.mgrid[0:ny, 0:nx]
        wavelength_grid, bandpass_grid = spectral_wcs.pixel_to_world(x_grid, y_grid)

        # Make interpolating functions for values at arbitrary pixel positions
        wavelength_interp = interp2d_like(np.asarray(wavelength_grid, dtype=float))
        bandpass_interp = interp2d_like(np.asarray(bandpass_grid, dtype=float))

        # Prepare data frame for storing photometry results
        if 'obs_log' not in matched_stars.columns:
            matched_stars['obs_log'] = [[] for _ in range(len(matched_stars))]

        # Convert sky to pixel positions and filter to sources w/in image
        coords = SkyCoord(
            ra=matched_stars['ra'],
            dec=matched_stars['dec'],
            unit=u.deg,
            frame='icrs'
        )
        x_pix, y_pix = wcs.world_to_pixel(coords)
        ny, nx = image_data.shape
        valid = (x_pix > -0.5) & (x_pix < nx) & (y_pix > -0.5) & (y_pix < ny)
        valid_idx = np.where(valid)[0]
        valid_xy = [(x_pix[i], y_pix[i]) for i in valid_idx]

        if len(valid_xy) == 0:
            continue

        # Get band from image header, if it exists. If not, default to 0,
        # which will trigger a constant SA map fallback in sa_mjy_scale_map.
        band = int(image_header.get("DETECTOR", 0))

        if flux_unit == u.MJy / u.sr:
            scale_map = sa_mjy_scale_map(
                band=band,
                ny=ny,
                nx=nx,
                SA_MAPS=SA_MAPS,
                SA_MAPS_interp=SA_MAPS_interp,
                default_sa_arcsec2_per_pix=DEFAULT_SA_ARCSEC2_PER_PIX,
            )
            image_for_phot = image_data * scale_map          # mJy/pix
            var_for_phot = var_img * (scale_map ** 2)        # (mJy/pix)^2
            flux_out_unit = "mJy"
            var_out_unit = "mJy^2"
        else:
            # If inputs are already flux-like, we do not guess conversions.
            # Photometry sums will be in the same units as the input image.
            image_for_phot = image_data
            var_for_phot = var_img
            flux_out_unit = str(flux_unit)
            var_out_unit = f"({flux_unit})^2"

        # Extract photometry
        P = extract_photometry(image_for_phot, valid_xy, aperture_radii, annulus_rad_in=7, annulus_rad_out=16)

        # AB mags only make sense if we are truly in flux density units.
        # In the MJy/sr branch, P[aper_sum_bkgsub] is mJy, so convert to Jy.
        for r in aperture_radii:
            r_str = f"{r:.1f}"
            flux_key = f'aper_sum_bkgsub_r{r_str}'
            abmag_key = f'abmag_r{r_str}'
            if flux_key in P.colnames:
                if flux_unit == u.MJy / u.sr:
                    P[abmag_key] = get_AB(1e-3 * P[flux_key])  # mJy -> Jy
                else:
                    P[abmag_key] = np.full(len(P[flux_key]), np.nan, dtype=float)

        # Grab flags
        dqf_tables = {}
        for flag in flag_names:
            mask = flag_masks[flag]  # 2D bool mask for that DQF
            dqf_tables[flag] = extract_dqf(mask, valid_xy, aperture_radii)

        # Grab variance
        var_bad = ~np.isfinite(var_for_phot)
        var_tables = {}
        for r in (aperture_radii):
            aperture = CircularAperture(valid_xy, r=r)
            phot = aperture_photometry(var_for_phot, aperture, mask=var_bad, method='subpixel', subpixels=5)
            var_tables[f"{r:.1f}"] = phot['aperture_sum']

        # Add photometry results for each source within image
        for i_local, i_global in enumerate(valid_idx):
            x, y = valid_xy[i_local]
            obs = {
                'wavelength': float(wavelength_interp(x, y)),
                'bandwidth': float(bandpass_interp(x, y)),
                'aperture_radii': aperture_radii,
                'aper_sum_bkgsub': [],
                'variance_sum': [],
                'abmag': [],
                'dqf_sum': {flag: [] for flag in flag_names},
                'image_id': current_image_id,
                'mjd_center': current_mjd
            }

            for r in aperture_radii:
                r_str = f"{r:.1f}"
                flux_key = f'aper_sum_bkgsub_r{r_str}'
                abmag_key = f'abmag_r{r_str}'

                obs['aper_sum_bkgsub'].append(float(P[flux_key][i_local]))
                obs['variance_sum'].append(float(var_tables[r_str][i_local]))
                obs['abmag'].append(float(P[abmag_key][i_local]) if abmag_key in P.colnames else np.nan)

                for flag in flag_names:
                    dqf_val = dqf_tables[flag][f'aperture_sum_r{r_str}'][i_local]
                    obs['dqf_sum'][flag].append(float(dqf_val))

            matched_stars.at[i_global, 'obs_log'].append(obs)

        with open(phot_path, 'wb') as f:
            pickle.dump(matched_stars, f)


def grab_star(obs_log_entry: list[dict], ap_rad: float, bad_flags: list[str]) -> tuple:
    """For each measurement of the requested aperture radius, filters out
    non-finite flux measurements and counts bad photometry flags.

    Parameters
    ----------
    obs_log_entry : `list` of `dict`
        List of observation log entries for a single star.
    ap_rad : `float`
        Aperture radius to consider.
    bad_flags : `list` of `str`
        List of bad photometry flags to check for.


    Returns
    -------
    `tuple` of `numpy.ndarray`
        Wavelengths (in microns), bandwidths (in microns), AB magnitudes,
        absolute AB magnitude errors, fluxes (in mJy), variances
        (in mJy^2), and bad flag counts for the specified aperture radius.
    """
    # Initialize lists to hold extracted data
    wavelengths = []
    abmags = []
    ab_errs = []
    bad_flags_aperture = []
    fluxes = []
    variances = []
    bandwidths = []

    # Iterate over measurements for this star
    for obs in obs_log_entry:
        if 'wavelength' in obs and 'abmag' in obs:
            # Fetch general values
            wavelengths.append(obs['wavelength'])
            bandwidths.append(obs['bandwidth'])

            # Fetch aperture-specific values
            idx_ap = obs['aperture_radii'].index(ap_rad)
            abmags.append(obs['abmag'][idx_ap])
            fluxes.append(obs['aper_sum_bkgsub'][idx_ap])
            variances.append(obs['variance_sum'][idx_ap])

            # Sum number of bad pixels in aperture
            bad_pixels = 0
            for flag in bad_flags:
                bad_pixels += obs['dqf_sum'][flag][idx_ap]
            bad_flags_aperture.append(bad_pixels)

    # Filter out non-finite flux measurements and convert to arrays
    valid = np.isfinite(fluxes)
    wavelengths = np.array(wavelengths)[valid]
    abmags = np.array(abmags)[valid]
    fluxes = np.array(fluxes)[valid]
    variances = np.array(variances)[valid]
    bandwidths = np.array(bandwidths)[valid]
    bad_flags_aperture = np.array(bad_flags_aperture)[valid]

    # Compute AB magnitude errors from fluxes and variances
    # AB magnitude error propagation:
    # sigma_mag = (2.5 / ln(10)) * (sigma_f / f)
    _AB_DMAG_FACTOR = 2.5 / np.log(10)  # 1.0857362047581296

    with np.errstate(divide='ignore', invalid='ignore'):
        ab_errs = _AB_DMAG_FACTOR * (np.sqrt(np.maximum(variances, 0.0)) / fluxes)
    ab_errs[~np.isfinite(ab_errs)] = np.nan

    return wavelengths, bandwidths, abmags, np.abs(ab_errs), fluxes, variances, bad_flags_aperture
