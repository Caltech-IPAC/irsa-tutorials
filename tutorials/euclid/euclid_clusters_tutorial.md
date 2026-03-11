---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Euclid Galaxy Clusters Analysis Tutorial

This tutorial explores galaxy clusters in the Euclid Q1 Multi-Epoch Release (MER) data to demonstrate cluster detection and validation techniques.
We select a cluster from this paper (https://arxiv.org/abs/2503.19196), identify a control field that is covered by Euclid Q1 and at least 15 arcmin from any known clusters.
We download multi-band images and galaxy catalogs, apply clustering algorithms to confirm the existence of galaxy overdensities and identify cluster members.
We analyze color-magnitude diagrams, extract spectra, and cross-match with external databases for validation.
This approach allows us to compare cluster and field galaxy properties and assess the reliability of cluster detections.

```{code-cell} ipython3
# Install required packages if needed
#!pip install astroquery requests aiohttp scikit-learn

# Import all necessary libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import shutil
import requests
import json
import s3fs
from tqdm import tqdm

# Astropy imports
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
from astropy.utils.data import download_file
from astropy.table import QTable

# External data access
from astroquery.ipac.irsa import Irsa
from astroquery.ipac.ned import Ned
from astroquery.exceptions import RemoteServiceError
from pyvo.dal import DALQueryError
from requests.exceptions import Timeout, ConnectionError

# Machine learning
from sklearn.cluster import DBSCAN

# Statistics and signal processing
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d

# Set up plotting style
plt.style.use('default')
sns.set_palette("husl")

Number = u.def_unit("Number")
u.add_enabled_units([Number])
```

## 1. Loading the Cluster Catalog
We begin by loading the Euclid Q1 cluster catalog. The catalog contains 35 galaxy clusters with photometric redshifts, coordinates, and richness estimates from the PZWav algorithm.

```{code-cell} ipython3
# Load the Euclid Q1 cluster catalog (https://arxiv.org/abs/2503.19196)
url = "https://ar5iv.labs.arxiv.org/html/2503.19196"
fname = "euclid_q1_clusters.csv"

# Read all tables from the arXiv HTML rendering
dfs = pd.read_html(url)

# Select the table containing the cluster catalog
df = dfs[1].copy()

# Drop the spurious row that contains units instead of data
df = df[df["ID"].notna()]

# Normalize Unicode minus signs and extract numeric RA/Dec values
# (HTML tables sometimes concatenate values or include formatting artifacts)
for col in ["RAPZWav", "DecPZWav"]:
    df[col] = (
        df[col].astype(str)
               .str.replace("−", "-", regex=False)
               .str.extract(r"([-+]?\d+(?:\.\d+)?)")[0]
               .astype(float)
    )

# Rename LaTeX-style column names to clean, code-friendly names
df = df.rename(columns={
    "zPZWavz_{\\mathrm{PZWav}}": "zPZWav",
    "zAMICOz_{\\mathrm{AMICO}}": "zAMICO",
    "λPmem\\lambda_{\\mathrm{Pmem}}": "lambdaPmem",
})

# Overwrite (or create) the cleaned CSV for downstream use
df.to_csv(fname, index=False)

print(f"Dataset shape: {df.shape}")
df.head(3)
```

## 2. A Cluster and Control Field Selection

We implement a systematic field selection process: (1) randomly select a target cluster from the catalog, and (2) identify a control field that maintains a minimum 20 arcmin separation from all known cluster locations.
This ensures the control field represents the general field population without contamination from known structures.

```{code-cell} ipython3
# Set random seed for reproducibility
np.random.seed(45)

# Function to check if a field has exactly 4 MER images (one per band)
def check_mer_tile_requirement(coord, search_radius=2.0):
    """Check whether a sky coordinate is covered by exactly 4 Euclid MER science images.

    The Euclid MER mosaic is organized into tiles. A position inside a single
    tile returns exactly 4 images, one per photometric band (VIS, Y, J, H).
    Positions near tile boundaries may overlap two tiles and return more than 4
    images; those are rejected to avoid partial-coverage complications.

    Parameters
    ----------
    coord : `~astropy.coordinates.SkyCoord`
        Sky coordinate to check.
    search_radius : float, optional
        Radius in arcminutes within which to search for MER images. Default is 2.0.

    Returns
    -------
    valid : bool
        True if exactly 4 science images are found; False otherwise.
    mer_images : `~astropy.table.Table` or None
        Table of matching MER image records, or None if the query failed.
    """
    try:
        mer_images = Irsa.query_sia(pos=(coord, search_radius * u.arcmin), collection='euclid_DpdMerBksMosaic')
        mer_images = mer_images[(mer_images['facility_name'] == 'Euclid') & (mer_images['dataproduct_subtype'] == 'science')]

        if len(mer_images) == 4:
            print(f"  ✓ Found exactly 4 images for coordinates")
            return True, mer_images
        else:
            print(f"  ✗ Found {len(mer_images)} images (expected 4) - trying next location")
            return False, mer_images
    except (DALQueryError, Timeout, ConnectionError, OSError) as e:
        print(f"  ✗ Error querying MER tiles: {e}")
        return False, None

# Function to find a random control field that avoids all cluster locations
def find_control_field_corrected(cluster_df, cluster_ra, cluster_dec, min_distance_arcmin=15, max_attempts=100):
    """Find a random control field offset from the known cluster catalog.

    Generates candidate control fields by randomly offsetting 25–35 arcmin from
    a randomly chosen catalog cluster in a random direction. Each candidate is
    accepted only if it lies at least ``min_distance_arcmin`` away from every
    cluster in the catalog, ensuring the control field is a genuine blank field
    representative of the general galaxy population.

    RA offsets are corrected for the cosine of the declination so that angular
    separations are accurate. If no valid field is found within ``max_attempts``,
    a fallback position directly offset from the target cluster is returned.

    Parameters
    ----------
    cluster_df : `pandas.DataFrame`
        Cluster catalog with columns ``RAPZWav`` and ``DecPZWav`` (degrees).
    cluster_ra : float
        RA of the target cluster in degrees (used for the fallback only).
    cluster_dec : float
        Dec of the target cluster in degrees (used for the fallback only).
    min_distance_arcmin : float, optional
        Minimum angular separation in arcminutes from all catalog clusters.
        Default is 15.
    max_attempts : int, optional
        Maximum number of random candidates to try before falling back.
        Default is 100.

    Returns
    -------
    control_ra : float
        RA of the selected control field in degrees.
    control_dec : float
        Dec of the selected control field in degrees.
    """
    for attempt in range(max_attempts):
        anchor_cluster = cluster_df.sample(n=1).iloc[0]
        anchor_ra, anchor_dec = anchor_cluster['RAPZWav'], anchor_cluster['DecPZWav']

        offset_arcmin = np.random.uniform(25, 35)
        offset_deg = offset_arcmin / 60.0
        direction = np.random.uniform(0, 2 * np.pi)

        ra_offset = offset_deg * np.cos(direction) / np.cos(np.radians(anchor_dec))
        dec_offset = offset_deg * np.sin(direction)

        control_ra = (anchor_ra + ra_offset) % 360.0
        control_dec = np.clip(anchor_dec + dec_offset, -90.0, 90.0)

        control_coord = SkyCoord(ra=control_ra, dec=control_dec, unit='deg')
        cluster_coords = SkyCoord(ra=cluster_df['RAPZWav'], dec=cluster_df['DecPZWav'], unit='deg')
        distances_arcmin = control_coord.separation(cluster_coords).to(u.arcmin).value

        if np.all(distances_arcmin > min_distance_arcmin):
            return control_ra, control_dec

    # Fallback
    fallback_offset = (min_distance_arcmin + 10) / 60.0
    fallback_ra = (cluster_ra + fallback_offset) % 360.0
    fallback_dec = np.clip(cluster_dec + fallback_offset, -90.0, 90.0)
    return fallback_ra, fallback_dec

# Function to find a valid cluster and control field combination
def find_valid_cluster_control_pair(cluster_df, max_attempts=50):
    """Find a cluster and matched control field each covered by a single MER tile.

    Iterates over randomly selected clusters and, for each, attempts to find a
    control field via :func:`find_control_field_corrected` such that both
    positions are covered by exactly 4 MER science images — i.e., neither falls
    across a tile boundary. The single-tile requirement keeps the data download
    straightforward and avoids mosaicking artefacts at tile edges.

    If no valid pair is found after ``max_attempts`` iterations, falls back to
    the first catalog entry with a best-effort control field.

    Parameters
    ----------
    cluster_df : `pandas.DataFrame`
        Cluster catalog with columns ``RAPZWav``, ``DecPZWav``, ``ID``,
        and ``zPZWav``.
    max_attempts : int, optional
        Maximum number of cluster/control combinations to try. Default is 50.

    Returns
    -------
    cluster : `pandas.Series`
        Row from ``cluster_df`` for the selected cluster.
    control_ra : float
        RA of the matched control field in degrees.
    control_dec : float
        Dec of the matched control field in degrees.
    cluster_mer_images : `~astropy.table.Table`
        MER image table for the cluster field.
    control_mer_images : `~astropy.table.Table`
        MER image table for the control field.
    """
    print("Searching for cluster and control field combination with single MER tiles...")

    for attempt in range(max_attempts):
        print(f"\nAttempt {attempt + 1}/{max_attempts}:")

        cluster = cluster_df.sample(n=1).iloc[0]
        cluster_coord = SkyCoord(ra=cluster['RAPZWav'], dec=cluster['DecPZWav'], unit='deg')
        print(f"Selected cluster: ID={cluster['ID']}, z={cluster['zPZWav']:.2f}")

        cluster_single_tile, cluster_mer_images = check_mer_tile_requirement(cluster_coord)
        if not cluster_single_tile:
            print("  Skipping cluster - requires multiple tiles")
            continue

        control_ra, control_dec = find_control_field_corrected(cluster_df, cluster['RAPZWav'], cluster['DecPZWav'])
        control_coord = SkyCoord(ra=control_ra, dec=control_dec, unit='deg')

        control_single_tile, control_mer_images = check_mer_tile_requirement(control_coord)
        if not control_single_tile:
            print("  Skipping control field - requires multiple tiles")
            continue

        print("  ✓ Both cluster and control field require single tiles!")
        return cluster, control_ra, control_dec, cluster_mer_images, control_mer_images

    # Fallback
    print(f"\nWarning: Could not find valid cluster/control pair after {max_attempts} attempts.")
    cluster = cluster_df.iloc[0]
    cluster_coord = SkyCoord(ra=cluster['RAPZWav'], dec=cluster['DecPZWav'], unit='deg')
    control_ra, control_dec = find_control_field_corrected(cluster_df, cluster['RAPZWav'], cluster['DecPZWav'])
    control_coord = SkyCoord(ra=control_ra, dec=control_dec, unit='deg')

    _, cluster_mer_images = check_mer_tile_requirement(cluster_coord)
    _, control_mer_images = check_mer_tile_requirement(control_coord)

    return cluster, control_ra, control_dec, cluster_mer_images, control_mer_images

# Find valid cluster and control field combination
cluster, control_ra, control_dec, cluster_mer_images, control_mer_images = find_valid_cluster_control_pair(df)

# Create coordinate objects for later use
cluster_coord = SkyCoord(ra=cluster['RAPZWav'], dec=cluster['DecPZWav'], unit='deg')
control_coord = SkyCoord(ra=control_ra, dec=control_dec, unit='deg')
```

## 3. Data Download and Caching

Download and cache Euclid Q1 MER mosaics for both the selected cluster and control field with 12 arcmin cutouts (time consuming).
The MER tile querying was already performed in Section 2 to ensure single-tile requirements.

```{code-cell} ipython3
# Define parameters for cutouts
im_cutout = 12.0 * u.arcmin  # 12 arcminutes for both cluster and control fields

# Create cache directory
cache_dir = 'data'
os.makedirs(cache_dir, exist_ok=True)

# use S3 cloud_access to cache a cutout FITS instead of downloading the full tile
s3 = s3fs.S3FileSystem(
    anon=True,
    default_block_size=256 * 1024 * 1024,  # fewer network roundtrips
    default_cache_type="readahead",
    default_fill_cache=True,
)

def download_and_cache_field(mer_images, field_name, field_coord, field_id):
    """Stream cutout FITS images from S3 and cache them locally.

    Uses the ``cloud_access`` column in the MER image table to open each band's
    FITS file directly from the AWS S3 mirror without downloading the full tile.
    A ``Cutout2D`` region of size ``im_cutout`` centred on ``field_coord`` is
    written to a local cache file. On subsequent calls, cached files are reused.

    After caching, the function re-opens each FITS, applies a second
    ``Cutout2D`` to confirm the footprint, and returns the 2-D pixel arrays
    together with the VIS-band WCS for downstream analysis.

    Parameters
    ----------
    mer_images : `~astropy.table.Table`
        MER image table as returned by :func:`check_mer_tile_requirement`,
        with columns ``energy_bandpassname``, ``access_url``, and
        ``cloud_access``.
    field_name : str
        Human-readable label (e.g. ``"cluster"`` or ``"control"``) used in
        progress messages.
    field_coord : `~astropy.coordinates.SkyCoord`
        Centre of the cutout.
    field_id : str
        Unique identifier used to name the cached FITS files on disk.

    Returns
    -------
    cutouts : dict
        Dictionary mapping band name (``'VIS'``, ``'Y'``, ``'J'``, ``'H'``)
        to the corresponding 2-D NumPy array.
    cutout_wcs : `~astropy.wcs.WCS`
        WCS object derived from the VIS-band cutout.
    """
    print(f"\nProcessing {field_name} field...")

    # Get URLs for each band (unchanged)
    vis_url = mer_images[mer_images['energy_bandpassname'] == 'VIS'][0]['access_url']
    y_url = mer_images[mer_images['energy_bandpassname'] == 'Y'][0]['access_url']
    j_url = mer_images[mer_images['energy_bandpassname'] == 'J'][0]['access_url']
    h_url = mer_images[mer_images['energy_bandpassname'] == 'H'][0]['access_url']

    # Download and cache images (ONLY this part is changed)
    cached_files = {}
    for band, url in [('VIS', vis_url), ('Y', y_url), ('J', j_url), ('H', h_url)]:
        cache_file = os.path.join(cache_dir, f'{band}_{field_id}.fits')
        if not os.path.exists(cache_file):
            print(f"  Downloading {band} band...")

            # Use cloud_access to avoid downloading the full FITS
            row = mer_images[mer_images['energy_bandpassname'] == band][0]
            cloud = json.loads(row['cloud_access'])
            s3_obj = f"{cloud['aws']['bucket_name']}/{cloud['aws']['key']}"

            # Open remote FITS and write only the cutout region to cache_file
            with s3.open(s3_obj, "rb") as f:
                with fits.open(f, memmap=False, lazy_load_hdus=True) as hdul:
                    hdu0 = hdul[0]
                    w0 = WCS(hdu0.header)
                    cut0 = Cutout2D(hdu0.section, position=field_coord, size=im_cutout, wcs=w0)
                    fits.writeto(cache_file, cut0.data, header=cut0.wcs.to_header(), overwrite=True)

            cached_files[band] = cache_file
        else:
            print(f"  Using cached {band} band")
            cached_files[band] = cache_file

    # Create cutouts and store WCS (unchanged)
    cutouts = {}
    cutout_wcs = None
    for band in ['VIS', 'Y', 'J', 'H']:
        hdu = fits.open(cached_files[band])
        cutout = Cutout2D(hdu[0].data, position=field_coord, size=im_cutout, wcs=WCS(hdu[0].header))
        cutouts[band] = cutout.data
        if band == 'VIS':  # Store WCS from VIS band
            cutout_wcs = cutout.wcs
        hdu.close()

    return cutouts, cutout_wcs

# Download and cache both fields
cluster_cutouts, cluster_cutout_wcs = download_and_cache_field(
    cluster_mer_images, "cluster", cluster_coord, cluster["ID"]
)

control_cutouts, control_cutout_wcs = download_and_cache_field(
    control_mer_images, "control", control_coord, f"CONTROL_{cluster['ID']}"
)

print(f"\nBoth fields processed successfully!")
print(f"Cluster field cutout size: {cluster_cutouts['VIS'].shape}")
print(f"Control field cutout size: {control_cutouts['VIS'].shape}")
```

## 4. Multi-band Image Visualization

Create and display VIS, Y, J, H bands plus RGB composite for both cluster and control fields.

```{code-cell} ipython3
# Improved normalization for consistent stretching between fields
def normalize_with_consistent_stretch(cluster_cutouts, control_cutouts, lower_percentile=1, upper_percentile=99):
    """Normalize cluster and control cutouts using a shared percentile stretch.
    The RGB composite is assembled as H→R, J→G, VIS→B, which maps redder
    (older) stellar populations toward red in the image.

    Parameters
    ----------
    cluster_cutouts : dict
        Band-keyed dictionary of 2-D arrays for the cluster field.
    control_cutouts : dict
        Band-keyed dictionary of 2-D arrays for the control field.
    lower_percentile : float, optional
        Lower percentile used for ``vmin``. Default is 1.
    upper_percentile : float, optional
        Upper percentile used for ``vmax``. Default is 99.

    Returns
    -------
    norm_cluster : dict
        Normalised band arrays for the cluster field.
    norm_control : dict
        Normalised band arrays for the control field.
    cluster_rgb : `~numpy.ndarray`
        (H, W, 3) RGB array for the cluster field.
    control_rgb : `~numpy.ndarray`
        (H, W, 3) RGB array for the control field.
    """
    norm_cluster = {}
    norm_control = {}

    for band in ['VIS', 'Y', 'J', 'H']:
        # Combine both fields to calculate global percentiles
        combined_data = np.concatenate([cluster_cutouts[band].flatten(), control_cutouts[band].flatten()])

        # Calculate global percentiles
        vmin = np.percentile(combined_data, lower_percentile)
        vmax = np.percentile(combined_data, upper_percentile)

        # Apply same normalization to both fields
        norm_cluster[band] = np.clip((cluster_cutouts[band] - vmin) / (vmax - vmin), 0, 1)
        norm_control[band] = np.clip((control_cutouts[band] - vmin) / (vmax - vmin), 0, 1)

        print(f"{band} band: vmin={vmin:.2e}, vmax={vmax:.2e}")

    # Create RGB composites
    cluster_rgb = np.dstack([norm_cluster['H'], norm_cluster['J'], norm_cluster['VIS']])
    control_rgb = np.dstack([norm_control['H'], norm_control['J'], norm_control['VIS']])

    return norm_cluster, norm_control, cluster_rgb, control_rgb
```

```{code-cell} ipython3
# Process both fields with consistent normalization

print("Using consistent stretching between cluster and control fields...")
cluster_norm_cutouts, control_norm_cutouts, cluster_rgb, control_rgb = normalize_with_consistent_stretch(
    cluster_cutouts, control_cutouts, lower_percentile=1, upper_percentile=99
)

# Plot both fields side by side with consistent stretching
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
bands = ['VIS', 'Y', 'J', 'H']
titles = ['VIS', 'Y', 'J', 'H', 'RGB']

# Cluster field (top row)
for i, (band, title) in enumerate(zip(bands, titles)):
    axes[0, i].imshow(cluster_norm_cutouts[band], cmap='gray', origin='lower')
    axes[0, i].set_title(f'Cluster - {title}')
    axes[0, i].axis('off')

axes[0, 4].imshow(cluster_rgb, origin='lower')
axes[0, 4].set_title('Cluster - RGB')
axes[0, 4].axis('off')

# Control field (bottom row)
for i, (band, title) in enumerate(zip(bands, titles)):
    axes[1, i].imshow(control_norm_cutouts[band], cmap='gray', origin='lower')
    axes[1, i].set_title(f'Control - {title}')
    axes[1, i].axis('off')

axes[1, 4].imshow(control_rgb, origin='lower')
axes[1, 4].set_title('Control - RGB')
axes[1, 4].axis('off')

plt.tight_layout()
plt.show()
```

**Figure 1. Euclid Q1 MER cutouts of the cluster candidate and control field.**
Top row: the cluster candidate field shown in the Euclid VIS band and the three NISP near-infrared bands (Y, J, H), followed by an RGB composite constructed as **R = H, G = J, B = VIS**.
Bottom row: a nearby control field displayed in the same set of bands and with the same RGB mapping.

```{code-cell} ipython3
# Query galaxies in both fields with BOX search
table_mer = 'euclid_q1_mer_catalogue'
table_phz = 'euclid_q1_phz_photo_z'

# Convert cutout size to degrees
cutout_deg = im_cutout.to(u.deg).value

# Function to query galaxies for a field
def query_galaxies_for_field(ra, dec, field_name, redshift_center, redshift_width=0.1):
    """Query galaxies within a redshift slice around a sky position.

    Submits a TAP/ADQL query to IRSA joining the Euclid Q1 MER photometric
    catalogue and the photo-z catalogue. Only objects satisfying all of the
    following criteria are returned:

    - Positive flux in all four bands (VIS, Y, J, H).
    - Photo-z classification = 2 (galaxy).
    - Fractional 90% credible interval width < 0.20, i.e.
      ``(phz_90_int2 - phz_90_int1) / (1 + phz_median) < 0.20``, selecting
      sources with well-constrained photometric redshifts.
    - ``phz_median`` within ``redshift_center ± redshift_width``.

    Parameters
    ----------
    ra : float
        Field centre RA in degrees.
    dec : float
        Field centre Dec in degrees.
    field_name : str
        Label used in progress messages.
    redshift_center : float
        Central photometric redshift of the cluster.
    redshift_width : float, optional
        Half-width of the redshift slice. Default is 0.1.

    Returns
    -------
    result : `~astropy.table.Table`
        Catalogue of galaxies matching all quality and spatial cuts.
    """
    print(f"\nQuerying galaxies for {field_name} field...")

    adql = (f"SELECT DISTINCT mer.object_id, mer.ra, mer.dec, "
            f"phz.flux_vis_unif, phz.flux_y_unif, phz.flux_j_unif, phz.flux_h_unif, "
            f"phz.phz_classification, phz.phz_median, phz.phz_90_int1, phz.phz_90_int2 "
            f"FROM {table_mer} AS mer "
            f"JOIN {table_phz} as phz "
            f"ON mer.object_id = phz.object_id "
            f"WHERE 1 = CONTAINS(POINT('ICRS', mer.ra, mer.dec), "
            f"BOX('ICRS', {ra}, {dec}, {cutout_deg/np.cos(np.radians(dec))}, {cutout_deg})) "
            f"AND phz.flux_vis_unif > 0 "
            f"AND phz.flux_y_unif > 0 "
            f"AND phz.flux_j_unif > 0 "
            f"AND phz.flux_h_unif > 0 "
            f"AND phz.phz_classification = 2 "
            f"AND ((phz.phz_90_int2 - phz.phz_90_int1) / (1 + phz.phz_median)) < 0.20 "
            f"AND phz.phz_median BETWEEN {redshift_center-redshift_width} AND {redshift_center+redshift_width}")

    result = Irsa.query_tap(adql).to_table()
    print(f"Found {len(result)} galaxies in {field_name} field (z = {redshift_center:.2f} ± {redshift_width:.2f})")

    return result

# Query galaxies for both fields in the cluster redshift slice
cluster_galaxies = query_galaxies_for_field(
    cluster['RAPZWav'], cluster['DecPZWav'], "cluster",
    cluster['zPZWav'], redshift_width=0.12
)

control_galaxies = query_galaxies_for_field(
    control_ra, control_dec, "control",
    cluster['zPZWav'], redshift_width=0.12
)

# Convert to pandas DataFrames for easier analysis
cluster_df_galaxies = cluster_galaxies.to_pandas()
control_df_galaxies = control_galaxies.to_pandas()

# Calculate Overdensity delta
n_cluster = len(cluster_galaxies)
n_control = len(control_galaxies)

# Avoid division by zero just in case
if n_control > 0:
    overdensity = (n_cluster - n_control) / n_control
    ratio = n_cluster / n_control
else:
    overdensity = np.inf
    ratio = np.inf

print(f"\nDensity Analysis:")
print(f"Cluster field: {n_cluster} galaxies")
print(f"Control field: {n_control} galaxies")
print(f"Density Ratio (N_clust/N_cont): {ratio:.2f}")
print(f"Measured Overdensity (delta): {overdensity:.2f}")

cluster_df_galaxies.head()
```

## 5. Cluster Finding Algorithm

We apply DBSCAN (Density-Based Spatial Clustering of Applications with Noise) to the galaxy catalogues to confirm the cluster detection and identify candidate member galaxies.
DBSCAN is well-suited to this problem because it finds arbitrarily shaped overdensities without requiring a prior on the number of clusters, and it naturally labels sparse galaxies as noise, useful for separating cluster members from the field population.

**How redshift enters the analysis**

Redshift is not an input to DBSCAN itself.
Instead, both the cluster and control catalogues were already filtered in Section 4 to a narrow photo-z slice centred on the target cluster redshift (Δz = ±0.12).
DBSCAN therefore operates on the projected 2-D sky distribution of galaxies.
A real cluster should produce a compact spatial overdensity in this slice; the control field should show mostly noise.

**Pixel coordinates vs. sky coordinates**

Clustering is performed in image pixel space rather than in sky (RA/Dec) coordinates.
The main practical reason is that `eps` can be directly interpreted as an angular scale without computing a great-circle distance matrix.
At the native Euclid VIS pixel scale of `~0.1 arcsec/pixel`, `eps = 500` pixels corresponds to `~50` arcsec on the sky, and at `z~0.4` this maps to roughly 300 kpc, comparable to the virial radius of a moderate-mass cluster.
This approach is accurate over the small field of view of a single 12-arcmin cutout where projection effects are negligible.

**Parameter choices and customisation**

The defaults `eps=500` and `min_samples=18` were tuned empirically for a 12-arcmin cutout at z~0.4 with Euclid Q1 data.
They are not universal, the right values depend on the cluster redshift, richness, cutout size, and the depth of the photo-z sample. In particular:

- The physical scale probed by `eps` changes with redshift.
  At higher redshift the same angular scale subtends a larger physical distance, so you may want to increase `eps` to keep the probed scale near the virial radius.
- `min_samples` sets the minimum richness threshold for a detection; increasing it suppresses spurious small groups, while decreasing it recovers lower-richness structures at the cost of more false positives.
- A more physically motivated approach is to convert the expected virial radius at the cluster redshift (in arcseconds) to pixels and use that as `eps`.

**Query/cutout mismatch and edge effects**

The galaxy catalogue is queried with a rectangular BOX in ICRS coordinates, while the image cutout is defined in pixel space.
Due to WCS projection distortions, a small fraction of galaxies returned by the query will map to pixel positions outside the image bounds and are excluded before clustering.
If the cluster or control field is close to a tile edge, this fraction may be larger, reducing the effective sample size and potentially biasing results.
The function prints the number of galaxies within bounds as a diagnostic, if a large fraction are lost, consider selecting a field better centred on the tile.

```{code-cell} ipython3
# Function to apply DBSCAN clustering with validity check (needed due to query/cutout mismatch)
def apply_dbscan_clustering(galaxy_df, wcs, rgb_image, field_name, eps=500, min_samples=18):
    """Apply DBSCAN to detect galaxy overdensities in a redshift-selected sample. DBSCAN operates on the 2-D projected
    spatial distribution of galaxies at approximately the same redshift.

    **Parameter choices**

    ``eps=500`` and ``min_samples=18`` were tuned empirically for a 12-arcmin
    cutout at z~0.4 with a photo-z slice width of Δz=0.12. These values are
    data-dependent — feel free to adjust them to suit your science case:

    - Increase ``eps`` for higher-redshift or lower-richness clusters where
      members are more widely separated in projection.
    - Increase ``min_samples`` to suppress spurious small groups; decrease it
      to detect lower-richness structures.
    - For a physically motivated approach, estimate the virial radius at the
      cluster redshift, convert to pixels, and use that as ``eps``.

    Parameters
    ----------
    galaxy_df : `pandas.DataFrame`
        Galaxies pre-filtered to the cluster redshift slice, with columns
        ``ra`` and ``dec`` in degrees.
    wcs : `~astropy.wcs.WCS`
        WCS of the image cutout, used to project RA/Dec to pixel coordinates.
    rgb_image : `~numpy.ndarray`
        The RGB image array; its shape defines the valid pixel bounds.
    field_name : str
        Label used in printed diagnostics.
    eps : float, optional
        Maximum distance in pixels between neighbours in DBSCAN. Default is
        500 px (~50 arcsec at VIS resolution, ~300 kpc at z~0.4).
    min_samples : int, optional
        Minimum number of galaxies required to form a core point. Default is 18.

    Returns
    -------
    labels : `~numpy.ndarray`
        DBSCAN cluster labels for the valid galaxy subset; -1 indicates noise.
    valid_galaxy_coords : `~numpy.ndarray`
        (N, 2) pixel coordinate array for the galaxies that were clustered.
    n_clusters : int
        Number of clusters found (noise points excluded).
    n_noise : int
        Number of noise (unassigned) galaxies.
    """
    print(f"\nApplying DBSCAN clustering to {field_name} field...")

    # Convert galaxy coordinates to pixel coordinates
    galaxy_pixels = wcs.world_to_pixel_values(galaxy_df['ra'], galaxy_df['dec'])

    # Filter galaxies to only show those within the image bounds (needed due to query/cutout mismatch)
    image_height, image_width = rgb_image.shape[:2]
    valid_mask = ((galaxy_pixels[0] >= 0) & (galaxy_pixels[0] < image_width) &
                  (galaxy_pixels[1] >= 0) & (galaxy_pixels[1] < image_height))

    valid_galaxy_coords = np.column_stack([galaxy_pixels[0][valid_mask], galaxy_pixels[1][valid_mask]])

    print(f"  Total galaxies: {len(galaxy_df)}")
    print(f"  Galaxies within image bounds: {valid_mask.sum()}")

    # Apply DBSCAN clustering only to valid galaxies
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(valid_galaxy_coords)
    labels = clustering.labels_

    # Count clusters (excluding noise points labeled as -1)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)

    print(f"{field_name} field: {n_clusters} clusters, {n_noise} noise points")

    return labels, valid_galaxy_coords, n_clusters, n_noise

# Apply clustering to both fields (with validity check)
cluster_labels, cluster_galaxy_coords, cluster_n_clusters, cluster_n_noise = apply_dbscan_clustering(
    cluster_df_galaxies, cluster_cutout_wcs, cluster_rgb, "Cluster"
)

control_labels, control_galaxy_coords, control_n_clusters, control_n_noise = apply_dbscan_clustering(
    control_df_galaxies, control_cutout_wcs, control_rgb, "Control"
)
```

```{code-cell} ipython3
# Plot clustering results for both fields in one figure (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8),
                               subplot_kw={'projection': cluster_cutout_wcs})

# Cluster field (left subplot)
ax1.imshow(cluster_rgb, origin='lower')
if len(cluster_labels) > 0:
    cluster_unique_labels = set(cluster_labels)
    cluster_colors = plt.cm.Spectral(np.linspace(0, 1, len(cluster_unique_labels)))

    for k, col in zip(cluster_unique_labels, cluster_colors):
        if k == -1:
            col = 'black'
            marker = ''
            size = 10
            alpha = 0.3
        else:
            marker = 'o'
            size = 30
            alpha = 0.7

        class_member_mask = (cluster_labels == k)
        xy = cluster_galaxy_coords[class_member_mask]
        if len(xy) > 0:  # Check if there are any points to plot
            ax1.scatter(xy[:, 0], xy[:, 1], c=[col], marker=marker, s=size, alpha=alpha)
else:
    print("No cluster data to plot")

ax1.set_xlabel('RA')
ax1.set_ylabel('Dec')
ax1.set_title(f'Cluster Field - {cluster_n_clusters} clusters')

# Control field (right subplot)
ax2.imshow(control_rgb, origin='lower')
if len(control_labels) > 0:
    control_unique_labels = set(control_labels)
    control_colors = plt.cm.Spectral(np.linspace(0, 1, len(control_unique_labels)))

    for k, col in zip(control_unique_labels, control_colors):
        if k == -1:
            col = 'black'
            marker = ''
            size = 10
            alpha = 0.3
        else:
            marker = 'o'
            size = 30
            alpha = 0.7

        class_member_mask = (control_labels == k)
        xy = control_galaxy_coords[class_member_mask]
        if len(xy) > 0:  # Check if there are any points to plot
            ax2.scatter(xy[:, 0], xy[:, 1], c=[col], marker=marker, s=size, alpha=alpha)
else:
    print("No control field data to plot")

ax2.set_xlabel('RA')
ax2.set_ylabel('Dec')
ax2.set_title(f'Control Field - {control_n_clusters} clusters')

plt.tight_layout()
plt.show()
```

**Figure 2. DBSCAN clustering of galaxy candidates in redshift slices: cluster field vs. control field.**
Each panel shows the RGB cutout (R = H, G = J, B = VIS) with points overplotted for sources selected within successive redshift slices and clustered using **DBSCAN** (density-based spatial clustering).
Points assigned to a DBSCAN cluster are shown as colored circular markers; points labeled as noise/outliers (DBSCAN label = −1) are shown with low-opacity markers.
The left panel (cluster field) contains multiple spatial overdensities identified by DBSCAN, consistent with the expectation that a real cluster field may include one or more galaxy concentrations within the scanned redshift range.
In the right panel (control field), no real clusters are found as expected. In some examples, a “cluster” can be identified around a bright star; this is interpreted as an **artifact-driven detection** (e.g., spurious sources near diffraction spikes/halos in Euclid Q1), rather than a genuine galaxy overdensity.

+++

## 6. Color-Magnitude Diagram Analysis

We analyze the color-magnitude properties of cluster and field galaxies to understand their stellar populations and star formation histories.
The Y-H color vs H magnitude diagram reveals differences in galaxy properties between cluster and field environments.

```{code-cell} ipython3
# Function to identify cluster members and field galaxies
def identify_cluster_members(galaxy_df, labels, galaxy_coords, field_name):
    """Separate a galaxy sample into cluster members and field galaxies.

    Re-applies the same image-bounds filter used in
    :func:`apply_dbscan_clustering` to align ``galaxy_df`` with the DBSCAN
    ``labels`` array, which only covers galaxies that fell within the image
    footprint. Galaxies with ``label != -1`` are cluster members; those with
    ``label == -1`` are noise (field galaxies).

    Parameters
    ----------
    galaxy_df : `pandas.DataFrame`
        Full galaxy catalogue for the field (before image-bounds filtering),
        with columns ``ra`` and ``dec`` in degrees.
    labels : `~numpy.ndarray`
        DBSCAN labels as returned by :func:`apply_dbscan_clustering`.
    galaxy_coords : `~numpy.ndarray`
        Pixel coordinates of the clustered galaxies (unused directly but kept
        for API consistency).
    field_name : str
        Either ``"Cluster"`` or ``"Control"``; selects the WCS and image used
        for the bounds filter.

    Returns
    -------
    cluster_members : `pandas.DataFrame`
        Galaxies assigned to a cluster group, with an added ``cluster_id``
        column containing the DBSCAN label.
    field_galaxies : `pandas.DataFrame`
        Galaxies classified as noise (``cluster_id = -1``).
    """
    print(f"\nAnalyzing {field_name} field galaxy populations...")

    # The labels array corresponds to the galaxies that were actually used in clustering
    # (those within image bounds). We need to create a mapping.

    # Get the WCS for coordinate conversion
    if field_name == "Cluster":
        wcs = cluster_cutout_wcs
        rgb_image = cluster_rgb
    else:
        wcs = control_cutout_wcs
        rgb_image = control_rgb

    # Convert galaxy coordinates to pixel coordinates
    galaxy_pixels = wcs.world_to_pixel_values(galaxy_df['ra'], galaxy_df['dec'])

    # Filter galaxies to only show those within the image bounds (same as in clustering)
    image_height, image_width = rgb_image.shape[:2]
    valid_mask = ((galaxy_pixels[0] >= 0) & (galaxy_pixels[0] < image_width) &
                  (galaxy_pixels[1] >= 0) & (galaxy_pixels[1] < image_height))

    # Get only the valid galaxies (those used in clustering)
    valid_galaxies = galaxy_df[valid_mask].copy()

    print(f"  Total galaxies in catalog: {len(galaxy_df)}")
    print(f"  Galaxies within image bounds: {len(valid_galaxies)}")
    print(f"  Labels array length: {len(labels)}")

    # Now the labels array should match the valid_galaxies DataFrame
    if len(valid_galaxies) != len(labels):
        print(f"  WARNING: Mismatch between valid galaxies ({len(valid_galaxies)}) and labels ({len(labels)})")
        # Use the minimum length to avoid errors
        min_len = min(len(valid_galaxies), len(labels))
        valid_galaxies = valid_galaxies.iloc[:min_len]
        labels = labels[:min_len]

    # Create a mask for cluster members (labels != -1)
    cluster_member_mask = labels != -1
    field_galaxy_mask = labels == -1

    # Get cluster members and field galaxies
    cluster_members = valid_galaxies[cluster_member_mask].copy()
    field_galaxies = valid_galaxies[field_galaxy_mask].copy()

    # Add cluster assignment information
    cluster_members['cluster_id'] = labels[cluster_member_mask]
    field_galaxies['cluster_id'] = -1  # Field galaxies

    print(f"  Cluster members: {len(cluster_members)} galaxies")
    print(f"  Field galaxies: {len(field_galaxies)} galaxies")

    # Count galaxies per cluster
    if len(cluster_members) > 0:
        cluster_counts = cluster_members['cluster_id'].value_counts().sort_index()

    return cluster_members, field_galaxies

# Analyze both fields
cluster_members_cluster_field, field_galaxies_cluster_field = identify_cluster_members(
    cluster_df_galaxies, cluster_labels, cluster_galaxy_coords, "Cluster"
)

cluster_members_control_field, field_galaxies_control_field = identify_cluster_members(
    control_df_galaxies, control_labels, control_galaxy_coords, "Control"
)

# Combine all field galaxies (from both cluster and control fields)
all_field_galaxies = pd.concat([field_galaxies_cluster_field, field_galaxies_control_field], ignore_index=True)

# Combine all cluster members
all_cluster_members = pd.concat([cluster_members_cluster_field], ignore_index=True)

print(f"\nOverall Summary:")
print(f"Total cluster members: {len(all_cluster_members)}")
print(f"Total field galaxies: {len(all_field_galaxies)}")
```

```{code-cell} ipython3
# Create simplified Y-H vs H color-magnitude diagram with density contours
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

# Calculate Y-H color and H magnitude
def calculate_color_magnitude(df):
    """Convert uniform-aperture fluxes to AB magnitudes and compute Y-H color.

    Magnitudes are computed as ``m = -2.5 * log10(flux) + 23.9``, where the
    zero-point 23.9 corresponds to fluxes in microjanskies — the unit of the
    ``flux_*_unif`` columns in the Euclid Q1 photo-z catalogue. The resulting
    magnitudes are not K-corrected or corrected for Milky Way extinction.

    Parameters
    ----------
    df : `pandas.DataFrame`
        DataFrame containing columns ``flux_y_unif`` and ``flux_h_unif`` in
        microjanskies. Both must be positive (no non-detections).

    Returns
    -------
    df : `pandas.DataFrame`
        Copy of the input with three additional columns: ``H_mag``, ``Y_mag``
        (AB magnitudes), and ``Y_H_color``.
    """
    df = df.copy()

    # Convert fluxes to magnitudes (using -2.5 * log10(flux))
    # Note: These are instrumental magnitudes, not absolute magnitudes
    df['H_mag'] = -2.5 * np.log10(df['flux_h_unif'])+23.9
    df['Y_mag'] = -2.5 * np.log10(df['flux_y_unif'])+23.9

    # Calculate Y-H color
    df['Y_H_color'] = df['Y_mag'] - df['H_mag']

    return df

# Calculate color-magnitude properties using only the needed fluxes
cluster_cmd = calculate_color_magnitude(all_cluster_members[['flux_y_unif', 'flux_h_unif']])
field_cmd = calculate_color_magnitude(all_field_galaxies[['flux_y_unif', 'flux_h_unif']])

# Remove outliers using plot boundaries (much simpler and more intuitive)
def remove_outliers_bounds(df, h_min=17, h_max=25, yh_min=-0.5, yh_max=1.5):
    """Filter a colour-magnitude table to a physically motivated range.

    Retains only rows whose H-band magnitude and Y-H colour fall within the
    specified boundaries. The defaults bracket the expected locus of cluster
    galaxies at z~0.4 in Euclid Q1 data: objects brighter than H=17 are likely
    stars or saturated; objects fainter than H=25 have unreliable photo-z and
    colours; Y-H outside [-0.5, 1.5] typically indicates bad photometry or
    non-galaxy contaminants.

    Parameters
    ----------
    df : `pandas.DataFrame`
        DataFrame with columns ``H_mag`` and ``Y_H_color``.
    h_min : float, optional
        Minimum H magnitude (brightest limit). Default is 17.
    h_max : float, optional
        Maximum H magnitude (faintest limit). Default is 25.
    yh_min : float, optional
        Minimum Y-H colour. Default is -0.5.
    yh_max : float, optional
        Maximum Y-H colour. Default is 1.5.

    Returns
    -------
    df_clean : `pandas.DataFrame`
        Filtered copy of the input DataFrame.
    """
    df_clean = df.copy()
    df_clean = df_clean[
        (df_clean['H_mag'] >= h_min) & (df_clean['H_mag'] <= h_max) &
        (df_clean['Y_H_color'] >= yh_min) & (df_clean['Y_H_color'] <= yh_max)
    ]
    return df_clean

# Remove outliers from both populations (using plot boundaries)
print("Removing outliers outside plot boundaries...")
cluster_cmd_clean = remove_outliers_bounds(cluster_cmd)
field_cmd_clean = remove_outliers_bounds(field_cmd)

print(f"Cluster galaxies: {len(cluster_cmd)} -> {len(cluster_cmd_clean)} (removed {len(cluster_cmd) - len(cluster_cmd_clean)})")
print(f"Field galaxies: {len(field_cmd)} -> {len(field_cmd_clean)} (removed {len(field_cmd) - len(field_cmd_clean)})")

# Plot comparison with lower alpha for better visibility
ax.scatter(field_cmd_clean['H_mag'], field_cmd_clean['Y_H_color'],
          c='blue', alpha=0.1, s=25, label=f'Field ({len(field_cmd_clean)})')
ax.scatter(cluster_cmd_clean['H_mag'], cluster_cmd_clean['Y_H_color'],
          c='red', alpha=0.1, s=30, label=f'Cluster ({len(cluster_cmd_clean)})')

# Add density contours to show the tightness of each population

# Create density contours for field galaxies (using cleaned data)
if len(field_cmd_clean) > 10:  # Need enough points for meaningful contours
    field_h = field_cmd_clean['H_mag'].values
    field_yh = field_cmd_clean['Y_H_color'].values
    field_xy = np.vstack([field_h, field_yh])
    field_density = gaussian_kde(field_xy)

    # Create grid for contour plot
    h_grid = np.linspace(17, 25, 100)
    yh_grid = np.linspace(-0.5, 1.5, 100)
    H_grid, YH_grid = np.meshgrid(h_grid, yh_grid)
    positions = np.vstack([H_grid.ravel(), YH_grid.ravel()])
    field_z = np.reshape(field_density(positions).T, H_grid.shape)

    # Plot field contours
    ax.contour(H_grid, YH_grid, field_z, levels=3, colors='blue', alpha=0.6, linestyles='--', linewidths=1)

# Create density contours for cluster galaxies (using cleaned data)
if len(cluster_cmd_clean) > 10:  # Need enough points for meaningful contours
    cluster_h = cluster_cmd_clean['H_mag'].values
    cluster_yh = cluster_cmd_clean['Y_H_color'].values
    cluster_xy = np.vstack([cluster_h, cluster_yh])
    cluster_density = gaussian_kde(cluster_xy)

    # Use same grid
    cluster_z = np.reshape(cluster_density(positions).T, H_grid.shape)

    # Plot cluster contours
    ax.contour(H_grid, YH_grid, cluster_z, levels=3, colors='red', alpha=0.8, linestyles='-', linewidths=2)

ax.set_xlabel('H Magnitude')
ax.set_ylabel('Y-H Color')
ax.set_title('Cluster vs Field Color-Magnitude Diagram\n(with density contours)')
ax.grid(True, alpha=0.3)
ax.legend()

# Set axis limits as requested
ax.set_xlim(17, 25)  # H magnitude range
ax.set_ylim(-0.5, 1.5)  # Y-H color range

plt.tight_layout()
plt.show()
```

**Figure 3. Y−H versus H colour–magnitude diagram for the cluster candidates compared to the control field.**
We construct an H-band magnitude and a Y−H colour from the Euclid Q1 uniform-aperture fluxes (`flux_h_unif`, `flux_y_unif`) by converting microjansky fluxes to AB magnitudes. Points show individual objects: cluster members (red) and field galaxies (blue).

To highlight the dominant loci beyond the sparse scatter, we overlay density contours derived from a 2D Gaussian kernel density estimate (KDE) computed in the \((H,\,Y-H)\) plane separately for the cluster and field samples.
The KDE is evaluated on a regular grid spanning the plotted limits, and three contour levels are drawn for each population (solid red for cluster; dashed blue for field).
A genuine cluster is expected to show a relatively **tighter and/or shifted colour–magnitude locus** compared to the general field population (e.g., a concentration consistent with a red-sequence-like population), while the field sample traces the broader distribution of galaxies along the line of sight.

+++

## 7. Spectral Analysis

We extract 1D spectra for cluster and field galaxies from the Euclid spectroscopic data to analyze their emission line properties and star formation activity.
The analysis includes spectral smoothing, flux normalization, and identification of nebular emission lines (Hα, Hβ, [OII], [OIII], etc.) that fall within the observed near-infrared wavelength range at the cluster's redshift. Note: This section is computationally intensive and may be skipped for initial analysis.

```{code-cell} ipython3
spectra_cache_dir = "data/irsa_spectra"
os.makedirs(spectra_cache_dir, exist_ok=True)

BUCKET_NAME = "nasa-irsa-euclid-q1"
table_1dspectra = "euclid.objectid_spectrafile_association_q1"

cluster_object_ids = all_cluster_members["object_id"].tolist()
field_object_ids   = all_field_galaxies["object_id"].tolist()

def get_n_spectra(obj_ids, n=10):
    """
    Fetch up to n spectra for the given object_ids, stopping early once n are found.
    Uses one TAP query, groups by FITS file, caches per-object results in ./data/irsa_spectra/.
    Skips TAP rows with missing (masked) HDU indices.
    """
    obj_ids = [int(x) for x in obj_ids]
    if len(obj_ids) == 0 or n <= 0:
        return {}

    # Cache-first
    spectra = {}
    remaining = []
    for oid in obj_ids:
        cache_npz = os.path.join(spectra_cache_dir, f"{oid}.npz")
        if os.path.exists(cache_npz):
            d = np.load(cache_npz)
            spectra[oid] = {
                "wave": d["wave"] * u.angstrom,
                "flux": d["flux"] * (u.erg / u.s / u.cm**2 / u.angstrom),
                "error": d["error"] * (u.erg / u.s / u.cm**2 / u.angstrom),
                "object_id": oid,
            }
            if len(spectra) >= n:
                return dict(list(spectra.items())[:n])
        else:
            remaining.append(oid)

    if len(remaining) == 0:
        return dict(list(spectra.items())[:n])

    # TAP query once
    id_list = ",".join(map(str, remaining))
    adql_query = f"""
    SELECT objectid, path, hdu
    FROM {table_1dspectra}
    WHERE objectid IN ({id_list})
    """

    try:
        assoc = Irsa.query_tap(adql_query).to_table()
    except (DALQueryError, requests.exceptions.RequestException) as e:
        print("TAP query failed:", e)
        return dict(list(spectra.items())[:n])

    if len(assoc) == 0:
        print("No spectrum associations returned.")
        return dict(list(spectra.items())[:n])

    # Group by FITS file, skipping masked/invalid hdu
    groups = {}
    skipped_hdu = 0

    for row in assoc:
        obj_id = int(row["objectid"])
        uri = str(row["path"])

        # HDU can be masked for some rows -> skip
        hdu_val = row["hdu"]
        try:
            # This will fail for masked values
            hdu_index = int(hdu_val)
        except Exception:
            skipped_hdu += 1
            continue

        spec_fpath_key = uri.replace("api/spectrumdm/convert/euclid/", "").split("?")[0]
        s3_uri = f"s3://{BUCKET_NAME}/{spec_fpath_key}"
        groups.setdefault(s3_uri, []).append((obj_id, hdu_index))

    if skipped_hdu > 0:
        print(f"Skipped {skipped_hdu} TAP rows with missing/invalid HDU indices.")

    if len(groups) == 0:
        print("No usable (path,hdu) associations after filtering.")
        return dict(list(spectra.items())[:n])

    # Open files with progress bar
    for s3_uri, items in tqdm(list(groups.items()), desc="Opening FITS files", unit="file"):
        if len(spectra) >= n:
            break

        try:
            with fits.open(s3_uri, fsspec_kwargs={"anon": True}, lazy_load_hdus=True) as hdul:
                for (obj_id, hdu_index) in items:
                    if len(spectra) >= n:
                        break

                    try:
                        spec = QTable.read(hdul[hdu_index], format="fits")
                        header = hdul[hdu_index].header

                        fscale = header.get("FSCALE", 1.0)
                        wave = np.asarray(spec["WAVELENGTH"]) * u.angstrom
                        signal = np.asarray(spec["SIGNAL"])
                        var = np.asarray(spec["VAR"])
                        mask = np.asarray(spec["MASK"])

                        if not (len(wave) == len(signal) == len(var) == len(mask)):
                            continue

                        valid = (mask % 2 == 0) & (mask < 64)
                        wave = wave[valid]
                        flux = signal[valid] * fscale * u.erg / u.s / u.cm**2 / u.angstrom
                        error = np.sqrt(var[valid]) * fscale * flux.unit

                        spectra[obj_id] = {
                            "wave": wave,
                            "flux": flux,
                            "error": error,
                            "object_id": obj_id,
                        }

                        np.savez(
                            os.path.join(spectra_cache_dir, f"{obj_id}.npz"),
                            wave=wave.value,
                            flux=flux.value,
                            error=error.value,
                        )

                    except Exception:
                        continue

        except Exception as e:
            print(f"Failed to open {s3_uri}: {e}")
            continue

    return dict(list(spectra.items())[:n])

# Run (it will stop as soon as it finds 10 real spectra)
cluster_spectra = get_n_spectra(cluster_object_ids, n=10)
field_spectra   = get_n_spectra(field_object_ids, n=10)

print(f"Cluster spectra retrieved: {len(cluster_spectra)}")
print(f"Field spectra retrieved:   {len(field_spectra)}")
print(f"Cache dir: {spectra_cache_dir}/")
```

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, median_filter

# --- your line lists (kept as-is, but note: duplicate keys get overwritten in dicts) ---
optical_lines = {
    'Hα': 6563,
    'Hβ': 4861,
    'Hγ': 4340,
    'Hδ': 4102,
    '[OII]': 3727,
    '[OIII]5007': 5007,
    '[OIII]4959': 4959,
    '[NII]': 6583,
    '[SII]6717': 6717,
    '[SII]6731': 6731,
    '[NeIII]': 3869,
    '[ArIII]': 7136
}

nir_lines = {
    'Paα': 18751,
    'Paβ': 12818,
    'Paγ': 10938,
    'Brα': 40512,
    'Brβ': 26252,
    '[SIII]': 9532,
    '[FeII]12570': 12570,
    '[FeII]16440': 16440
}

emission_lines = {**optical_lines, **nir_lines}

# --- parameters you can tune quickly ---
sigma_smooth = 1.5        # 0 to disable smoothing
cont_window = 151          # running-median window in pixels (odd recommended)
norm_percentile = 95       # robust scale from residuals
n_grid = 900               # common wavelength grid points for stacking

# If you have these from selection, use them. Otherwise keep your old ±0.12 approach:
cluster_z = float(cluster['zPZWav']) if 'zPZWav' in cluster else float(cluster.get('zPZ', np.nan))
redshift_width = 0.12
z_min, z_max = cluster_z - redshift_width, cluster_z + redshift_width

def preprocess_spectrum(spec):
    """Continuum-remove + robust-normalize a spectrum for visualization."""
    w = np.asarray(spec['wave'].value, float)
    f = np.asarray(spec['flux'].value, float)
    e = np.asarray(spec['error'].value, float) if 'error' in spec else None

    good = np.isfinite(w) & np.isfinite(f)
    if e is not None:
        good &= np.isfinite(e) & (e > 0)

    w = w[good]
    f = f[good]
    if e is not None:
        e = e[good]

    if w.size < 30:
        return None

    # sort by wavelength (just in case)
    s = np.argsort(w)
    w, f = w[s], f[s]
    if e is not None:
        e = e[s]

    # light smoothing to suppress pixel-to-pixel noise
    f_s = gaussian_filter1d(f, sigma=sigma_smooth) if (sigma_smooth and sigma_smooth > 0) else f

    # running median continuum
    k = int(cont_window)
    if k % 2 == 0:
        k += 1
    if k >= f_s.size:
        k = max(5, (f_s.size // 2) * 2 - 1)

    cont = median_filter(f_s, size=k, mode='nearest')
    resid = f_s - cont

    # robust normalization so objects are comparable in amplitude
    scale = np.nanpercentile(np.abs(resid), norm_percentile)
    if not np.isfinite(scale) or scale <= 0:
        scale = np.nanstd(resid)
    if not np.isfinite(scale) or scale <= 0:
        return None

    y = resid / scale
    return w, y

def build_stack(spectra_dict, w_grid):
    """Interpolate processed spectra onto a common grid and return matrix [nobj, ngrid]."""
    Ys = []
    for obj_id, spec in spectra_dict.items():
        out = preprocess_spectrum(spec)
        if out is None:
            continue
        w, y = out
        # interpolate to common grid; outside range -> NaN so it doesn't bias stats
        y_i = np.interp(w_grid, w, y, left=np.nan, right=np.nan)
        Ys.append(y_i)

    if len(Ys) == 0:
        return None

    Y = np.vstack(Ys)
    med = np.nanmedian(Y, axis=0)
    p16 = np.nanpercentile(Y, 16, axis=0)
    p84 = np.nanpercentile(Y, 84, axis=0)
    return Y, med, p16, p84

def lines_in_window(wmin, wmax, z):
    """Return dict of lines whose observed wavelength falls in [wmin, wmax]."""
    out = {}
    for name, rest in emission_lines.items():
        obs = rest * (1 + z)
        if wmin <= obs <= wmax:
            out[name] = obs
    return out

# --- determine a common observed-frame wavelength grid from available spectra ---
all_w = []
for d in [cluster_spectra, field_spectra]:
    for spec in d.values():
        w = np.asarray(spec['wave'].value, float)
        w = w[np.isfinite(w)]
        if w.size:
            all_w.append([np.nanmin(w), np.nanmax(w)])

if len(all_w) == 0:
    raise RuntimeError("No valid wavelengths found in cluster_spectra/field_spectra.")

wmin = min(a for a, b in all_w)
wmax = max(b for a, b in all_w)
w_grid = np.linspace(wmin, wmax, n_grid)

# --- build stacks ---
cluster_stack = build_stack(cluster_spectra, w_grid)
field_stack   = build_stack(field_spectra,   w_grid)

# --- plot: two panels, close to your original, but readable ---
def label_emission_line(ax, x, label, y=0.92, rotation=90, color='k'):
    """Place a vertical emission-line label at wavelength x."""
    ax.text(
        x, y, label,
        transform=ax.get_xaxis_transform(),
        rotation=rotation,
        color=color,
        fontsize=10,
        ha='right',
        va='top',
        bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=1)
    )

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

def plot_panel(ax, spectra_dict, stack, color, title):
    # individual processed spectra (faint)
    for obj_id, spec in spectra_dict.items():
        out = preprocess_spectrum(spec)
        if out is None:
            continue
        w, y = out
        ax.plot(w, y, color=color, alpha=0.35, lw=1)

    # median + spread (this is the key improvement)
    if stack is not None:
        Y, med, p16, p84 = stack
        ax.plot(w_grid, med, color=color, lw=2.5, alpha=0.9, label=f"median (n={Y.shape[0]})")
        ax.fill_between(w_grid, p16, p84, color=color, alpha=0.15, label="16–84%")

    ax.axhline(0, color="k", lw=1, alpha=0.25)
    ax.set_ylabel("Continuum-subtracted, robust-normalized flux")
    ax.set_title(title)
    ax.grid(True, alpha=0.35)
    ax.legend(loc="upper right")

plot_panel(ax1, cluster_spectra, cluster_stack, "tab:blue", "Cluster spectra (continuum-subtracted + robust-normalized)")
plot_panel(ax2, field_spectra,   field_stack,   "tab:green", "Control/field spectra (continuum-subtracted + robust-normalized)")

# --- line markers + z-slice shading (same concept you had) ---
lines_here = lines_in_window(wmin, wmax, cluster_z)
colors = ['red', 'orange', 'purple', 'brown', 'magenta', 'gray', 'cyan', 'goldenrod']

for i, (line_name, obs_wave) in enumerate(lines_here.items()):
    c = colors[i % len(colors)]
    rest = emission_lines[line_name]

    wave_min = rest * (1 + z_min)
    wave_max = rest * (1 + z_max)

    for ax in (ax1, ax2):
        ax.axvspan(wave_min, wave_max, alpha=0.07, color=c)
        ax.axvline(obs_wave, color=c, ls="--", lw=1.2, alpha=0.8)
        label_emission_line(
            ax1,
            obs_wave,
            f"{line_name}",
            color=c
        )

# cosmetic
ax2.set_xlabel("Observed wavelength (Å)")
plt.tight_layout()
plt.show()
```

**Figure 4 — Continuum-subtracted, normalized Euclid Q1 spectra for cluster (top) and control (bottom) samples.**
Thin lines show individual galaxy spectra; thick curves indicate the median, with shaded regions marking the 16–84 percentile range.
Dashed vertical lines mark the expected observed wavelengths of prominent nebular emission lines at the cluster redshift, while shaded bands indicate the wavelength range allowed by the finite redshift slice of the sample.
Spectra are not stacked in redshift; the figure is intended to highlight relative excess emission in physically motivated regions.
Euclid Q1 spectra contain known instrumental artifacts that are addressed in DR1, so better results to be seen in the future.

+++

## 8. NED Database Search

We search the NASA/IPAC Extragalactic Database (NED) for information about our cluster center, field center, and cluster member galaxies.
This might provide additional information on cluster members. We search a smaller radius of 3 arcmin to avoid NED timeout.

```{code-cell} ipython3
# Get coordinates from our analysis
cluster_ra = cluster['RAPZWav']  # 60.4686 degrees
cluster_dec = cluster['DecPZWav']  # -50.4780 degrees
cluster_z = cluster['zPZWav']  # 0.43

# Search NED for cluster center
print("=== SEARCHING NED FOR CLUSTER CENTER ===")
try:
    cluster_center_results = Ned.query_region(SkyCoord(ra=cluster_ra, dec=cluster_dec, unit='deg'),
                                            radius=5*u.arcsec)
    print(f"Found {len(cluster_center_results)} objects within 5 arcsec of cluster center")
    if len(cluster_center_results) > 0:
        print("\\nCluster center objects:")
        for i, obj in enumerate(cluster_center_results[:5]):  # Show first 5
            print(f"  {i+1}. {obj['Object Name']} - Type: {obj['Type']} - z: {obj['Redshift']}")
    else:
        print("No objects found in NED at cluster center")
except (RemoteServiceError, Timeout, ConnectionError) as e:
    print(f"Error searching cluster center: {e}")

print()

# Search NED for control field center
print("=== SEARCHING NED FOR CONTROL FIELD CENTER ===")
try:
    control_center_results = Ned.query_region(SkyCoord(ra=control_ra, dec=control_dec, unit='deg'),
                                            radius=5*u.arcsec)
    print(f"Found {len(control_center_results)} objects within 5 arcsec of control field center")
    if len(control_center_results) > 0:
        print("\\nControl field center objects:")
        for i, obj in enumerate(control_center_results[:5]):  # Show first 5
            print(f"  {i+1}. {obj['Object Name']} - Type: {obj['Type']} - z: {obj['Redshift']}")
    else:
        print("No objects found in NED at control field center")
except (RemoteServiceError, Timeout, ConnectionError) as e:
    print(f"Error searching control field center: {e}")

print()
```

```{code-cell} ipython3
# Cluster center + redshift slice
cluster_ra  = cluster['RAPZWav']
cluster_dec = cluster['DecPZWav']
cluster_z   = cluster['zPZWav']
z_min, z_max = cluster_z - 0.06, cluster_z + 0.06

# Try NED search with retry logic
max_retries = 3
for attempt in range(max_retries):
    try:
        print(f"NED search attempt {attempt + 1}/{max_retries}...")

        # Search NED within 3 arcmin of cluster center
        results = Ned.query_region(
            SkyCoord(ra=cluster_ra, dec=cluster_dec, unit='deg'),
            radius=3 * u.arcmin
        )

        print(f"\nNED search results:")
        print(f"  Total objects found: {len(results)}")

        # Count objects with redshift information
        has_z = np.isfinite(np.asarray(results['Redshift'], float))
        objects_with_z = results[has_z]
        print(f"  Objects with redshift: {len(objects_with_z)}")

        # Count objects in cluster redshift range
        in_z_range = (objects_with_z['Redshift'] >= z_min) & (objects_with_z['Redshift'] <= z_max)
        cluster_objects = objects_with_z[in_z_range]
        print(f"  Objects in z={z_min:.2f}-{z_max:.2f}: {len(cluster_objects)}")

        # --- robustly choose RA/Dec columns (NED output varies) ---
        print("\nNED columns:", results.colnames)

        ra_candidates  = ['RA(deg)', 'RA_deg', 'RA', 'RAJ2000']
        dec_candidates = ['DEC(deg)', 'DEC_deg', 'DEC', 'DEJ2000']

        def pick_col(table, candidates):
            for c in candidates:
                if c in table.colnames:
                    return c
            return None

        ra_col  = pick_col(results, ra_candidates)
        dec_col = pick_col(results, dec_candidates)

        if ra_col is None or dec_col is None:
            raise KeyError(f"Could not find RA/Dec columns in NED table. Columns are: {results.colnames}")

        print(f"Using RA column: {ra_col}")
        print(f"Using Dec column: {dec_col}")

        # Decide whether RA/Dec are in degrees or sexagesimal strings
        # (If the column name explicitly indicates degrees, treat as deg/deg; otherwise assume RA is hourangle, Dec is deg.)
        use_deg = (('deg' in ra_col.lower()) or ra_col.lower().endswith('_deg')) and (('deg' in dec_col.lower()) or dec_col.lower().endswith('_deg'))
        coord_unit = ('deg', 'deg') if use_deg else (u.hourangle, u.deg)

        if len(cluster_objects) > 0:
            print("\nObjects in cluster redshift range:")

            cluster_coord = SkyCoord(ra=cluster_ra, dec=cluster_dec, unit='deg')

            for obj in cluster_objects:
                obj_coord = SkyCoord(ra=obj[ra_col], dec=obj[dec_col], unit=coord_unit)
                sep = cluster_coord.separation(obj_coord).to(u.arcmin).value
                print(f"  {obj['Object Name']} - {obj['Type']} - z={obj['Redshift']:.3f} - {sep:.1f}'")

            # Type counts (Astropy-friendly)
            types = np.array([str(t) for t in cluster_objects['Type']], dtype=str)
            unique, counts = np.unique(types, return_counts=True)
            print("\nTypes:", dict(zip(unique, counts)))
        else:
            print("No objects found in cluster redshift range")

        break  # Success, exit retry loop

    except (Timeout, ConnectionError) as e:
        print(f"Timeout/connection error (attempt {attempt + 1}): {e}")
        if attempt < max_retries - 1:
            print("Retrying in 5 seconds...")
            time.sleep(5)
        else:
            print("NED search failed after all retries")
    except (RemoteServiceError, OSError, KeyError) as e:
        print(f"NED search error: {e}")
        break
```

```{code-cell} ipython3
# Try NED search with retry logic for control field
max_retries = 3
for attempt in range(max_retries):
    try:
        # Search NED within 5 arcmin of control field center
        control_results = Ned.query_region(SkyCoord(ra=control_ra, dec=control_dec, unit='deg'), radius=3*u.arcmin)

        # Filter for objects with redshift in cluster range
        has_z = ~np.isnan(control_results['Redshift'])
        in_z_range = (control_results['Redshift'] >= z_min) & (control_results['Redshift'] <= z_max)
        control_objects = control_results[has_z & in_z_range]

        print(f"Control field NED search: {len(control_results)} total objects, {len(control_objects)} in z={z_min:.2f}-{z_max:.2f}")

        if len(control_objects) > 0:
            print("\\nObjects in control field redshift range:")
            for obj in control_objects:
                obj_coord = SkyCoord(ra=obj['RA(deg)'], dec=obj['DEC(deg)'], unit='deg')
                control_coord = SkyCoord(ra=control_ra, dec=control_dec, unit='deg')
                sep = control_coord.separation(obj_coord).to(u.arcmin).value
                print(f"  {obj['Object Name']} - {obj['Type']} - z={obj['Redshift']:.3f} - {sep:.1f}'")

            print(f"\\nControl field types: {dict(control_objects['Type'].value_counts())}")
        else:
            print("No objects found in control field redshift range")

        break  # Success, exit retry loop

    except (Timeout, ConnectionError) as e:
        if attempt < max_retries - 1:
            time.sleep(5)
        else:
            print("Control field NED search failed after all retries")
    except (RemoteServiceError, OSError) as e:
        print(f"Control field NED search error: {e}")
        break
```

```{code-cell} ipython3
# Overlay NED sources on cluster member visualization
if 'cluster_objects' in locals() and len(cluster_objects) > 0:
    # Create the same plot as before but with NED sources added
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8),
                                   subplot_kw={'projection': cluster_cutout_wcs})

    # Cluster field (left subplot) - make background much more transparent
    ax1.imshow(cluster_rgb, origin='lower', alpha=0.3)
    cluster_unique_labels = set(cluster_labels)
    cluster_colors = plt.cm.Spectral(np.linspace(0, 1, len(cluster_unique_labels)))

    for k, col in zip(cluster_unique_labels, cluster_colors):
        if k == -1:
            col = 'black'
            marker = 'o'
            size = 10
            alpha = 0.3
        else:
            marker = 'o'
            size = 20
            alpha = 0.8

        class_member_mask = (cluster_labels == k)
        xy = cluster_galaxy_coords[class_member_mask]

        if len(xy) > 0:
            ax1.scatter(xy[:, 0], xy[:, 1], c='red', marker=marker, s=size, alpha=alpha,
                       edgecolors='white', linewidth=0.5)

    # Add NED sources in cluster field
    for obj in cluster_objects:
        # Find coordinate columns
        ra_col = None
        dec_col = None

        for col in cluster_objects.colnames:
            col_lower = col.lower()
            if ('ra' in col_lower and ('deg' in col_lower or 'degree' in col_lower)) or col_lower == 'ra':
                ra_col = col
            elif ('dec' in col_lower and ('deg' in col_lower or 'degree' in col_lower)) or col_lower == 'dec':
                dec_col = col
            elif 'right' in col_lower and 'ascension' in col_lower:
                ra_col = col
            elif 'declination' in col_lower:
                dec_col = col

        if ra_col and dec_col:
            ned_coord = SkyCoord(ra=obj[ra_col], dec=obj[dec_col], unit='deg')
            ned_pixel = cluster_cutout_wcs.world_to_pixel(ned_coord)

            # Check if NED source is within the image bounds
            if (0 <= ned_pixel[0] < cluster_rgb.shape[1] and
                0 <= ned_pixel[1] < cluster_rgb.shape[0]):
                ax1.scatter(ned_pixel[0], ned_pixel[1], c='none', marker='o', s=50,
                           alpha=0.9, edgecolors='blue', linewidth=3,
                           label=f"NED {obj['Type']} (z={obj['Redshift']:.3f})")
        else:
            # Try to use first two numeric columns as fallback
            numeric_cols = []
            for col in cluster_objects.colnames:
                try:
                    float(obj[col])
                    numeric_cols.append(col)
                except:
                    pass
            if len(numeric_cols) >= 2:
                try:
                    ned_coord = SkyCoord(ra=obj[numeric_cols[0]], dec=obj[numeric_cols[1]], unit='deg')
                    ned_pixel = cluster_cutout_wcs.world_to_pixel(ned_coord)
                    if (0 <= ned_pixel[0] < cluster_rgb.shape[1] and
                        0 <= ned_pixel[1] < cluster_rgb.shape[0]):
                        ax1.scatter(ned_pixel[0], ned_pixel[1], c='none', marker='o', s=200,
                                   alpha=0.9, edgecolors='blue', linewidth=3,
                                   label=f"NED {obj['Type']} (z={obj['Redshift']:.3f})")
                except:
                    pass

    ax1.set_title('Cluster Field with NED Sources', fontsize=14, fontweight='bold')
    ax1.set_xlabel('RA')
    ax1.set_ylabel('Dec')

    # Control field (right subplot) - make background more transparent for comparison
    ax2.imshow(control_rgb, origin='lower', alpha=0.3)
    control_unique_labels = set(control_labels)
    control_colors = plt.cm.Spectral(np.linspace(0, 1, len(control_unique_labels)))

    for k, col in zip(control_unique_labels, control_colors):
        if k == -1:
            col = 'black'
            marker = 'o'
            size = 10
            alpha = 0.3
        else:
            marker = 'o'
            size = 20
            alpha = 0.8

        class_member_mask = (control_labels == k)
        xy = control_galaxy_coords[class_member_mask]

        if len(xy) > 0:
            ax2.scatter(xy[:, 0], xy[:, 1], c='red', marker=marker, s=size, alpha=alpha,
                       edgecolors='white', linewidth=0.5)

    # Add NED sources in control field if any were found
    if 'control_objects' in locals() and len(control_objects) > 0:
        for obj in control_objects:
            # Find coordinate columns
            ra_col = None
            dec_col = None

            for col in control_objects.colnames:
                col_lower = col.lower()
                if ('ra' in col_lower and ('deg' in col_lower or 'degree' in col_lower)) or col_lower == 'ra':
                    ra_col = col
                elif ('dec' in col_lower and ('deg' in col_lower or 'degree' in col_lower)) or col_lower == 'dec':
                    dec_col = col
                elif 'right' in col_lower and 'ascension' in col_lower:
                    ra_col = col
                elif 'declination' in col_lower:
                    dec_col = col

            if ra_col and dec_col:
                ned_coord = SkyCoord(ra=obj[ra_col], dec=obj[dec_col], unit='deg')
                ned_pixel = control_cutout_wcs.world_to_pixel(ned_coord)

                # Check if NED source is within the image bounds
                if (0 <= ned_pixel[0] < control_rgb.shape[1] and
                    0 <= ned_pixel[1] < control_rgb.shape[0]):
                    ax2.scatter(ned_pixel[0], ned_pixel[1], c='none', marker='o', s=200,
                               alpha=0.9, edgecolors='blue', linewidth=3,
                               label=f"NED {obj['Type']} (z={obj['Redshift']:.3f})")
            else:
                # Try to use first two numeric columns as fallback
                numeric_cols = []
                for col in control_objects.colnames:
                    try:
                        float(obj[col])
                        numeric_cols.append(col)
                    except:
                        pass
                if len(numeric_cols) >= 2:
                    try:
                        ned_coord = SkyCoord(ra=obj[numeric_cols[0]], dec=obj[numeric_cols[1]], unit='deg')
                        ned_pixel = control_cutout_wcs.world_to_pixel(ned_coord)
                        if (0 <= ned_pixel[0] < control_rgb.shape[1] and
                            0 <= ned_pixel[1] < control_rgb.shape[0]):
                            ax2.scatter(ned_pixel[0], ned_pixel[1], c='none', marker='o', s=200,
                                       alpha=0.9, edgecolors='blue', linewidth=3,
                                       label=f"NED {obj['Type']} (z={obj['Redshift']:.3f})")
                    except:
                        pass

    ax2.set_title('Control Field', fontsize=14, fontweight='bold')
    ax2.set_xlabel('RA')
    ax2.set_ylabel('Dec')

    # Add legend for NED sources
    if len(cluster_objects) > 0:
        # Get unique NED types for legend
        ned_types = set(obj['Type'] for obj in cluster_objects)
        handles = []
        labels = []

        handles.append(plt.Line2D([0], [0], marker='o', color='none',
                                linestyle='None', markersize=3,
                                markeredgecolor='blue', markeredgewidth=2))
        labels.append(f"NED Galaxy")


        handles.append(plt.Line2D([0], [0], marker='o', color='red',
                                linestyle='None', markersize=1, alpha=0.3))
        labels.append(r'Euclid Galaxy in $\Delta$z')
        # Add cluster member legend

        handles.append(plt.Line2D([0], [0], marker='o', color='red',
                                linestyle='None', markersize=8, alpha=0.8))
        labels.append('Euclid Cluster Members')

        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.95),
                  ncol=len(handles), fontsize=10)

    plt.tight_layout()
    plt.show()
```

**Figure 5- Euclid and NED galaxy distributions in cluster and control fields.**
Left: Cluster field showing Euclid-selected galaxies in the cluster redshift slice (red) and NED galaxies with spectroscopic redshift (blue).
Right: matched control field showing no confirmed NED detection in the redshift slice of the cluster.

```{code-cell} ipython3
# Cross-match Euclid and NED sources for redshift comparison

# Get Euclid sources in cluster field
euclid_coords = SkyCoord(ra=cluster_members_cluster_field['ra'],
                        dec=cluster_members_cluster_field['dec'], unit='deg')
euclid_redshifts = cluster_members_cluster_field['phz_median'].values

# Get NED sources if available
if 'cluster_objects' in locals() and len(cluster_objects) > 0:
    # Find coordinate columns for NED
    ra_col = None
    dec_col = None

    for col in cluster_objects.colnames:
        col_lower = col.lower()
        if ('ra' in col_lower and ('deg' in col_lower or 'degree' in col_lower)) or col_lower == 'ra':
            ra_col = col
        elif ('dec' in col_lower and ('deg' in col_lower or 'degree' in col_lower)) or col_lower == 'dec':
            dec_col = col
        elif 'right' in col_lower and 'ascension' in col_lower:
            ra_col = col
        elif 'declination' in col_lower:
            dec_col = col

    if ra_col and dec_col:
        # Create NED coordinates
        ned_coords = SkyCoord(ra=cluster_objects[ra_col].data,
                             dec=cluster_objects[dec_col].data, unit='deg')
        ned_redshifts = cluster_objects['Redshift'].data

        # Cross-match within 1 arcsec
        idx_ned, idx_euclid, d2d, d3d = search_around_sky(ned_coords, euclid_coords, 2*u.arcsec)

        if len(idx_ned) > 0:
            print(f"Found {len(idx_ned)} matches between NED and Euclid sources within 1 arcsec")

            # Get matched redshifts
            ned_z_matched = ned_redshifts[idx_ned]
            euclid_z_matched = euclid_redshifts[idx_euclid]

            # Simple redshift comparison plot
            plt.figure(figsize=(8, 6))
            plt.scatter(euclid_z_matched, ned_z_matched, alpha=0.7, s=50)
            plt.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='Perfect match')

            # Add cluster redshift reference lines
            cluster_z = cluster['zPZWav']
            plt.axvline(cluster_z, color='orange', linestyle=':', alpha=0.7, label=f'Cluster z={cluster_z:.3f}')
            plt.axhline(cluster_z, color='orange', linestyle=':', alpha=0.7)

            plt.xlabel('Euclid Redshift')
            plt.ylabel('NED Redshift')
            plt.title('Redshift Comparison: Euclid vs NED')
            plt.legend()
            plt.grid(True, alpha=0.3)

            # Calculate and display statistics
            z_diff = ned_z_matched - euclid_z_matched
            mean_diff = np.mean(z_diff)
            std_diff = np.std(z_diff)

            plt.text(0.05, 0.95, f'Mean difference: {mean_diff:.4f}\nStd deviation: {std_diff:.4f}',
                    transform=plt.gca().transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            plt.tight_layout()
            plt.xlim(z_min-0.1,z_max+0.1)
            plt.ylim(z_min-0.1,z_max+0.1)
            plt.show()

```

**Figure 6- Euclid–NED redshift comparison for cross-matched galaxies within 1″**

```{code-cell} ipython3

```
