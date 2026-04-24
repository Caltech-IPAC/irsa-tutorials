---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
  language: python
---

# Quickstart: Accessing OpenUniverse2024 Data

## Learning Goals

By the end of this tutorial, you will be able to:

- Browse the OpenUniverse2024 data directories on S3
- Explore the structure of Roman and Rubin FITS image files
- Read the OpenUniverse2024 parquet catalogs (transient, galaxy, and galaxy flux)
- Query Roman and Rubin images covering a sky position using the IRSA SIA service

## Introduction

The [OpenUniverse2024](https://arxiv.org/abs/2501.05632) simulation suite delivers ~70 deg² of matched optical/infrared imagery for both the LSST Wide-Fast-Deep (WFD) and the Nancy Grace Roman Space Telescope high-latitude survey, producing roughly 400 TB of publicly available synthetic imaging and catalogs. All data are stored in the cloud (AWS S3) and can be accessed anonymously without any credentials.

This tutorial is a focused introduction to data access only. It covers the three main categories:

1. **Directory structure for FITS images** — Roman and Rubin simulated science images stored in S3
2. **Parquet catalogs** — transient (SNANA), galaxy, and galaxy-flux tables, indexed by HEALPix sky region
3. **Image search via SIA** — querying which images cover a given sky position using astroquery and the IRSA Simple Image Access service

No astrophysical analysis is performed here. For science workflows that build on these access patterns, see the [TDE Light Curve](TDE_light_curve) and [SED Fitting](SED_fit) tutorials in this repository.

### Instructions

This notebook is designed to be run sequentially from top to bottom. All code is self-contained and relies on publicly accessible data.

### Input

- OpenUniverse2024 Roman and Rubin images and catalogs on AWS S3 (`s3://nasa-irsa-simulations/`)

### Output

- A gallery of example Roman FITS images
- Summary of parquet catalog structure and contents
- A table of image files overlapping a chosen sky position

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy astropy s3fs photutils matplotlib pyarrow hpgeom astroquery
```

```{code-cell} ipython3
import numpy as np
import s3fs
from matplotlib import pyplot as plt
import pyarrow.fs
import pyarrow.parquet as pq
import hpgeom
import json
from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
```

## 1. Explore Directory Structure for FITS images

The OpenUniverse2024 data live on the cloud in a public AWS S3 bucket and can be accessed anonymously using `s3fs`. This section shows how to establish that connection, navigate the directory tree, and inspect the contents of a FITS image file.

In the path below, `simple_model` refers to the simulated images with noise and realistic instrument effects, as opposed to `truth` images which are noise-free. The `full` simulation covers the complete survey footprint; a smaller `preview` subset is also available. See the [OpenUniverse2024 paper](https://arxiv.org/abs/2501.05632) for details on the differences. A `pointing` is a unique Roman observation visit — each pointing corresponds to one placement of the 18-detector focal plane on the sky, producing up to 18 individual FITS files (one per detector).

```{code-cell} ipython3
# Create an anonymous (public read-only) connection to the NASA IRSA S3 bucket.
s3 = s3fs.S3FileSystem(anon=True)

# Top-level path components
BUCKET_NAME = "nasa-irsa-simulations"
OU_PREFIX = "openuniverse2024"
ROMAN_TDS_PREFIX = "roman/full/RomanTDS/images/simple_model"

# Pick one band to explore
BAND = "J129"
band_directory = f"{BUCKET_NAME}/{OU_PREFIX}/{ROMAN_TDS_PREFIX}/{BAND}"
```

The pointings available for a given band can be listed by calling `s3.ls` on the band directory.

```{code-cell} ipython3
# List all pointings available for the chosen band
all_pointings = [p.split("/")[-1] for p in s3.ls(band_directory)]
print(f"Found {len(all_pointings)} pointings in band {BAND}:")
print(all_pointings[:10], "...")
```

We pick one of these pointings to explore further.

```{code-cell} ipython3
# Select one pointing and list the files it contains
POINTING = "10190"
image_directory = f"{band_directory}/{POINTING}"

files = [f"s3://{f}" for f in s3.ls(image_directory)]
print(f"Found {len(files)} files in pointing {POINTING}")
```

```{code-cell} ipython3
# Open one FITS file and inspect its extensions
fname = files[0]
with fits.open(fname, use_fsspec=True, fsspec_kwargs={"anon": True}, memmap=False) as hdul:
    print(f"File: {fname}")
    print(f"Number of extensions: {len(hdul)}\n")
    hdul.info()
```

Each Roman TDS FITS file contains four extensions: a `primary` header with no data, followed by three 4088×4088 pixel planes — `SCI` (science image), `ERR` (per-pixel uncertainty), and `DQ` (data quality mask).

+++

Let's display a gallery of example images to get a sense of the data. 
Note this gallery can take about a minute to build.

```{code-cell} ipython3
def show_gallery(files, max_images=9):
    """
    Display a gallery of FITS images.

    Parameters
    ----------
    files : list of str
        List of S3 URIs to FITS files.
    max_images : int, optional
        Maximum number of images to display (default: 9).
    """
    n_images = min(len(files), max_images)
    ncols = n_images if n_images < 4 else 3
    nrows = (n_images + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for i, f in enumerate(files[:n_images]):
        with fits.open(f, fsspec_kwargs={"anon": True}, memmap=False) as hdul:
            data = hdul[1].data
            vmin, vmax = np.nanpercentile(data, [5, 99])
            axes[i].imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
            axes[i].set_title(f.split("/")[-1], fontsize=8)
            axes[i].axis("off")

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    plt.show()
```

```{code-cell} ipython3
# Display up to 3 images from the selected directory.
show_gallery(files, max_images=3)
```

## 2. Access the Parquet Catalogs

The OpenUniverse2024 catalogs are stored as [Apache Parquet](https://parquet.apache.org/) files, partitioned by [HEALPix](https://healpix.sourceforge.io/) sky region (nside=32, RING ordering). Each region has three file types:

1. `snana_{region}.parquet` — one row per simulated transient event (supernovae, TDEs, etc.), with event type (`model_name`) and host galaxy ID (`host_id`)
2. `galaxy_{region}.parquet` — host galaxy positions and physical properties
3. `galaxy_flux_{region}.parquet` — multi-band Roman and Rubin photometry for each galaxy

We first look for the correct catalog for the center of the Roman Time-Domain Survey(TDS).  The `region` number in each filename is the HEALPix pixel index. Because we know that the catalogs were built with nside=32 and RING ordering, we can convert sky coordinates to a region index using [`hpgeom`](https://hpgeom.readthedocs.io/en/latest/).

```{code-cell} ipython3
# The Roman Time-Domain Survey is centered near the LSST ELAIS-S1 Deep Drilling Field.
ra = 9.45
dec = -44.02

# Convert sky coordinates to a HEALPix region index (nside=32, RING ordering)
nside = 32
region = hpgeom.angle_to_pixel(nside, ra, dec, lonlat=True, nest=False)
print(f"HEALPix region for RA={ra}, Dec={dec}: {region}")
```

```{code-cell} ipython3
# Build the S3 paths for this region's catalog files
CATALOG_NAME = "roman_rubin_cats_v1.1.2_faint"
catalog_prefix = f"{BUCKET_NAME}/{OU_PREFIX}/roman/full/{CATALOG_NAME}"

snana_path    = f"{catalog_prefix}/snana_{region}.parquet"
galaxy_path   = f"{catalog_prefix}/galaxy_{region}.parquet"
gal_flux_path = f"{catalog_prefix}/galaxy_flux_{region}.parquet"

print("SNANA file:       ", snana_path)
print("Galaxy info file: ", galaxy_path)
print("Galaxy flux file: ", gal_flux_path)
```

### 2.1 Inspect the SNANA Transient Catalog

`inspect_parquet_columns()` reads only the Parquet metadata footer to print the row count and column names — no data is loaded into memory.
We use it here for the SNANA catalog and repeat it for the galaxy info and flux catalogs below.

```{code-cell} ipython3
def inspect_parquet_files(s3_path):
    """
    Print the structure of a Parquet file on S3 without reading its data.

    Reads only the Parquet metadata footer (row count, column names and types),
    which is fast regardless of file size.

    Parameters
    ----------
    s3_path : str
        S3 path to the Parquet file (without the s3:// prefix).
    """
    fs = pyarrow.fs.S3FileSystem(anonymous=True)
    meta = pq.read_metadata(s3_path, filesystem=fs)
    schema = pq.read_schema(s3_path, filesystem=fs)

    print(f"Rows: {meta.num_rows}  |  Columns: {len(schema.names)}")
    print("\nColumn names:")
    for name in schema.names:
        print("  ", name)
```

```{code-cell} ipython3
inspect_parquet_files(snana_path)
```

```{code-cell} ipython3
# Read just the model_name column to see what transient types are in this region
fs = pyarrow.fs.S3FileSystem(anonymous=True)
model_names = pq.read_table(snana_path, filesystem=fs, columns=["model_name"]).to_pandas()
model_names["model_name"].unique()
```

### 2.2 Inspect the Galaxy Info Catalog

```{code-cell} ipython3
inspect_parquet_files(galaxy_path)
```

### 2.3 Inspect the Galaxy Flux Catalog

```{code-cell} ipython3
inspect_parquet_files(gal_flux_path)
```

### 2.4 Join Transient Events to Their Host Galaxies

A common operation is to take a transient from the SNANA file and look up its host galaxy's sky position from the galaxy info file.
The two files share a common key: `host_id` in the SNANA file corresponds to `galaxy_id` in the galaxy info file.
We read the full SNANA catalog here, then use a filter to fetch only the matching row from the galaxy file without loading the entire galaxy catalog.

Note: The next cell takes ~45s to run

```{code-cell} ipython3
fs = pyarrow.fs.S3FileSystem(anonymous=True)

# Read the full SNANA catalog for this region
df_snana = pq.read_table(snana_path, filesystem=fs).to_pandas()

# Pick one transient — here we grab the first row as an example
example_transient = df_snana.iloc[0]
print("Example transient:")
print(example_transient[["model_name", "host_id", "start_mjd", "end_mjd"]])

# Look up its host galaxy by matching host_id (SNANA) to galaxy_id (galaxy info file)
host = pq.read_table(
    galaxy_path,
    filesystem=fs,
    filters=[("galaxy_id", "==", example_transient["host_id"])]
).to_pandas()

print("\nHost galaxy info:")
host
```

## 3. Image Search

Given a sky position (e.g., the host galaxy coordinates from Section 2), we can search for all Roman or Rubin images that cover that position using the IRSA Simple Image Access (SIA) service via `astroquery`.
First we set up the connection to the SIA service and list the available catalogs, then we query by position to get a table of matching images, and finally we extract the cloud locations (S3 URIs) so the files can be opened directly.

```{code-cell} ipython3
# Point the astroquery IRSA client to the correct locations.
Irsa.sia_url = "https://irsa.ipac.caltech.edu/SIA"
Irsa.tap_url = "https://irsa.ipac.caltech.edu/TAP"

# List all available simulated image collections
Irsa.list_collections(servicetype='SIA')
```

```{code-cell} ipython3
# Collection names for OpenUniverse2024
OU_ROMAN_SIA_COLLECTION = 'simulated_roman_openuniverse2024'
OU_RUBIN_SIA_COLLECTION = 'simulated_rubin_openuniverse2024'
```

```{code-cell} ipython3
def get_s3_fpath(cloud_access):
    """Extract the S3 URI from the cloud_access JSON string in an SIA result."""
    cloud_info = json.loads(cloud_access)
    bucket_name = cloud_info['aws']['bucket_name']
    key = cloud_info['aws']['key']
    return f's3://{bucket_name}/{key}'
```

```{code-cell} ipython3
# Use the host galaxy position from Section 2 (or set any RA/Dec you want to query).
host_ra  = float(host.iloc[0]["ra"])
host_dec = float(host.iloc[0]["dec"])
search_radius = 1 * u.arcsec  # small radius: we just need images that contain this point

#convert ra, dec to SkyCoords for ease of use
coords = SkyCoord(host_ra, host_dec, unit='deg')

# Query Roman TDS images in the J129 band
sia_results = Irsa.query_sia(pos=(coords, search_radius.to(u.deg)),
                             collection=OU_ROMAN_SIA_COLLECTION)

# We first choose to look at the J129 band and the simple_model images
bandname = "J129"
roman_images = sia_results[
    ['TDS_simple_model' in r['obs_id'] and bandname in r['energy_bandpassname']
     for r in sia_results]
]
roman_images['s3_uri'] = [get_s3_fpath(r['cloud_access']) for r in roman_images]

print(f"Found {len(roman_images)} Roman {bandname} images at RA={host_ra:.4f}, Dec={host_dec:.4f}")
roman_images['obs_id', 't_min', 't_max', 's3_uri']
```

```{code-cell} ipython3
# The same search works for Rubin images — just swap the collection name and band filter.
# Unlike the Roman collection, the Rubin collection contains only one image type (calexp),
# so no obs_id filter is needed beyond selecting the desired band.
rubin_band = "r"
rubin_results = Irsa.query_sia(pos=(coords, search_radius.to(u.deg)),
                               collection=OU_RUBIN_SIA_COLLECTION)

rubin_images = rubin_results[
    [rubin_band in r['energy_bandpassname'] for r in rubin_results]
]
rubin_images['s3_uri'] = [get_s3_fpath(r['cloud_access']) for r in rubin_images]

print(f"Found {len(rubin_images)} Rubin {rubin_band}-band images at RA={host_ra:.4f}, Dec={host_dec:.4f}")
rubin_images['obs_id', 't_min', 't_max', 's3_uri']
```

You now have S3 URIs for all Roman and Rubin images covering your target position. To open any of these images, pass the URI to `astropy.io.fits.open` with `fsspec_kwargs={"anon": True}` as shown in Section 1.

## Acknowledgements

- [IPAC-IRSA](https://irsa.ipac.caltech.edu/)
- This work made use of Astropy:\footnote{http://www.astropy.org} a community-developed core Python package and an ecosystem of tools and resources for astronomy.

## About this notebook

**Authors:** Jessica Krick, Jaladh Singhal, Brigitta Sipőcz

**Updated:** 2026-04-22

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions
or problems.

**Runtime:** As of the date above, this notebook takes about 2 minutes to run to completion on a machine with 8GB RAM and 4 CPU.

**AI Acknowledgement:**

This tutorial was developed with the assistance of AI tools

**References:**

- [Robitaille et al., 2013](https://www.aanda.org/articles/aa/full_html/2013/10/aa22068-13/aa22068-13.html)

- [Astropy Collaboration et al., 2018](https://arxiv.org/abs/1801.02634)

- [Astropy Collaboration et al., 2022](https://arxiv.org/abs/2206.14220)

- [OpenUniverse et al., 2025](https://arxiv.org/abs/2501.05632)
