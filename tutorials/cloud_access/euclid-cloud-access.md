---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Euclid Q1: cloud access

+++

## Learning Goals
By the end of this tutorial, you will:
- Learn where Euclid Q1 data are stored in the cloud.
- Retrieve an image cutout from the cloud.
- Retrieve a spectrum from the cloud.

+++

## 1. Introduction

Euclid launched in July 2023 as a European Space Agency (ESA) mission with involvement by NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission is providing space-based imaging and spectroscopy as well as supporting ground-based imaging to achieve these primary goals.
These data will be archived by multiple global repositories, including IRSA, where they will support transformational work in many areas of astrophysics.

Euclid Quick Release 1 (Q1) consists of consists of ~30 TB of imaging, spectroscopy, and catalogs covering four non-contiguous fields:
Euclid Deep Field North (22.9 sq deg), Euclid Deep Field Fornax (12.1 sq deg), Euclid Deep Field South (28.1 sq deg), and LDN1641.

IRSA maintains copies of the Euclid Q1 data products both on premises at IPAC and on the cloud via Amazon Web Services (AWS).
This notebook provides an introduction to accessing Euclid Q1 data from the cloud.
If you have questions, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## 2. Imports
- `s3fs` for browsing S3 buckets
- `astropy` for handling coordinates, units, FITS I/O, tables, images, etc.
- `astroquery>=0.4.10` for querying Euclid data products from IRSA
- `matplotlib` for visualization
- `json` for decoding JSON strings

```{important}
We rely on ``astroquery`` features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install s3fs astropy 'astroquery>=0.4.10' matplotlib
```

```{code-cell} ipython3
import s3fs
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.table import Table
from astroquery.ipac.irsa import Irsa
from matplotlib import pyplot as plt
import json
```

## 3. Browse Euclid Q1 cloud-hosted data

```{code-cell} ipython3
BUCKET_NAME = 'nasa-irsa-euclid-q1'
```

[s3fs](https://s3fs.readthedocs.io/en/latest/) provides a filesystem-like python interface for AWS S3 buckets. First we create a s3 client:

```{code-cell} ipython3
s3 = s3fs.S3FileSystem(anon=True)
```

Then we list the `q1` directory that contains Euclid Q1 data products:

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1')
```

Let's navigate to MER images (available as FITS files):

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/MER')[:10] # ls only top 10 to limit the long output
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/MER/102018211') # pick any tile ID from above
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/MER/102018211/VIS') # pick any instrument from above
```

As per "Browsable Directories" section in [user guide](https://irsa.ipac.caltech.edu/data/Euclid/docs/euclid_archive_at_irsa_user_guide.pdf), we need `MER/{tile_id}/{instrument}/EUC_MER_BGSUB-MOSAIC*.fits` for displaying background-subtracted mosiac images. But these images are stored under TILE IDs so first we need to find TILE ID for a coordinate search we are interested in. We will use astroquery (in next section) to retrieve FITS file paths for our coordinates by doing spatial search.

+++

## 4. Do a spatial search for MER mosaics

Pick a target and search radius:

```{code-cell} ipython3
target_name = 'TYC 4429-1677-1'
coord = SkyCoord.from_name(target_name)
search_radius = 10 * u.arcsec
```

% List all Simple Image Access (SIA) collections for IRSA with names containing "euclid":
%
% ```{code-cell} ipython3
% collections = Irsa.list_collections(servicetype='SIA', filter='euclid')
% collections
% ```

As per "Data Products Overview" in [user guide](https://irsa.ipac.caltech.edu/data/Euclid/docs/euclid_archive_at_irsa_user_guide.pdf) and above table, we identify that MER Mosiacs are available as the following collection:

```{code-cell} ipython3
img_collection = 'euclid_DpdMerBksMosaic'
```

Now query this collection for our target's coordinates and search radius:

```{code-cell} ipython3
img_tbl = Irsa.query_sia(pos=(coord, search_radius), collection=img_collection)
img_tbl
```

Let's narrow it down to the images with science dataproduct subtype and Euclid facility:

```{code-cell} ipython3
euclid_sci_img_tbl = img_tbl[[row['facility_name']=='Euclid'
                              and row['dataproduct_subtype']=='science'
                              for row in img_tbl]]
euclid_sci_img_tbl
```

We can see there's a `cloud_access` column that gives us the location info of the image files we are interested in. So let's extract the S3 bucket file path from it:

```{code-cell} ipython3
def get_s3_fpath(cloud_access):
    cloud_info = json.loads(cloud_access) # converts str to dict
    bucket_name = cloud_info['aws']['bucket_name']
    key = cloud_info['aws']['key']

    return f'{bucket_name}/{key}'
```

```{code-cell} ipython3
[get_s3_fpath(row['cloud_access']) for row in euclid_sci_img_tbl]
```

Let's also extract filter names to use when displaying the images:

```{code-cell} ipython3
def get_filter_name(instrument, bandpass):
    return f'{instrument}_{bandpass}' if instrument!=bandpass else instrument
```

```{code-cell} ipython3
[get_filter_name(row['instrument_name'], row['energy_bandpassname']) for row in euclid_sci_img_tbl]
```

## 5. Efficiently retrieve mosaic cutouts
These image files are very big (~1.4GB), so we use astropy's lazy-loading capability of FITS for better performance. (See [Obtaining subsets from cloud-hosted FITS files](https://docs.astropy.org/en/stable/io/fits/usage/cloud.html#fits-io-cloud).)

```{code-cell} ipython3
cutout_size = 1 * u.arcmin
```

```{code-cell} ipython3
cutouts = []
filters = []

for row in euclid_sci_img_tbl:
    s3_fpath = get_s3_fpath(row['cloud_access'])
    filter_name = get_filter_name(row['instrument_name'], row['energy_bandpassname'])

    with fits.open(f's3://{s3_fpath}', fsspec_kwargs={"anon": True}) as hdul:
        print(f'Retrieving cutout for {filter_name} ...')
        cutout = Cutout2D(hdul[0].section,
                          position=coord,
                          size=cutout_size,
                          wcs=WCS(hdul[0].header))
        cutouts.append(cutout)
        filters.append(filter_name)
```

```{code-cell} ipython3
fig, axes = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2), subplot_kw={'projection': cutouts[0].wcs})

for idx, ax in enumerate(axes.flat):
    norm = ImageNormalize(cutouts[idx].data, interval=PercentileInterval(99), stretch=AsinhStretch())
    ax.imshow(cutouts[idx].data, cmap='gray', origin='lower', norm=norm)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.text(0.95, 0.05, filters[idx], color='white', fontsize=14, transform=ax.transAxes, va='bottom', ha='right')

plt.tight_layout()
```

## 6. Find the MER catalog for a given tile
Let's navigate to MER catalog in the Euclid Q1 bucket:

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/catalogs')
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/catalogs/MER_FINAL_CATALOG')[:10] # ls only top 10 to limit the long output
```

```{code-cell} ipython3
mer_tile_id = 102160339 # from the image paths for the target we picked
s3.ls(f'{BUCKET_NAME}/q1/catalogs/MER_FINAL_CATALOG/{mer_tile_id}')
```

As per "Browsable Directiories" section in [user guide](https://irsa.ipac.caltech.edu/data/Euclid/docs/euclid_archive_at_irsa_user_guide.pdf), we can use `catalogs/MER_FINAL_CATALOG/{tile_id}/EUC_MER_FINAL-CAT*.fits` for listing the objects catalogued. We can read the identified FITS file as table and do filtering on ra, dec columns to find object ID(s) only for the target we picked. But it will be an expensive operation so we will instead use astroquery (in next section) to do a spatial search in the MER catalog provided by IRSA.

```{note}
Once the catalogs are available as Parquet files in the cloud, we can efficiently do spatial filtering directly on the cloud-hosted file to identify object ID(s) for our target. But for the time being, we can use catalog VO services through astroquery to do the same.
```

+++

## 7. Find the MER Object ID for our target
First, list the Euclid catalogs provided by IRSA:

```{code-cell} ipython3
catalogs = Irsa.list_catalogs(full=True, filter='euclid')
catalogs
```

From this table, we can extract the MER catalog name. We also see several other interesting catalogs, let's also extract spectral file association catalog for retrieving spectra later.

```{code-cell} ipython3
euclid_mer_catalog = 'euclid_q1_mer_catalogue'
euclid_spec_association_catalog = 'euclid.objectid_spectrafile_association_q1'
```

Now, we do a region search within a cone of 5 arcsec around our target to pinpoint its object ID in Euclid catalog:

```{code-cell} ipython3
search_radius = 5 * u.arcsec

mer_catalog_tbl = Irsa.query_region(coordinates=coord, spatial='Cone',
                                    catalog=euclid_mer_catalog, radius=search_radius)
mer_catalog_tbl
```

```{code-cell} ipython3
object_id = int(mer_catalog_tbl['object_id'][0])
object_id
```

## 8. Find the spectrum of an object in the MER catalog
Using the object ID(s) we extracted above, we can narrow down the spectral file association catalog to identify spectra file path(s). So we do the following TAP search:

```{code-cell} ipython3
adql_query = f"SELECT * FROM {euclid_spec_association_catalog} \
    WHERE objectid = {object_id}"

spec_association_tbl = Irsa.query_tap(adql_query).to_table()
spec_association_tbl
```

```{warning}
If you picked a target other than what this notebook uses, it's possible that there is no spectrum associated for your target's object ID. In that case, `spec_association_tbl` will contain 0 rows.
```

In above table, we can see that the `path` column gives us a url that can be used to call the SpectrumDM service to get the spectrum of our object. We can map it to an S3 bucket key to retrieve a spectra file from the cloud. This is a very big FITS spectra file with multiple extensions where each extension contains spectrum of one object. The `hdu` column gives us the extension number for our object. So let's extract both of these.

```{code-cell} ipython3
spec_fpath_key = spec_association_tbl['path'][0].replace('api/spectrumdm/convert/euclid/', '').split('?')[0]
spec_fpath_key
```

```{code-cell} ipython3
object_hdu_idx = int(spec_association_tbl['hdu'][0])
object_hdu_idx
```

Again, we use astropy's lazy-loading capability of FITS to only retrieve the spectrum table of our object from the S3 bucket.

```{code-cell} ipython3
with fits.open(f's3://{BUCKET_NAME}/{spec_fpath_key}', fsspec_kwargs={'anon': True}) as hdul:
    spec_hdu = hdul[object_hdu_idx]
    spec_tbl = Table.read(spec_hdu)
    spec_header = spec_hdu.header
```

```{code-cell} ipython3
spec_tbl
```

```{code-cell} ipython3
# The signal needs to be multiplied by the scale factor in the header.
plt.plot(spec_tbl['WAVELENGTH'], spec_header['FSCALE'] * spec_tbl['SIGNAL'])
plt.xlabel(spec_tbl['WAVELENGTH'].unit.to_string('latex_inline'))
plt.ylabel(spec_tbl['SIGNAL'].unit.to_string('latex_inline'))

plt.title(f'Spectrum of Target: {target_name}\n(Euclid Object ID: {object_id})');
```

## About this Notebook

**Author:** Jaladh Singhal (IRSA Developer) in conjunction with Vandana Desai, Brigitta Sipőcz, Tiffany Meshkat, Troy Raen, and the IRSA Data Science Team

**Updated:** 2025-09-23

**Contact:** the [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
