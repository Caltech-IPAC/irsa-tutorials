---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.3
kernelspec:
  display_name: irsa-tutorials
  language: python
  name: python3
---

# Euclid Quick Release 1: Cloud Access

+++

## Learning Goals
- Learn where Euclid QR1 data is present in the cloud
- Find an image and retireve its cutout from cloud
- Find an object and retrieve its spectrum from cloud

+++

## Introduction
TODO: fill with links of docs/specs?

+++

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install s3fs astropy astroquery matploltib
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

## Browse Euclid QR1 Bucket

```{code-cell} ipython3
BUCKET_NAME='nasa-irsa-euclid-q1' # internal to IPAC until public release (use LAN or VPN w/ Tunnel-all)
```

```{code-cell} ipython3
s3 = s3fs.S3FileSystem(anon=True)
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1')
```

## Find images for a coordinate search

+++

### Locate MER images in the bucket

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/MER')
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/MER/102018211')
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/MER/102018211/VIS')
```

As per doc specification, we need `MER/{tile_id}/{instrument}/EUC_MER_BGSUB-MOSAIC*.fits` for displaying background-subtracted mosiac images. But these images are stored under TILE IDs so first we need to find TILE ID for a coordinate search we are interested in. We will use astroquery (in next section) to retrieve FITS file paths for our coordinates.

+++

### Get image file paths for a coordinate search of interest

```{code-cell} ipython3
coord = SkyCoord.from_name("TYC 4429-1677-1")
search_radius = 10 * u.arcsec
```

List all Simple Image Access (SIA) collections for IRSA.

```{code-cell} ipython3
collections = Irsa.list_collections(servicetype='SIA')
len(tbl)
```

Filter to only those containing "euclid":

```{code-cell} ipython3
tbl[['euclid' in v for v in tbl['collection']]]
```

We identify the collection we need for MER images:

```{code-cell} ipython3
img_collection = 'euclid_DpdMerBksMosaic'
```

```{code-cell} ipython3
img_tbl = Irsa.query_sia(pos=(coord, search_radius), collection=img_collection).to_table()
img_tbl
```

Now we narrow it down to the images with science dataproduct subtype and Euclid facility:

```{code-cell} ipython3
euclid_sci_img_tbl = img_tbl[[row['facility_name']=='Euclid' and row['dataproduct_subtype']=='science' for row in img_tbl]]
euclid_sci_img_tbl
```

We can see there's a `cloud_access` column that gives us the location info of the image files we are interested in. So let's extract the S3 bucket file path from it.

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

## Retrieve image cutouts from the cloud
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

## Find objects for the coordinates of our interest

+++

### Locate MER catalogs in the bucket

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/catalogs')
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/catalogs/MER_FINAL_CATALOG')
```

```{code-cell} ipython3
s3.ls(f'{BUCKET_NAME}/q1/catalogs/MER_FINAL_CATALOG/102018211')
```

As per doc specification, we need `catalogs/MER_FINAL_CATALOG/{tile_id}/EUC_MER_FINAL-CAT*.fits` for listing the objects catalogued. But we only need to find objects in our coordinates of interest so we will use astroquery to do a spatial search in MER catalog (combined for all tiles).

+++

### Get object IDs for the coordinates of our interest

```{code-cell} ipython3
tbl_catalogs = Irsa.list_catalogs(full=True).to_table()
len(tbl_catalogs)
```

```{code-cell} ipython3
tbl_catalogs[['euclid' in v for v in tbl_catalogs['schema_name']]]
```

From this table, we can extract the MER catalog name. We also see several other interesting catalogs, let's also extract spectral file association catalog for retrieving spectra later.

```{code-cell} ipython3
euclid_mer_catalog = 'euclid_q1_mer_catalogue'
euclid_spec_association_catalog = 'euclid.objectid_spectrafile_association_q1'
```

Now, we do a TAP search with spatial constraints for our coordinates. We use cone of 5 arcsec around our source to pinpoint its object ID in Euclid catalog.

```{code-cell} ipython3
search_radius = (5 * u.arcsec).to('deg')

adql_query = f"SELECT * \
    FROM {euclid_mer_catalog} \
    WHERE CONTAINS(POINT('ICRS', ra, dec), \
        CIRCLE('ICRS', {ra.value}, {dec.value}, {search_radius.value})) = 1"

mer_catalog_tbl = Irsa.query_tap(query=adql_query).to_table()
mer_catalog_tbl
```

```{code-cell} ipython3
object_id = int(mer_catalog_tbl['object_id'][0])
object_id
```

## Find spectra for the coordinates of our interest
Using the object ID(s) we extracted above, we can narrow down the spectral file association catalog to identify spectra file path(s).

```{code-cell} ipython3
adql_query = f"SELECT * FROM {euclid_spec_association_catalog} \
    WHERE objectid = {object_id}"

spec_association_tbl = Irsa.query_tap(adql_query).to_table()
spec_association_tbl
```

We can see the `uri` column that gives us location of spectra file on IBE, we can map it to S3 bucket key to retrieve spectra file from the cloud. This is a very big FITS spectra file with multiple extensions where each extension contains spectrum of one object. The `hdu` column gives us the extension number for our object. So let's extract both of these.

```{code-cell} ipython3
spec_fpath_key = spec_association_tbl['uri'][0].replace('ibe/data/euclid/', '')
spec_fpath_key
```

```{code-cell} ipython3
object_hdu_idx = int(spec_association_tbl['hdu'][0])
object_hdu_idx
```

## Retrieve spectrum from the cloud
Again we use astropy's lazy-loading capability of FITS to only retrieve the spectrum table of our object from the S3 bucket.

```{code-cell} ipython3
with fits.open(f's3://{BUCKET_NAME}/{spec_fpath_key}', fsspec_kwargs={'anon': True}) as hdul:
    spec_hdu = hdul[object_hdu_idx]
    spec_tbl = Table.read(spec_hdu)
```

```{code-cell} ipython3
spec_tbl
```

```{code-cell} ipython3
plt.plot(spec_tbl['WAVELENGTH'], spec_tbl['SIGNAL'])
plt.xlabel(spec_tbl['WAVELENGTH'].unit.to_string('latex_inline'))
plt.ylabel(spec_tbl['SIGNAL'].unit.to_string('latex_inline'))

plt.title(f'Euclid Object ID: {object_id}');
```

## About this Notebook

**Author:** Jaladh Singhal (IRSA Developer) in conjunction with Tiffany Meshkat, Vandana Desai, Brigitta Sipőcz, and the IPAC Science Platform team

**Updated:** 2025-03-13

**Contact:** the [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
