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

# Euclid Q1: PHZ catalogs

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- Understand the basic characteristics of Euclid Q1 photo-z catalog and how to match it with MER mosaics.
- Understand what PHZ catalogs are available and how to view the columns in those catalogs.
- How to query with ADQL in the PHZ catalog to find galaxies between a redshift of 1.4 and 1.6.
- Pull and plot a spectrum of one of the galaxies in that catalog.
- Cutout an image of the galaxy to view it close up.
- Learn how to upload images and catalogs to Firefly to inspect individual sources in greater detail.

+++

## Introduction

+++

Euclid launched in July 2023 as a European Space Agency (ESA) mission with involvement by NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission is providing space-based imaging and spectroscopy as well as supporting ground-based imaging to achieve these primary goals.
These data will be archived by multiple global repositories, including IRSA, where they will support transformational work in many areas of astrophysics.

Euclid Quick Release 1 (Q1) consists of consists of ~30 TB of imaging, spectroscopy, and catalogs covering four non-contiguous fields:
Euclid Deep Field North (22.9 sq deg), Euclid Deep Field Fornax (12.1 sq deg), Euclid Deep Field South (28.1 sq deg), and LDN1641.

Among the data products included in the Q1 release are multiple catalogs created by the PHZ Processing Function.
This notebook provides an introduction to the main PHZ catalog, which contains 61 columns describing the photometric redshift probability distribution, fluxes, and classification for each source.
If you have questions about this notebook, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## Imports

```{important}
We rely on astroquery features that have been recently added, so please make sure you have version v0.4.10 or newer installed.
```

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install matplotlib 'astropy>=5.3' 'astroquery>=0.4.10' fsspec firefly_client
```

```{code-cell} ipython3
import os
import re
import urllib

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import QTable
from astropy import units as u
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch, LogStretch, quantity_support
from astropy.wcs import WCS

from firefly_client import FireflyClient
from astroquery.ipac.irsa import Irsa
```

## 1. Find the MER Tile ID that corresponds to a given RA and Dec

In this case, choose random coordinates to show a different MER mosaic image. Search a radius around these coordinates.

```{code-cell} ipython3
ra = 268
dec = 66
search_radius= 10 * u.arcsec

coord = SkyCoord(ra, dec, unit='deg', frame='icrs')
```

### Use IRSA to search for all Euclid data on this target

This searches specifically in the euclid_DpdMerBksMosaic "collection" which is the MER images and catalogs.

+++

```{note}
This table lists all MER mosaic images available in this position. These mosaics include the Euclid VIS, Y, J, H images, as well as ground-based telescopes which have been put on the same pixel scale. For more information, see the [Euclid documentation at IPAC](https://euclid.caltech.edu/page/euclid-faq-tech/). We use the ``facility`` argument below to query for Euclid images only.
```

```{code-cell} ipython3
image_table = Irsa.query_sia(pos=(coord, search_radius), collection='euclid_DpdMerBksMosaic', facility='Euclid')
```

Note that there are various image types are returned as well, we filter out the science images from these:

```{code-cell} ipython3
science_images = image_table[image_table['dataproduct_subtype'] == 'science']
science_images
```

Choose the VIS image and pull the filename and tileID

```{code-cell} ipython3
filename = science_images[science_images['energy_bandpassname'] == 'VIS']['access_url'][0]
tileID = science_images[science_images['energy_bandpassname'] == 'VIS']['obs_id'][0][:9]

print(f'The MER tile ID for this object is : {tileID}')
```

## 2. Download PHZ catalog from IRSA

Use IRSA's TAP to search catalogs

```{code-cell} ipython3
Irsa.list_catalogs(filter='euclid')
```

```{code-cell} ipython3
table_mer = 'euclid_q1_mer_catalogue'
table_phz = 'euclid_q1_phz_photo_z'
table_1dspectra = 'euclid.objectid_spectrafile_association_q1'
```

### Learn some information about the photo-z catalog:

- How many columns are there?
- List the column names

```{code-cell} ipython3
columns_info = Irsa.list_columns(catalog=table_phz)
print(len(columns_info))
```

```{tip}
The PHZ catalog contains 67 columns, below are a few highlights:

- object_id
- flux_vis_unif, flux_y_unif, flux_j_unif, flux_h_unif
- median redshift (phz_median)
- phz_classification
- phz_90_int1,  phz_90_int2 (The phz PDF interval containing 90% of the probability, upper and lower values)
```

```{code-cell} ipython3
# Full list of columns and their description
columns_info
```

```{note}
The phz_catalog on IRSA has more columns than it does on the ESA archive.
This is because the ESA catalog stores some information in one column (for example, phz_90_int is stored as [lower, upper], rather than in two separate columns).

The fluxes are different from the fluxes derived in the MER catalog.
The _unif fluxes are: "Unified flux recomputed after correction from galactic extinction and filter shifts".
```

+++

### Find some galaxies between 1.4 and 1.6 at a selected RA and Dec

We specify the following conditions on our search:
- We select just the galaxies where the flux is greater than zero, to ensure the appear in all four of the Euclid MER images.
- Select only objects in a circle (search radius selected below) around our selected RA and Dec
- `phz_classification = 2` means we select only galaxies
- Using the `phz_90_int1` and `phz_90_int2`, we select just the galaxies where the error on the photometric redshift is less than 20%
- Select just the galaxies between a median redshift of 1.4 and 1.6
- We search just a 5 arcminute box around an RA and Dec

+++

Search based on ``tileID``:

```{code-cell} ipython3
######################## User defined section ############################
## How large do you want the image cutout to be?
im_cutout= 5 * u.arcmin

## What is the center of the cutout?
ra_cutout = 267.8
dec_cutout =  66

coords_cutout = SkyCoord(ra_cutout, dec_cutout, unit='deg', frame='icrs')
size_cutout = im_cutout.to(u.deg).value
```

```{code-cell} ipython3
adql = ("SELECT DISTINCT mer.object_id, mer.ra, mer.dec, "
        "phz.flux_vis_unif, phz.flux_y_unif, phz.flux_j_unif, phz.flux_h_unif, "
        "phz.phz_classification, phz.phz_median, phz.phz_90_int1, phz.phz_90_int2 "
        f"FROM {table_mer} AS mer "
        f"JOIN {table_phz} as phz "
        "ON mer.object_id = phz.object_id "
        "WHERE 1 = CONTAINS(POINT('ICRS', mer.ra, mer.dec), "
                            f"BOX('ICRS', {ra_cutout}, {dec_cutout}, {size_cutout/np.cos(coords_cutout.dec)}, {size_cutout})) "
        "AND phz.flux_vis_unif> 0 "
        "AND  phz.flux_y_unif > 0 "
        "AND phz.flux_j_unif > 0 "
        "AND phz.flux_h_unif > 0 "
        "AND phz.phz_classification = 2 "
        "AND ((phz.phz_90_int2 - phz.phz_90_int1) / (1 + phz.phz_median)) < 0.20 "
        "AND phz.phz_median BETWEEN 1.4 AND 1.6")


## Use TAP with this ADQL string
result_galaxies = Irsa.query_tap(adql).to_table()
result_galaxies[:5]
```

```{warning}
Note that we use `to_table` above rather than `to_qtable`. While astropy's `QTable` is more powerful than its `Table`, as it e.g. handles the column units properly, we cannot use it here due to a known bug; it mishandles the large integer numbers in the `object_id` column and recast them as float during which process some precision is being lost.

Once the bug is fixed, we plan to update the code in this notebook and simplify some of the approaches below.
```

+++

## 3. Read in a cutout of the MER image from IRSA directly

+++

Due to the large field of view of the MER mosaic, let's cut out a smaller section (5'x5') of the MER mosaic to inspect the image.

```{code-cell} ipython3
## Use fsspec to interact with the fits file without downloading the full file
hdu = fits.open(filename, use_fsspec=True)

## Store the header
header = hdu[0].header

## Read in the cutout of the image that you want
cutout_image = Cutout2D(hdu[0].section, position=coords_cutout, size=im_cutout, wcs=WCS(header))
```

```{code-cell} ipython3
cutout_image.data.shape
```

```{code-cell} ipython3
norm = ImageNormalize(cutout_image.data, interval=PercentileInterval(99.9), stretch=AsinhStretch())
_ = plt.imshow(cutout_image.data, cmap='gray', origin='lower', norm=norm)
```

## 4. Overplot the catalog on the MER mosaic image

+++

```{tip}
We can rely on astropy's WCSAxes framework for making plots of Astronomical data in Matplotlib. Please note the usage of `projection` and `transform` arguments in the code example below.

For more info, please visit the [WCSAxes documentation](https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html).
```

```{code-cell} ipython3
ax = plt.subplot(projection=cutout_image.wcs)

ax.imshow(cutout_image.data, cmap='gray', origin='lower',
          norm=ImageNormalize(cutout_image.data, interval=PercentileInterval(99.9), stretch=LogStretch()))
plt.scatter(result_galaxies['ra'], result_galaxies['dec'], s=36, facecolors='none', edgecolors='red',
            transform=ax.get_transform('world'))

_ = plt.title('Galaxies between z = 1.4 and 1.6')
```

## 5. Pull the spectra on the top brightest source based on object ID

```{code-cell} ipython3
result_galaxies.sort(keys='flux_h_unif', reverse=True)
```

```{code-cell} ipython3
result_galaxies[:3]
```

Let's pick one of these galaxies. Note that the table has been sorted above, we can use the same index here and below to access the data for this particular galaxy.

```{code-cell} ipython3
index = 2

obj_id = result_galaxies['object_id'][index]
redshift = result_galaxies['phz_median'][index]
```

We will use TAP and an ASQL query to find the spectral data for this particular galaxy.

```{code-cell} ipython3
adql_object = f"SELECT * FROM {table_1dspectra} WHERE objectid = {obj_id}"

## Pull the data on this particular galaxy
result_spectra = Irsa.query_tap(adql_object).to_table()
result_spectra
```

Pull out the file name from the ``result_spectra`` table:

```{code-cell} ipython3
file_uri = urllib.parse.urljoin(Irsa.tap_url, result_spectra['uri'][0])
file_uri
```

```{code-cell} ipython3
with fits.open(file_uri) as hdul:
    spectrum = QTable.read(hdul[result_spectra['hdu'][0]], format='fits')
    spectrum_header = hdul[result_spectra['hdu'][0]].header
```

### Now the data are read in, plot the spectrum

```{tip}
As we use astropy.visualization’s quantity_support, matplotlib automatically picks up the axis units from the quantitites we plot.
```

```{code-cell} ipython3
quantity_support()
```

```{code-cell} ipython3
plt.plot(spectrum['WAVELENGTH'].to(u.micron), spectrum['SIGNAL'])

plt.xlim(1.25, 1.85)
plt.ylim(-0.5, 0.5)
_ = plt.title(f"Object {obj_id} with phz_median={redshift}")
```

Let's cut out a very small patch of the MER image to see what this galaxy looks like. Remember that we sorted the table above, so can reuse the same index to pick up the coordinates for the galaxy. Otherwise we could filter on the object ID.

```{code-cell} ipython3
result_galaxies[index]
```

```{code-cell} ipython3
## How large do you want the image cutout to be?
size_galaxy_cutout = 2.0 * u.arcsec
```

Use the `ra` and `dec` columns for the galaxy to create a `SkyCoord`.

```{code-cell} ipython3
coords_galaxy = SkyCoord(result_galaxies['ra'][index], result_galaxies['dec'][index], unit='deg')
```

```{code-cell} ipython3
coords_galaxy
```

We haven't closed the image file above, so use `Cutout2D` again to cut out a section around the galaxy.

```{code-cell} ipython3
cutout_galaxy = Cutout2D(hdu[0].section, position=coords_galaxy, size=size_galaxy_cutout, wcs=WCS(header))
```

Plot to show the cutout on the galaxy

```{code-cell} ipython3
ax = plt.subplot(projection=cutout_galaxy.wcs)

ax.imshow(cutout_galaxy.data, cmap='gray', origin='lower',
          norm=ImageNormalize(cutout_galaxy.data, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
```

## 6. Load the image on Firefly to be able to interact with the data directly

+++

Save the data locally if you have not already done so, in order to upload to IRSA viewer.

```{code-cell} ipython3
download_path = "data"
if os.path.exists(download_path):
    print("Output directory already created.")
else:
    print("Creating data directory.")
    os.mkdir(download_path)
```

### Vizualize the image with Firefly

First initialize the client, then set the path to the image, upload it to firefly, load it and align with WCS.

Note this can take a while to upload the full MER image.

```{code-cell} ipython3
fc = FireflyClient.make_client('https://irsa.ipac.caltech.edu/irsaviewer')

fc.show_fits(url=filename)

fc.align_images(lock_match=True)
```

### Save the table as a CSV for Firefly upload

```{code-cell} ipython3
csv_path = os.path.join(download_path, "mer_df.csv")
result_galaxies.write(csv_path, format="csv")
```

### Upload the CSV table to Firefly and display as an overlay on the FITS image

```{code-cell} ipython3
uploaded_table = fc.upload_file(csv_path)
print(f"Uploaded Table URL: {uploaded_table}")

fc.show_table(uploaded_table)
```

## About this Notebook

**Author**: Tiffany Meshkat, Anahita Alavi, Anastasia Laity, Andreas Faisst, Brigitta Sipőcz, Dan Masters, Harry Teplitz, Jaladh Singhal, Shoubaneh Hemmati, Vandana Desai

**Updated**: 2025-04-10

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
