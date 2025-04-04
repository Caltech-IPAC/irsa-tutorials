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

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install matplotlib pandas 'astropy>=5.3' 'astroquery>=0.4.10' fsspec firefly_client
```

```{code-cell} ipython3
import os
import re
import urllib

import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy import units as u
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch, LogStretch
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

```{code-cell} ipython3
image_table = Irsa.query_sia(pos=(coord, search_radius), collection='euclid_DpdMerBksMosaic')
```

```{note}
This table lists all MER mosaic images available in this position. These mosaics include the Euclid VIS, Y, J, H images, as well as ground-based telescopes which have been put on the same pixel scale. For more information, see the [Euclid documentation at IPAC](https://euclid.caltech.edu/page/euclid-faq-tech/).
```

```{code-cell} ipython3
# Convert the table to pandas dataframe
df_im_irsa=image_table.to_pandas()
```

```{code-cell} ipython3
df_im_euclid=df_im_irsa[ (df_im_irsa['dataproduct_subtype']=='science') &  (df_im_irsa['facility_name']=='Euclid')]

df_im_euclid.head()
```

Choose the VIS image and pull the filename and tileID

```{code-cell} ipython3
filename=df_im_euclid[df_im_euclid['energy_bandpassname']=='VIS']['access_url'].to_list()[0]

tileID=re.search(r'TILE\s*(\d{9})', filename).group(1)
print('The MER tile ID for this object is :',tileID)
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
- phz_classification =2 means we select only galaxies
- Using the phz_90_int1 and phz_90_int2, we select just the galaxies where the error on the photometric redshift is less than 20%
- Select just the galaxies between a median redshift of 1.4 and 1.6
- We search just a 3 arcminute box around an RA and Dec

+++

Search based on ``tileID``:

```{code-cell} ipython3
######################## User defined section ############################
## How large do you want the image cutout to be?
im_cutout= 5 * u.arcmin

## What is the center of the cutout?
ra_cutout = 267.8
dec_cutout =  66

coords_cutout = SkyCoord(ra_cutout, dec_cutout, unit=(u.deg, u.deg), frame='icrs')
##########################################################################
im_cutout_deg=im_cutout.to(u.deg).value

## Create ra and dec bounds for the box
ra0=ra_cutout-im_cutout_deg/2
ra1=ra_cutout+im_cutout_deg/2

dec0=dec_cutout-im_cutout_deg/2
dec1=dec_cutout+im_cutout_deg/2
```

```{code-cell} ipython3
adql = f"SELECT DISTINCT mer.object_id,mer.ra, mer.dec, phz.flux_vis_unif, phz.flux_y_unif, \
phz.flux_j_unif, phz.flux_h_unif, phz.phz_classification, phz.phz_median, phz.phz_90_int1, phz.phz_90_int2 \
FROM {table_mer} AS mer \
JOIN {table_phz} as phz \
ON mer.object_id = phz.object_id \
WHERE 1 = CONTAINS(POINT('ICRS', mer.ra, mer.dec), CIRCLE('ICRS', {ra_cutout}, {dec_cutout}, {im_cutout_deg/2})) \
AND phz.flux_vis_unif> 0 \
AND  phz.flux_y_unif > 0 \
AND phz.flux_j_unif > 0 \
AND phz.flux_h_unif > 0 \
AND phz.phz_classification = 2 \
AND ((phz.phz_90_int2 - phz.phz_90_int1) / (1 + phz.phz_median)) < 0.20 \
AND phz.phz_median BETWEEN 1.4 AND 1.6 \
"
adql


## Use TAP with this ADQL string
result = Irsa.query_tap(adql)


## Convert table to pandas dataframe
df_g_irsa = result.to_table().to_pandas()

# Display first few rows
df_g_irsa.head()
```

```{code-cell} ipython3
len(df_g_irsa)
```

## 3. Read in a cutout of the MER image from IRSA directly

```{code-cell} ipython3
## Use fsspec to interact with the fits file without downloading the full file
hdu = fits.open(filename, use_fsspec=True)

## Store the header
header = hdu[0].header

## Read in the cutout of the image that you want
cutout_data = Cutout2D(hdu[0].section, position=coords_cutout, size=im_cutout, wcs=WCS(hdu[0].header))

## Define a new fits file based on this smaller cutout, with accurate WCS based on the cutout size
new_hdu = fits.PrimaryHDU(data=cutout_data.data, header=header)
new_hdu.header.update(cutout_data.wcs.to_header())
```

```{code-cell} ipython3
im_mer_irsa = new_hdu.data
```

```{code-cell} ipython3
im_mer_irsa.shape
```

```{code-cell} ipython3
norm = ImageNormalize(im_mer_irsa, interval=PercentileInterval(99.9), stretch=AsinhStretch())
plt.imshow(im_mer_irsa, cmap='gray', origin='lower', norm=norm)
```

## 4. Overplot the catalog on the MER mosaic image

```{code-cell} ipython3
## Use the WCS package to extract the coordinates from the header of the image
head_mer_irsa=new_hdu.header
wcs_irsa=WCS(head_mer_irsa) # VIS
```

```{code-cell} ipython3
## Convert the catalog to match the pixels in the image
xy_irsa = wcs_irsa.all_world2pix(df_g_irsa["ra"],df_g_irsa["dec"],0)
```

```{code-cell} ipython3
df_g_irsa['x_pix']=xy_irsa[0]
df_g_irsa['y_pix']=xy_irsa[1]
```

Due to the large field of view of the MER mosaic, let's cut out a smaller section (3'x3')of the MER mosaic to inspect the image

+++

Plot MER catalog sources on the Euclid VIS image 3'x3' cutout

```{code-cell} ipython3
plt.imshow(im_mer_irsa, cmap='gray', origin='lower', norm=ImageNormalize(im_mer_irsa, interval=PercentileInterval(99.9), stretch=LogStretch()))
colorbar = plt.colorbar()
plt.scatter(df_g_irsa['x_pix'], df_g_irsa['y_pix'], s=36, facecolors='none', edgecolors='red')

plt.title('Galaxies between z = 1.4 and 1.6')
plt.show()
```

Pull the spectra on the top brightest source based on object ID

```{code-cell} ipython3
df_g_irsa_sort=df_g_irsa.sort_values(by='flux_h_unif',ascending=False)
```

```{code-cell} ipython3
df_g_irsa_sort.iloc[0:3]
```

```{code-cell} ipython3
obj_id=df_g_irsa_sort['object_id'].iloc[2]
redshift = df_g_irsa_sort['phz_median'].iloc[2]

## Pull the data on these objects
adql_object = f"SELECT * \
FROM {table_1dspectra} \
WHERE objectid = {obj_id}"

## Pull the data on this particular galaxy
result2 = Irsa.query_tap(adql_object)
df2=result2.to_table().to_pandas()
df2
```

Pull out the file name from the ``result`` table:

```{code-cell} ipython3
file_uri = urllib.parse.urljoin(Irsa.tap_url, result2['uri'][0])
file_uri
```

```{code-cell} ipython3
with fits.open(file_uri) as hdul:
    hdu = hdul[df2['hdu'].iloc[0]]
    dat = Table.read(hdu, format='fits', hdu=1)
    df_obj_irsa = dat.to_pandas()
```

### Now the data are read in, plot the spectrum

Divide by 10000 to convert from Angstrom to micron

```{code-cell} ipython3
plt.plot(df_obj_irsa['WAVELENGTH']/10000., df_obj_irsa['SIGNAL'])

plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux (erg / (s cm2))')
plt.xlim(1.25, 1.85)
plt.ylim(-0.5,0.5)
plt.title('Object ID is '+str(obj_id)+'with phz_median='+str(redshift))
```

Let's cut out a very small patch of the MER image to see what this galaxy looks like

```{code-cell} ipython3
## How large do you want the image cutout to be?
im_cutout= 2.0 * u.arcsec

## Use the ra and dec of the galaxy
ra = df_g_irsa[df_g_irsa['object_id']==obj_id]['ra'].iloc[0]
dec =  df_g_irsa[df_g_irsa['object_id']==obj_id]['dec'].iloc[0]

coords_cutout = SkyCoord(ra, dec, unit=(u.deg,u.deg), frame='icrs')
```

Use ``fsspec`` to obtain a cutout of the fits file

```{code-cell} ipython3
hdu = fits.open(filename, use_fsspec=True)

header = hdu[0].header
```

```{code-cell} ipython3
## Read in the cutout of the image that you want
cutout_data = Cutout2D(hdu[0].section, position=coords_cutout, size=im_cutout, wcs=WCS(hdu[0].header))
```

```{code-cell} ipython3
new_hdu = fits.PrimaryHDU(data=cutout_data.data, header=header)
new_hdu.header.update(cutout_data.wcs.to_header())
```

```{code-cell} ipython3
## Plot a quick simple plot to show the cutout on the galaxy

plt.imshow(new_hdu.data, cmap='gray', origin='lower',
           norm=ImageNormalize(new_hdu.data, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
colorbar = plt.colorbar()
```

## 5. Load the image on Firefly to be able to interact with the data directly

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
df_g_irsa.to_csv(csv_path, index=False)
```

### Upload the CSV table to Firefly and display as an overlay on the FITS image

```{code-cell} ipython3
uploaded_table = fc.upload_file(csv_path)
print(f"Uploaded Table URL: {uploaded_table}")

fc.show_table(uploaded_table)
```

## About this Notebook

**Author**: Tiffany Meshkat, Anahita Alavi, Anastasia Laity, Andreas Faisst, Brigitta Sipőcz, Dan Masters, Harry Teplitz, Jaladh Singhal, Shoubaneh Hemmati, Vandana Desai

**Updated**: 2025-03-31

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
