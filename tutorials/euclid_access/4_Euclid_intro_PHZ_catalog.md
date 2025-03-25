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

# Euclid Quick Release 1: PHZ catalog matching to MER images

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

Euclid is a European Space Agency (ESA) space mission with NASA participation, to study the geometry and nature of the dark Universe. The Quick Data Release 1 (Q1) are the first data release from the Euclid mission after the Early Release Observations (ERO). On March 19, 2025 the data will be available on the ESA archive (https://easidr.esac.esa.int/sas/) and on the IRSA archive (https://irsa.ipac.caltech.edu).

These notebooks focus on how to access, download, and process Euclid Q1 data from the IRSA archive. At the end of the notebook, we also include some information for how to access the Q1 data from the ESA archive. If you have any issues accessing data from the archives, please contact the helpdesk directly: IRSA (irsasupport@ipac.caltech.edu) and ESA (https://support.cosmos.esa.int/euclid).

The photometry of every source is processed through a photometric redshift fitting pipeline, producing several different catalogs. This notebook provides an introduction to photo-z catalog released as part of Euclid Q1. Other Euclid notebooks show how to use other data products released as part of Euclid Q1.

+++

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install requests matplotlib pandas astropy pyvo fsspec firefly_client
```

```{code-cell} ipython3
from io import BytesIO
import os
import re

import requests
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
import pyvo as vo
```

# Introduction to Euclid Q1 PHZ catalog

+++

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
irsa_service= vo.dal.sia2.SIA2Service('https://irsa.ipac.caltech.edu/SIA')

image_table = irsa_service.search(pos=(coord, search_radius), collection='euclid_DpdMerBksMosaic').to_table()
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

## Choose the VIS image and pull the filename and tileID

```{code-cell} ipython3
filename=df_im_euclid[df_im_euclid['energy_bandpassname']=='VIS']['access_url'].to_list()[0]

tileID=re.search(r'TILE\s*(\d{9})', filename).group(1)
print('The MER tile ID for this object is :',tileID)
```

## 2. Download PHZ catalog from IRSA directly to this notebook

```{code-cell} ipython3
## Use IRSA to search for catalogs

service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")


## Search for all tables in IRSA labled as euclid_q1
tables = service.tables
for tablename in tables.keys():
    if "tap_schema" not in tablename and "euclid_q1" in tablename:
            tables[tablename].describe()
```

```{code-cell} ipython3
table_mer= 'euclid_q1_mer_catalogue'
table_phz= 'euclid_q1_phz_photo_z'
table_1dspectra= 'euclid.objectid_spectrafile_association_q1'
```

### Learn some information about the table:
- How many columns are there?
- List the column names

```{code-cell} ipython3
columns = tables[table_phz].columns
print(len(columns))
```

```{code-cell} ipython3
for col in columns:
    print(f'{f"{col.name}":30s}  {col.unit}  {col.description}') ## Currently no descriptions
```

## Note that the phz catalog contains 67 columns, below are a few highlights:

- object_id
- flux_vis_unif, flux_y_unif, flux_j_unif, flux_h_unif
- median redshift (phz_median)
- phz_classification
- phz_90_int1,  phz_90_int2 (The phz PDF interval containing 90% of the probability, upper and lower values)

We note that the phz_catalog on IRSA has more columns than it does on the ESA archive. This is because the ESA catalog stores some information in one column (for example, phz_90_int is stored as [lower, upper], rather than in two separate columns)

The fluxes are different from the fluxes derived in the MER catalog. The _unif fluxes are: "Unified flux recomputed after correction from galactic extinction and filter shifts"

+++

## Find some galaxies between 1.4 and 1.6 at a selected RA and Dec

We specify the following conditions on our search:
- We select just the galaxies where the flux is greater than zero, to ensure the appear in all four of the Euclid MER images.
- Select only objects in a circle (search radius selected below) around our selected RA and Dec
- phz_classification =2 means we select only galaxies
- Using the phz_90_int1 and phz_90_int2, we select just the galaxies where the error on the photometric redshift is less than 20%
- Select just the galaxies between a median redshift of 1.4 and 1.6

+++

### Search based on tileID

```{code-cell} ipython3
adql = f"SELECT DISTINCT mer.object_id,mer.ra, mer.dec, phz.flux_vis_unif, phz.flux_y_unif, \
phz.flux_j_unif, phz.flux_h_unif, phz.phz_classification, phz.phz_median, phz.phz_90_int1, phz.phz_90_int2 \
FROM {table_mer} AS mer \
JOIN {table_phz} as phz \
ON mer.object_id = phz.object_id \
WHERE phz.flux_vis_unif> 0 \
AND  phz.flux_y_unif > 0 \
AND phz.flux_j_unif > 0 \
AND phz.flux_h_unif > 0 \
AND phz.phz_classification = 2 \
AND mer.tileid = {tileID} \
AND ((phz.phz_90_int2 - phz.phz_90_int1) / (1 + phz.phz_median)) < 0.20 \
AND phz.phz_median BETWEEN 1.4 AND 1.6 \
"
adql


## Use TAP with this ADQL string using pyvo
result = service.search(adql)


## Convert table to pandas dataframe
df_g_irsa = result.to_table().to_pandas()

# Display first few rows
df_g_irsa.head()
```

## 3. Read in the MER image from IRSA directly

```{code-cell} ipython3
print(filename)
```

```{code-cell} ipython3
##Download the MER image -- note this file is about 1.46 GB

fname = download_file(filename, cache=True)
hdu_mer_irsa = fits.open(fname)
head_mer_irsa = hdu_mer_irsa[0].header

print(hdu_mer_irsa.info())
```

#### Now you've downloaded this large file, if you would like to save it to disk, uncomment the following cell

```{code-cell} ipython3
# download_path='/yourlocalpath/'
# hdu_mer_irsa.writeto(download_path+'./MER_image_VIS.fits', overwrite=True)
```

## 4. Overplot the catalog on the MER mosaic image

```{code-cell} ipython3
df_g_irsa.head()
```

```{code-cell} ipython3
## Use the WCS package to extract the coordinates from the header of the image
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

```{code-cell} ipython3
df_g_irsa
```

## Pull the spectra on the top brightest source based on object ID

```{code-cell} ipython3
df_g_irsa_sort=df_g_irsa.sort_values(by='flux_vis_unif',ascending=False)
```

```{code-cell} ipython3
df_g_irsa_sort[0:3]
```

```{code-cell} ipython3
obj_id=df_g_irsa_sort['object_id'].iloc[1]

## Pull the data on these objects
adql_object = f"SELECT * \
FROM {table_1dspectra} \
WHERE objectid = {obj_id} \
AND uri IS NOT NULL "

## Pull the data on this particular galaxy
result2 = service.search(adql_object)
df2=result2.to_table().to_pandas()
df2
```

```{code-cell} ipython3
## Create the full filename/url
irsa_url='https://irsa.ipac.caltech.edu/'

file_url=irsa_url+df2['uri'].iloc[0]
file_url
```

```{code-cell} ipython3
## Open the large FITS file without loading it entirely into memory
## pulling out just the extension we want for the 1D spectra of our object
response = requests.get(file_url)

with fits.open(BytesIO(response.content), memmap=True) as hdul:
    hdu = hdul[df2['hdu'].iloc[0]]
    dat = Table.read(hdu, format='fits', hdu=1)
    df_obj_irsa = dat.to_pandas()
```

```{code-cell} ipython3
## Now the data are read in, show an image

plt.plot(df_obj_irsa['WAVELENGTH'], df_obj_irsa['SIGNAL'])

plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux (erg / (Angstrom s cm2))')
# plt.ylim(10,50)
plt.title('Object ID is '+str(obj_id))
```

## Lets cut out a very small patch of the MER image to see what this galaxy looks like

```{code-cell} ipython3
## How large do you want the image cutout to be?
im_cutout= 2.0 * u.arcsec

## Use the ra and dec of the galaxy
ra = df_g_irsa[df_g_irsa['object_id']==obj_id]['ra'].iloc[0]
dec =  df_g_irsa[df_g_irsa['object_id']==obj_id]['dec'].iloc[0]

coords_cutout = SkyCoord(ra, dec, unit=(u.deg,u.deg), frame='icrs')
```

### Use fsspec to download a cutout of the fits file

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

# 5. Load the image on Firefly to be able to interact with the data directly

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

**Author**: Tiffany Meshkat (IPAC Scientist)

**Updated**: 2025-03-19

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
