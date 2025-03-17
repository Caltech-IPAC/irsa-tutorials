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

# Euclid Quick Release 1: SPE catalog

+++

## Learning Goals

+++

By the end of this tutorial, you will: 
- Understand the basic characteristics of Euclid Q1 SPE catalogs.
- Understand what SPE catalogs are available and how to view the columns in those catalogs.
- How to query with ADQL in the SPE lines catalog to find strong H-alpha detections.
- How to make a plot the detected line features over the 1D spectra.

+++

## Introduction

+++

Euclid is a European Space Agency (ESA) space mission with NASA participation, to study the geometry and nature of the dark Universe. The Quick Data Release 1 (Q1) are the first data release from the Euclid mission after the Early Release Observations (ERO). On March 19, 2025 the data will be available on the ESA archive (https://easidr.esac.esa.int/sas/) and on the IRSA archive (https://irsa.ipac.caltech.edu).

These notebooks focus on how to access, download, and process Euclid Q1 data from the IRSA archive. At the end of the notebook, we also include some information for how to access the Q1 data from the ESA archive. If you have any issues accessing data from the archives, please contact the helpdesk directly: IRSA (irsasupport@ipac.caltech.edu) and ESA (https://support.cosmos.esa.int/euclid).

Every one dimensional spectrum is processed through a template and line fitting pipeline, producing several different 'SPE' catalogs. This notebook provides an introduction to the SPE catalogs released as part of Euclid Q1. Other Euclid notebooks show how to use other data products released as part of Euclid Q1.

+++

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# !pip install matplotlib pandas astropy pyvo
```

```{code-cell} ipython3
from io import BytesIO
import re

import matplotlib.pyplot as plt
import pandas as pd
import requests

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch

import pyvo as vo
```

## 1. Find the MER Tile ID that corresponds to a given RA and Dec

In this case, choose the coordinates from the first notebook to save time downloading the MER mosaic. Search a radius of 1.5 arcminutes around these coordinates.

```{code-cell} ipython3
search_radius = 10 * u.arcsec
coord = SkyCoord.from_name('HD 168151')
```

### Use IRSA to search for all Euclid data on this target

This searches specifically in the euclid_DpdMerBksMosaic "collection" which is the MER images and catalogs.

```{code-cell} ipython3
irsa_service= vo.dal.sia2.SIA2Service('https://irsa.ipac.caltech.edu/SIA')

im_table = irsa_service.search(pos=(coord, search_radius), collection='euclid_DpdMerBksMosaic')

## Convert the table to pandas dataframe
df_im_irsa=im_table.to_table().to_pandas()
```

```{code-cell} ipython3
## Change the settings so we can see all the columns in the dataframe and the full column width 
## (to see the full long URL)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
```

#### This dataframe contains other non-Euclid datasets that have been "Euclidized", meaning they have been put on the same pixel scale as the Euclid data. For this example we just want to look at the Euclid data, so select Euclid for the facility name, and choose science as the data product subtype.

```{code-cell} ipython3
df_im_euclid=df_im_irsa[ (df_im_irsa['dataproduct_subtype']=='science') &  (df_im_irsa['facility_name']=='Euclid')]

df_im_euclid.head()
```

## Choose the VIS image and pull the filename:

```{code-cell} ipython3
filename=df_im_euclid[df_im_euclid['energy_bandpassname']=='VIS']['access_url'].to_list()[0]

# ## Extract the tileID from the filename
tileID=re.search(r'TILE\s*(\d{9})', filename).group(1)

print('The MER tile ID for this object is :',tileID)
```

## 2. Read in the MER image from IRSA directly

```{code-cell} ipython3
######### TEMP
######## Note to testers, for now we need to replace the irsa.ipac.caltech.edu url with irsa
######## This will not be the same after the data are made public so this cell will be deleted at that time
def add_dev_to_domain(domain):
    parts = domain.split('.', 1)  # Split at the first dot
    if len(parts) == 2:
        return f"{parts[0]}dev.{parts[1]}"
    return domain

filename_dev = add_dev_to_domain(filename)
print(filename_dev)  

#####################
```

### Download the MER image to this notebook
Note this file is about 1.46 GB

```{code-cell} ipython3
fname = download_file(filename_dev, cache=True)
hdu_mer_irsa = fits.open(fname)
head_mer_irsa = hdu_mer_irsa[0].header

print(hdu_mer_irsa.info())
```

#### Extract just the primary image

```{code-cell} ipython3
im_mer_irsa=hdu_mer_irsa[0].data
```

#### Make a quick and simple plot to show the full MER image, with its large FOV

```{code-cell} ipython3
plt.imshow(im_mer_irsa, cmap='gray', origin='lower', 
           norm=ImageNormalize(im_mer_irsa, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
colorbar = plt.colorbar()
```

## 3. Download SPE catalog from IRSA directly to this notebook

Search for all tables in IRSA labeled as euclid

```{code-cell} ipython3
service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

tables = service.tables
for tablename in tables.keys():
    if "tap_schema" not in tablename and "euclid" in tablename:
            tables[tablename].describe()
            
```

```{code-cell} ipython3
table_mer= 'euclid_q1_mer_catalogue'
table_galaxy_candidates= 'euclid_q1_spectro_zcatalog_spe_galaxy_candidates'
table_1dspectra= 'euclid.objectid_spectrafile_association_q1'
table_spe= 'euclid_q1_spe_lines_line_features'
```

### Learn some information about the table:
- How many columns are there?
- List the column names

```{code-cell} ipython3
columns = tables[table_galaxy_candidates].columns
print(len(columns))
```

```{code-cell} ipython3
for col in columns:
    print(f'{f"{col.name}":30s}  {col.unit}  {col.description}') ## Currently no descriptions
```

```{code-cell} ipython3
## Change the settings so we can see all the columns in the dataframe and the full column width 
## (to see the full long URL)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)


## Can use the following lines to reset the max columns and column width of pandas
# pd.reset_option('display.max_columns')
# pd.reset_option('display.max_colwidth')
```

## Find some objects with spectra in our tileID

We specify the following conditions on our search: 
- The two signal to noise ratio columns (spe_line_snr_gf and spe_line_snr_di) should be greater than 5
- We want to detect H-alpha.
- We choose in which tileID to search, usign the tileID from the first notebook.
- Choose spectroscopic redshift (spe_z) beween 1.4 and 1.6 and spe_z_prob greater than 0.999

Finally we sort the data by descending spe_line_snr_gf to have the largest SNR H-alpha lines detected at the top.

```{code-cell} ipython3
adql = f"SELECT DISTINCT mer.object_id,mer.ra, mer.dec, mer.tileid, mer.flux_y_templfit, \
spe.spe_line_snr_gf,spe.spe_line_snr_di, spe.spe_line_name, spe.spe_line_central_wl_gf,\
spe.spe_line_ew_gf, galaxy.spe_z_err, galaxy.spe_z,galaxy.spe_z_prob \
FROM {table_mer} AS mer \
JOIN {table_spe} AS spe \
ON mer.object_id = spe.object_id \
JOIN {table_galaxy_candidates} AS galaxy \
ON mer.object_id = galaxy.object_id \
WHERE spe.spe_line_snr_gf >5 \
AND spe.spe_line_snr_di > 5 \
AND spe.spe_line_name = 'Halpha' \
AND mer.tileid = {tileID} \
AND galaxy.spe_z_prob > 0.999 \
AND galaxy.spe_z BETWEEN 1.4 AND 1.6 \
ORDER BY spe.spe_line_snr_gf DESC \
"

# Use TAP with this ADQL string using pyvo
result = service.search(adql)

# Convert table to pandas dataframe and drop duplicates
result_table = result.to_qtable()
```

### Choose an object of interest, lets look at an object with a strong Halpha line detected with high SNR.

```{code-cell} ipython3
obj_id = 2739401293646823742

obj_2739401293646823742 = result_table[(result_table['object_id'] == obj_id)]

obj_2739401293646823742
```

### Pull the spectrum of this object

```{code-cell} ipython3
adql_object = f"SELECT *  FROM {table_1dspectra}  WHERE objectid = {obj_id} AND uri IS NOT NULL "

result2 = service.search(adql_object)
df2 = result2.to_table().to_pandas()
df2
```

### The following steps to read in the spectrum follows the 3_Euclid_intro_1D_spectra notebook.

This involves reading in the spectrum without readin in the full FITS file, just pulling the extension we want.

```{code-cell} ipython3
irsa_url = 'https://irsa.ipac.caltech.edu/'

file_url = irsa_url + df2['uri'].iloc[0]
file_url

response = requests.get(file_url)

with fits.open(BytesIO(response.content), memmap=True) as hdul:
    hdu = hdul[df2['hdu'].iloc[0]]
    dat = Table.read(hdu, format='fits', hdu=1)
    df_obj_irsa = dat.to_pandas()
    
```

### Now the data are read in, plot the spectrum with the H-alpha line labeled

Divide by 10000 to convert from Angstrom to micron

```{code-cell} ipython3
wavelengths = df_obj['spe_line_central_wl_gf']/10000.
line_names = df_obj['spe_line_name']
snr_gf=df_obj['spe_line_snr_gf']

plt.plot(df_obj_irsa['WAVELENGTH']/10000., df_obj_irsa['SIGNAL'])

for wl, name,snr in zip(wavelengths, line_names,snr_gf):
    plt.axvline(wl, color='b', linestyle='--', alpha=0.3)
    plt.text(wl+0.02, .1, name+' SNR='+str(round(snr)), rotation=90, ha='center', va='bottom', fontsize=10)


plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux (erg / (Angstrom s cm2))')
plt.title(obj_id)
```

## About this Notebook

**Author**: Tiffany Meshkat (IPAC Scientist)

**Updated**: March 19, 2025

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
