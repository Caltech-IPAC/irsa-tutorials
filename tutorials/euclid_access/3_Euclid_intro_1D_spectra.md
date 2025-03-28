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

# Euclid Q1: 1D spectra

+++

## Learning Goals

+++

By the end of this tutorial, you will:
- Understand the basic characteristics of Euclid Q1 SIR 1D spectra.
- What columns are available in the MER catalog.
- How to query with ADQL in the MER catalog.
- How to make a simple color-magnitude diagram with the data.

+++

## Introduction

Euclid launched in July 2023 as a European Space Agency (ESA) mission with involvement by NASA.
The primary science goals of Euclid are to better understand the composition and evolution of the dark Universe.
The Euclid mission is providing space-based imaging and spectroscopy as well as supporting ground-based imaging to achieve these primary goals.
These data will be archived by multiple global repositories, including IRSA, where they will support transformational work in many areas of astrophysics.

Euclid Quick Release 1 (Q1) consists of consists of ~30 TB of imaging, spectroscopy, and catalogs covering four non-contiguous fields:
Euclid Deep Field North (22.9 sq deg), Euclid Deep Field Fornax (12.1 sq deg), Euclid Deep Field South (28.1 sq deg), and LDN1641.

Among the data products included in the Q1 release are the 1D spectra created by the SIR Processing Function.
This notebook provides an introduction to these SIR 1D spectra.
If you have questions about it, please contact the [IRSA helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html).

+++

## Imports

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed
# !pip install matplotlib pandas requests astropy pyvo
```

```{code-cell} ipython3
from io import BytesIO

import matplotlib.pyplot as plt
import pandas as pd
import requests

from astropy.io import fits
from astropy.table import Table

import pyvo as vo
```

## 1. Download 1D spectra from IRSA directly to this notebook

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
table_1dspectra= 'euclid.objectid_spectrafile_association_q1'
table_phz= 'euclid_q1_phz_photo_z'
table_galaxy_candidates= 'euclid_q1_spectro_zcatalog_spe_galaxy_candidates'
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

## 2. Search for the spectrum of a specific galaxy in the 1D spectra table

```{code-cell} ipython3
obj_id=2739401293646823742

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

### Create the full filename/url

```{code-cell} ipython3
irsa_url='https://irsa.ipac.caltech.edu/'

file_url=irsa_url+df2['uri'].iloc[0]
file_url
```

## 3. Read in the spectrum using the file_url and the extension just for this object

Currently IRSA has the spectra stored in very large files containing multiple (14220) extensions with spectra of many targets within one tile. You can choose to read in the big file below to see what it looks like (takes a few mins to load) or skip this step and just read in the specific extension we want for the 1D spectra (recommended).

```{code-cell} ipython3
#### Code to read in the large file with many extensions and spectra from a tile
#### Currently commented out

# ## Complete file url with the irsa url at the start
# url = file_url
# response = requests.get(url)

# hdul = fits.open(BytesIO(response.content))  # Open FITS file from memory
# hdul.info()  # Show file info
```

### Open the large FITS file without loading it entirely into memory, pulling out just the extension we want for the 1D spectra of our object

```{code-cell} ipython3
response = requests.get(file_url)

with fits.open(BytesIO(response.content), memmap=True) as hdul:
    hdu = hdul[df2['hdu'].iloc[0]]
    dat = Table.read(hdu, format='fits', hdu=1)
    df_obj_irsa = dat.to_pandas()
```

### Plot the image of the extracted spectrum

- Convert the wavelength to microns

```{code-cell} ipython3
## Now the data are read in, show an image

## Converting from Angstrom to microns
plt.plot(df_obj_irsa['WAVELENGTH']/10000., df_obj_irsa['SIGNAL'])

plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux'+dat['SIGNAL'].unit.to_string('latex_inline'))
plt.title(obj_id)
```

## About this Notebook

**Author**: Tiffany Meshkat (IPAC Scientist)

**Updated**: 2025-03-19

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
