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

# Introduction to Euclid Q1 1D spectra

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

+++

Euclid is a European Space Agency (ESA) space mission with NASA participation, to study the geometry and nature of the dark Universe. The Quick Data Release 1 (Q1) are the first data release from the Euclid mission after the Early Release Observations (ERO). On March 19, 2025 the data will be available on the ESA archive (https://easidr.esac.esa.int/sas/) and on the IRSA archive (https://irsa.ipac.caltech.edu).

These notebooks focus on how to access, download, and process Euclid Q1 data from the IRSA archive. At the end of the notebook, we also include some information for how to access the Q1 data from the ESA archive. If you have any issues accessing data from the archives, please contact the helpdesk directly: IRSA (irsasupport@ipac.caltech.edu) and ESA (https://support.cosmos.esa.int/euclid).

For the Euclid Wide Survey standard operating mode, the telescope undertakes a 4-point dither pattern. At each position VIS and NISP each take a 570s exposure, consisting of a direct visible image and a red grism exposure. This is followed by further NISP exposures in the Y, J, and H band filters (112 seconds each). The telescope is then dithered, and the sequence is repeated starting with a different grism position angle. There are actually two operational grisms oriented 180 degrees from each other. Each grism which will be used twice in this sequence, but with slight angular offsets (+/- 4 degrees), effectively creating the four different grism angles (Scaramella et al. 2022, A&A 662, A112).	

Data which can be obtained for SIR include: SIR "images", which effectively show the full image of objects with the spectral traces overlapping, and SIR 1D spectra for individual objects. Below we will describe how to access and process the 1D spectra products. For most users, simply accessing th 1D spectra is probably the preferred option, unless they would like to extract the spectrum again, or inspect the images to see if there is any artifact which might add noise to the spectrum.

This notebook provides an introduction to the SIR 1D spectra released as part of Euclid Q1. Other Euclid notebooks show how to use other data products released as part of Euclid Q1.

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
service = vo.dal.TAPService("https://irsadev.ipac.caltech.edu/TAP")

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
irsa_url='https://irsadev.ipac.caltech.edu/'

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

## Exercise

+++

### Optional -- Access the data from the ESA archive website directly

+++

#### 1. Download the data from the ESA archive

- Go to https://easidr.esac.esa.int/sas/ and sign in with your credentials.
- Go to the ADQL form and do the following search:

```{raw-cell}
SELECT TOP 10 spec.source_id 
FROM sedm.spectra_source as spec
JOIN catalogue.phz_classification AS phz_class
ON phz_class.object_id=spec.source_id
WHERE phz_class.phz_classification = 2
```

- This shows the first 10 sources with spectra in the list.
- Click the link/chain
- Make sure "source_id" is selected under "IDs columns"
- Click "show data" then click the download button next to "Spectra -- source_id"

NOTE: You need to unzip the file to get a fits file with the combined spectra. On a mac you can do this by just double clicking the file to unzip.

+++

## About this Notebook

**Author**: Tiffany Meshkat (IPAC Scientist)

**Updated**: March 19, 2025

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
