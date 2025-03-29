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
# !pip install matplotlib astropy 'astroquery>=0.4.10'
```

```{code-cell} ipython3
import urllib
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import QTable
from astropy.visualization import quantity_support

from astroquery.ipac.irsa import Irsa
```

## 1. Search for the spectra of a specific galaxy

First, explore what Euclid catalogs are available. Note that we need to use the object ID for our targets to be able to download their spectrum.

Search for all tables in IRSA labeled as "euclid".

```{code-cell} ipython3
Irsa.list_catalogs(filter='euclid')
```

```{code-cell} ipython3
table_1dspectra = 'euclid.objectid_spectrafile_association_q1'
```

## 2. Search for the spectrum of a specific galaxy in the 1D spectra table

```{code-cell} ipython3
obj_id = 2739401293646823742
```

We will use TAP and an ASQL query to find the spectral data for our galaxy. (ADQL is the [IVOA Astronomical Data Query Language](https://www.ivoa.net/documents/latest/ADQL.html) and is based on SQL.)

```{code-cell} ipython3
adql_object = f"SELECT * FROM {table_1dspectra} WHERE objectid = {obj_id} AND uri IS NOT NULL "

# Pull the data on this particular galaxy
result = Irsa.query_tap(adql_object).to_table()
```

Pull out the file name from the ``result`` table:

```{code-cell} ipython3
file_uri = urllib.parse.urljoin(Irsa.tap_url, result['uri'][0])
file_uri
```

## 3. Read in the spectrum for only our specific object

Currently IRSA has the spectra stored in very large files containing multiple (14220) extensions with spectra of many targets within one tile. You can choose to read in the big file below to see what it looks like (takes a few mins to load) or skip this step and just read in the specific extension we want for the 1D spectra (recommended).

```{code-cell} ipython3
# hdul = fits.open(file_uri)
# hdul.info()
```

Open the large FITS file without loading it entirely into memory, pulling out just the extension we want for the 1D spectra of our object

```{code-cell} ipython3
with fits.open(file_uri) as hdul:
    spectra = QTable.read(hdul[result['hdu'][0]], format='fits')
```

```{code-cell} ipython3
spectra
```

## 4. Plot the image of the extracted spectrum

```{tip}
As we use astropy.visualization's ``quantity_support``, matplotlib automatically picks up the axis units from the quantitites we plot.
```

```{code-cell} ipython3
quantity_support()
```

```{code-cell} ipython3
plt.plot(spectra['WAVELENGTH'], spectra['SIGNAL'])
plt.title(obj_id)
```

## About this Notebook

**Author**: Tiffany Meshkat (IPAC Scientist)

**Updated**: 2025-03-19

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
