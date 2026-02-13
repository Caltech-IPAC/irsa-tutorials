---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Searching for contributed COSMOS images

+++

This notebook tutorial demonstrates the process of querying IRSA's Simple Image Access (SIA) service for the COSMOS images, making a cutout image (thumbnail), and displaying the cutout.


+++

## Learning Goals

By the end of this tutorial, you will:

* Learn how to access IRSA's COSMOS images via the Simple Image Access (SIA) service.
* Use the Python pyvo package to identify which of IRSA's COSMOS images cover a specified coordinate.
* Download one of the identified images.
* Create and display a cutout of the downloaded image.

+++ {"jp-MarkdownHeadingCollapsed": true}

## Introduction

The COSMOS Archive serves data taken for the Cosmic Evolution Survey with HST (COSMOS) project, using IRSA's general search service, Atlas. COSMOS is an HST Treasury Project to survey a 2 square degree equatorial field with the ACS camera. For more information about COSMOS, see:

https://irsa.ipac.caltech.edu/Missions/cosmos.html

The [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu) at Caltech is one of the archives for COSMOS images and catalogs. The COSMOS images that are the subject of this tutorial are made accessible via the [International Virtual Observatory Alliance (IVOA)](https://ivoa.net) [Simple Image Access (SIA)](https://wiki.ivoa.net/internal/IVOA/SiaInterface/SIA-V2-Analysis.pdf) protocol. IRSA's SEIP SIA service is registered in the NASA Astronomical Virtual Observatory (NAVO) [Directory](https://vao.stsci.edu). Based on the registered information, the Python package [pyvo](https://pyvo.readthedocs.io) can be used to query the SIA service for a list of images that meet specified criteria, and standard Python libraries can be used to download and manipulate the images.
Other datasets at IRSA are available through other SIA services:

https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html

```{note}
IRSA supports both SIA v1 and SIA v2 protocols. The version used depends on the specific dataset. The IRSA website at https://irsa.ipac.caltech.edu/ibe/sia.html provides information on which version each service uses and how to access them.
```


+++

## Imports

- `pyvo` for querying IRSA's COSMOS SIA service
- `astropy.coordinates` for defining coordinates
- `astropy.nddata` for creating an image cutout
- `astropy.wcs` for interpreting the World Coordinate System header keywords of a fits file
- `astropy.units` for attaching units to numbers passed to the SIA service
- `matplotlib.pyplot` for plotting
- `astropy.utils.data` for downloading files
- `astropy.io` to manipulate FITS files

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install matplotlib astropy pyvo
```

```{code-cell} ipython3
import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
```

## Section 1 - Setup

+++

Set images to display in the notebook

```{code-cell} ipython3
%matplotlib inline
```

Define coordinates of a bright source

```{code-cell} ipython3
ra = 149.99986
dec = 2.24875
pos = SkyCoord(ra=ra, dec=dec, unit='deg')
```

## Section 2 - Define a service for COSMOS images

+++

IRSA provides Simple Image Access (SIA) services for various datasets. A list of available datasets and their access URLs can be found at:

https://irsa.ipac.caltech.edu/ibe/sia.html

This tutorial uses SIA v1 for COSMOS images.

The COSMOS images service URL is:

https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas?mission=COSMOS&hdr_location=%5CCOSMOSDataPath%5C&collection_desc=Cosmic+Evolution+Survey+with+HST+%28COSMOS%29&SIAP_ACTIVE=1&

```{code-cell} ipython3
cosmos_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas?mission=COSMOS&hdr_location=%5CCOSMOSDataPath%5C&collection_desc=Cosmic+Evolution+Survey+with+HST+%28COSMOS%29&SIAP_ACTIVE=1&")
```

## Section 3 - Search the service

+++

Search for images covering within 1 arcsecond of the star

```{code-cell} ipython3
im_table = cosmos_service.search(pos=pos, size=1.0*u.arcsec)
```

Inspect the table of images that is returned

```{code-cell} ipython3
im_table
```

```{code-cell} ipython3
im_table.to_table().colnames
```

View the first ten entries of the table

```{code-cell} ipython3
im_table.to_table()[:10]
```

## Section 4 - Locate and download an image of interest

+++

Locate the first image in the band_name of i+

```{code-cell} ipython3
for i in range(len(im_table)):
    if im_table[i]['band_name'] == 'i+':
        break
print(im_table[i].getdataurl())
```

Download the image

```{code-cell} ipython3
fname = download_file(im_table[i].getdataurl(), cache=True)
image1 = fits.open(fname)
```

## Section 5 - Extract a cutout and plot it

```{code-cell} ipython3
wcs = WCS(image1[0].header)
```

Make a cutout centered on the position

```{code-cell} ipython3
cutout = Cutout2D(image1[0].data, pos, (60, 60), wcs=wcs)
wcs = cutout.wcs
```

```{code-cell} ipython3
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1, projection=wcs)
ax.imshow(cutout.data, cmap='gray_r', origin='lower')
ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=500, edgecolor='red', facecolor='none')
```

***

+++

## About this notebook

+++

**Author:** David Shupe, IRSA Scientist, and the IRSA Science Team

**Updated:** 2022-02-14

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.

+++

## Citations

+++

If you use `astropy` for published research, please cite the authors. Follow these links for more information about citing `astropy`:

* [Citing `astropy`](https://www.astropy.org/acknowledging.html)

+++

If you use COSMOS ACS imaging data in published research,  please cite the dataset Digital Object Identifier (DOI): [10.26131/IRSA178](https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA178).
