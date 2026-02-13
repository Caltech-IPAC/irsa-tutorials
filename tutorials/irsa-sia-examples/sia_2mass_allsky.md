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

# Searching for 2MASS All-Sky Atlas Images

This notebook tutorial demonstrates the process of querying IRSA's Simple Image Access (SIA) service for the 2MASS All-Sky Atlas, making a cutout image (thumbnail), and displaying the cutout.

+++

***

+++ {"jp-MarkdownHeadingCollapsed": true}

## Learning Goals

By the end of this tutorial, you will:

* Learn how to access IRSA's 2MASS images via the Simple Image Access (SIA) service.
* Use the Python pyvo package to identify which of IRSA's 2MASS images cover a specified coordinate.
* Download one of the identified images.
* Create and display a cutout of the downloaded image.

+++

## Introduction

The Two Micron All Sky Survey (2MASS) project uniformly scanned the entire sky in three near-infrared bands to detect and characterize point sources brighter than about 1 mJy in each band, with signal-to-noise ratio (SNR) greater than 10. More information about 2MASS can be found at:

https://irsa.ipac.caltech.edu/Missions/2mass.html

The [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu) at Caltech is the archive for 2MASS images and catalogs. The 2MASS images that are the subject of this tutorial are made accessible via the [International Virtual Observatory Alliance (IVOA)](https://ivoa.net) [Simple Image Access (SIA)](https://wiki.ivoa.net/internal/IVOA/SiaInterface/SIA-V2-Analysis.pdf) protocol. IRSA's 2MASS SIA service is registered in the NASA Astronomical Virtual Observatory (NAVO) [Directory](https://vao.stsci.edu). Based on the registered information, the Python package [pyvo](https://pyvo.readthedocs.io) can be used to query the 2MASS SIA service for a list of images that meet specified criteria, and standard Python libraries can be used to download and manipulate the images.
Other datasets at IRSA are available through other SIA services:

https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html

```{note}
IRSA supports both SIA v1 and SIA v2 protocols. The version used depends on the specific dataset. The IRSA website at https://irsa.ipac.caltech.edu/ibe/sia.html provides information on which version each service uses and how to access them.
```

+++

## Imports

- `pyvo` for querying IRSA's 2MASS SIA service
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

Define coordinates of a bright star

```{code-cell} ipython3
ra = 314.30417
dec = 77.595559
pos = SkyCoord(ra=ra, dec=dec, unit='deg')
```

## Section 2 - Define a service for 2MASS images

+++

IRSA provides Simple Image Access (SIA) services for various datasets. A list of available datasets and their access URLs can be found at:

https://irsa.ipac.caltech.edu/ibe/sia.html

This tutorial uses SIA v1 for 2MASS All-Sky Atlas images.

The 2MASS All-Sky Atlas images service URL is:

https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&

```{code-cell} ipython3
twomass_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?type=at&ds=asky&")
```

## Section 3 - Search the service

+++

Search for images covering within 1 arcsecond of the star

```{code-cell} ipython3
im_table = twomass_service.search(pos=pos, size=1.0*u.arcsec)
```

Examine the table of images that is returned

```{code-cell} ipython3
im_table.to_table()
```

## Section 4 - Locate and download an image of interest

+++

Locate the first H-band image and display its URL

```{code-cell} ipython3
for i in range(len(im_table)):
    if im_table[i]['band'] == 'H':
        break
print(im_table[i].getdataurl())
```

Download the image and open it in Astropy

```{code-cell} ipython3
fname = download_file(im_table[i].getdataurl(), cache=True)
image1 = fits.open(fname)
```

## Section 5 - Extract a cutout and plot it

```{code-cell} ipython3
wcs = WCS(image1[0].header)
```

```{code-cell} ipython3
cutout = Cutout2D(image1[0].data, pos, (60, 60), wcs=wcs)
wcs = cutout.wcs
```

```{code-cell} ipython3
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1, projection=wcs)
ax.imshow(cutout.data, cmap='gray_r', origin='lower',
          vmax = 1000)
ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=500, edgecolor='red', facecolor='none')
```

## Exercise

+++

Repeat the steps above to retrieve a cutout from the AllWISE Atlas images

```{code-cell} ipython3

```

***

+++

## About this notebook

+++

**Author:** David Shupe, IRSA Scientist, and the IRSA Science Team

**Updated:** 2023-02-16

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.


+++

## Citations

+++

If you use `astropy` for published research, please cite the authors. Follow these links for more information about citing `astropy`:

* [Citing `astropy`](https://www.astropy.org/acknowledging.html)

+++

If you use 2MASS data in published research, please cite the canonical paper [Skrutskie et al (2006)](http://adsabs.harvard.edu/abs/2006AJ....131.1163S), and include the following standard acknowledgment:

*"This publication makes use of data products from the Two Micron All Sky Survey, which is a joint project of the University of Massachusetts and the Infrared Processing and Analysis Center/California Institute of Technology, funded by the National Aeronautics and Space Administration and the National Science Foundation."*

Please also cite the doi for the 2MASS All-Sky Atlas Image Service at https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA121
