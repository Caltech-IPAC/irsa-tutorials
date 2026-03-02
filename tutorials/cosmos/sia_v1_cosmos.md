---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  name: python3
  display_name: python3
  language: python
authors:
  - name: IRSA Data Science Team
  - name: Troy Raen
  - name: Brigitta Sipőcz
  - name: Jessica Krick
  - name: Andreas Faisst
  - name: Jaladh Singhal
  - name: Vandana Desai
  - name: Dave Shupe
---

# Searching for contributed COSMOS images

+++

This notebook tutorial demonstrates the process of querying IRSA's Simple Image Access (SIA) service for the COSMOS images, making a cutout image (thumbnail), and displaying the cutout.

+++

## Learning Goals

By the end of this tutorial, you will:

* Learn how to search the NASA Astronomical Virtual Observatory Directory web portal for a service that provides access to IRSA's COSMOS images.
* Use the Python pyvo package to identify which of IRSA's COSMOS images cover a specified coordinate.
* Download one of the identified images.
* Create and display a cutout of the downloaded image.

+++ {"jp-MarkdownHeadingCollapsed": true}

## Introduction

The COSMOS Archive serves data taken for the Cosmic Evolution Survey with HST (COSMOS) project, using IRSA's general search service, Atlas. COSMOS is an HST Treasury Project to survey a 2 square degree equatorial field with the ACS camera. For more information about COSMOS, see:

https://irsa.ipac.caltech.edu/Missions/cosmos.html

The [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu) at Caltech is one of the [archives](https://irsa.ipac.caltech.edu/Missions/cosmos.html) for COSMOS images and catalogs. The COSMOS images that are the subject of this tutorial are made accessible via the [International Virtual Observatory Alliance (IVOA)](https://ivoa.net) [Simple Image Access (SIA)](https://wiki.ivoa.net/internal/IVOA/SiaInterface/SIA-V2-Analysis.pdf) protocol.

```{note}
IRSA supports both SIA v1 and SIA v2 protocols. The version used depends on the specific dataset. This IRSA [website](https://irsa.ipac.caltech.edu/ibe/sia.html) provides information on which version each service uses and how to access them. Further information on how to access IRSA data with different techniques is available [here](https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html). This tutorial uses SIA v1 for COSMOS images.
This SIA v1 service is based on an older set of SIA protocols and is limited to the COSMOS, WISE, 2MASS, and PTF datasets. It allows for only position-based searches to a single table. The IRSA SIA v1 search service has been superseded by the SIA v2 service for datasets other than COSMOS and PTF.
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
!pip -q install matplotlib astropy pyvo jupyter_firefly_extensions
```

```{code-cell} ipython3
import pyvo as vo
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from firefly_client import FireflyClient
```

## 1. Define the target
Define coordinates of a bright star

```{code-cell} ipython3
ra = 149.99986
dec = 2.24875
pos = SkyCoord(ra=ra, dec=dec, unit='deg')
```

## 2. Discover COSMOS images

```{code-cell} ipython3
cosmos_service = vo.dal.SIAService("https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas?mission=COSMOS&hdr_location=%5CCOSMOSDataPath%5C&collection_desc=Cosmic+Evolution+Survey+with+HST+%28COSMOS%29&SIAP_ACTIVE=1&")
```

## 3. Search for images
Which images in the COSMOS dataset include our target of interest?

```{code-cell} ipython3
# Get a table of all images that cover this position
# This service actually returns a cutout of whatever size you choose
im_table = cosmos_service.search(pos=pos, size=150*u.arcsec)
```

```{code-cell} ipython3
# Inspect the top of the table that is returned
im_table.to_table()[:10]
```

```{code-cell} ipython3
# Look at a list of the column names included in this table
im_table.to_table().colnames
```

```{code-cell} ipython3
# Let's look at the unique values in one of the columns
print(np.unique(im_table['band_name']))
```

##

+++

## 4.Locate and visualize an image of interest

We start by filtering the image results for the first IRAC1 band images.
Then look at the header of one of the resulting image of our target star.
Finally, we create an interactive FITS display of the IRAC1 image by using [Firefly](https://caltech-ipac.github.io/firefly_client/index.html), an open-source interactive visualization tool for astronomical data.
To understand how to open the Firefly viewer in a new tab from your Python notebook, refer to [this documentation](https://caltech-ipac.github.io/firefly_client/usage/initializing-vanilla.html) on how to initialize FireflyClient.

```{code-cell} ipython3
# You can put the URL from the column "sia_url" into a browser to download the file.
# Or you can work with it in Python, as shown below.

im_table_astropy = im_table.to_table()

irac1_rows = im_table_astropy[
    im_table_astropy['band_name'] == 'IRAC1'
]
```

```{code-cell} ipython3
# Lets look at the data access url in the column named 'sia_url'.
# We will focus on the first image for now.
image_url = irac1_rows[0]['sia_url']
print(image_url)
```

```{code-cell} ipython3
#Use Astropy to examine the header of the URL from the previous step.
hdulist = fits.open(image_url)
hdulist.info()
```

```{code-cell} ipython3
# Uncomment when opening a Firefly viewer in a tab within Jupyter Lab with jupyter_firefly_extensions installed
#fc = FireflyClient.make_lab_client()

# Uncomment when opening Firefly viewer in contexts other than the above
fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")

# Visualize an image by sending its URL to the viewer.
fc.show_fits_image(file_input=image_url,
             plot_id="image",
             Title="Image"
             )

#Try use the interactive tools in the viewer to explore the data.
```

## 5. Extract a cutout and plot it
If you want to see just a cutout of a certain region around the target, we do that below using astropy's Cutout2D.

```{code-cell} ipython3
data = hdulist[0].data
wcs = WCS(hdulist[0].header)

# make 0.5' x 0.5' cutout
cutout = Cutout2D(data, position=pos, size=0.5 * u.arcmin, wcs=wcs)

# display
plt.figure()
plt.imshow(cutout.data, origin='lower')
plt.colorbar()
```

***

+++

## Acknowledgements

- [Caltech/IPAC-IRSA](https://irsa.ipac.caltech.edu/)

+++


## About this notebook

**Updated:** 2 March 2026

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.

**Runtime:** As of the date above, this notebook takes about 20 seconds to run to completion on a machine with 8GB RAM and 4 CPU.
This runtime is dependent on archive servers which means runtime will vary for users.

+++

## Citations

**Astropy:**
To see the Bibtex references for this, uncomment the below cell

**COSMOS:**
If you use COSMOS ACS imaging data in published research,  please cite the dataset Digital Object Identifier (DOI): [10.26131/IRSA178](https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA178).

```{code-cell} ipython3
#import astropy

#astropy.__citation__
```
