---
authors:
- name: IRSA Data Science Team
- name: Troy Raen
- name: "Brigitta Sip\u0151cz"
- name: Jessica Krick
- name: Andreas Faisst
- name: Jaladh Singhal
- name: Vandana Desai
- name: Dave Shupe
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


## Introduction

The COSMOS Archive at IRSA serves data taken for the Cosmic Evolution Survey with HST (COSMOS) project. COSMOS is an HST Treasury Project to survey a 2 square degree equatorial field with the ACS camera. For more information about COSMOS, see:

https://irsa.ipac.caltech.edu/Missions/cosmos.html

The [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu) at Caltech is one of the [archives](https://irsa.ipac.caltech.edu/Missions/cosmos.html) for COSMOS images and catalogs. The COSMOS images that are the subject of this tutorial are made accessible via the [International Virtual Observatory Alliance (IVOA)](https://ivoa.net) [Simple Image Access (SIA)](https://wiki.ivoa.net/internal/IVOA/SiaInterface/SIA-V2-Analysis.pdf) protocol.

```{note}
IRSA supports both SIA v1 and SIA v2 protocols. The version used depends on the specific dataset. This IRSA [website](https://irsa.ipac.caltech.edu/ibe/sia.html) provides information on which version each service uses and how to access them. Further information on how to access IRSA data with different techniques is available [here](https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html). This tutorial uses SIA v1 for COSMOS images.
This SIA v1 service is based on an older set of SIA protocols and is limited to the COSMOS, WISE, 2MASS, and PTF datasets. It allows for only position-based searches to a single table. The IRSA SIA v1 search service has been superseded by the SIA v2 service for datasets other than COSMOS and PTF.
```

+++

## Imports

- `pyvo` for discovering and querying IRSA's COSMOS SIA service
- `numpy` for working with tables
- `astropy.coordinates` for defining coordinates
- `astropy.nddata` for creating an image cutout
- `astropy.wcs` for interpreting the World Coordinate System header keywords of a fits file
- `astropy.units` for attaching units to numbers passed to the SIA service
- `matplotlib.pyplot` for plotting
- `astropy.io` to manipulate FITS files
- `firefly_client` for visualizing data

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
!pip -q install matplotlib astropy pyvo jupyter_firefly_extensions
```

```{code-cell} ipython3
from pyvo import regsearch
import numpy as np
import astropy
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from firefly_client import FireflyClient
from astropy.visualization import simple_norm
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
# Search the Virtual Observatory Registry for image services at IRSA associated with the COSMOS survey.
image_services = regsearch(
    servicetype='sia1',
    keywords=['cosmos', 'irsa']
)

# Make sure we got what we were looking for
for i, r in enumerate(image_services):
    print(f"{i:2d}  {r.short_name:20s}  {r.res_title}")

# Turn the result into a usable image access service
resource = image_services[0] 
cosmos_service = resource.get_service("sia1")
```

## 3. Search for images
Which images in the COSMOS dataset include our target of interest?

```{code-cell} ipython3
# Get a table of all images that cover this position.
# Choose the size of the returned image. 
im_results = cosmos_service.search(pos=pos, size=150*u.arcsec)

# Convert the PyVO result to an Astropy Table
im_table = im_results.to_table()
```

```{code-cell} ipython3
# Inspect the top of the table that is returned
im_table[:10]
```

```{code-cell} ipython3
# Look at a list of the column names included in this table
im_table.colnames
```

```{code-cell} ipython3
# Let's look at the unique values in one of the columns
print(np.unique(im_results['band_name']))
```

##

+++

## 4. Locate and visualize an image of interest

We start by filtering the image results for the first IRAC1 band images.
Then look at the header of one of the resulting image of our target star.
Finally, we create an interactive FITS display of the IRAC1 image by using [Firefly](https://caltech-ipac.github.io/firefly_client/index.html), an open-source interactive visualization tool for astronomical data.
To understand how to open the Firefly viewer in a new tab from your Python notebook, refer to [this documentation](https://caltech-ipac.github.io/firefly_client/usage/initializing-vanilla.html) on how to initialize FireflyClient.

```{code-cell} ipython3
# You can put the URL from the column "sia_url" into a browser to download the file.
# Or you can work with it in Python, as shown below.

irac1_rows = im_table[im_table['band_name'] == 'IRAC1']
```

```{code-cell} ipython3
# Let's look at the data access url in the column named 'sia_url'.
# We will focus on the first image for now.
image_url = irac1_rows[0]['sia_url']
print(image_url)
```

```{code-cell} ipython3
# Use Astropy to examine the header of the URL from the previous step,
# and grab the data and wcs from the header.
with fits.open(image_url, memmap=False) as hdul:
    hdul.info()           
    data = hdul[0].data
    wcs = WCS(hdul[0].header)
    
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

#Try using the interactive tools in the viewer to explore the data.
```

## 5. Extract a cutout and plot it
If you want to see just a cutout of a certain region around the target, we do that below using astropy's Cutout2D.

```{code-cell} ipython3
# make 0.5' x 0.5' cutout
cutout = Cutout2D(data, position=pos, size=0.5 * u.arcmin, wcs=wcs)

#add quick normalization/stretch
norm = simple_norm(cutout.data, stretch="sqrt", percent=99)

# display
plt.figure()
plt.imshow(cutout.data, origin='lower', norm = norm)
plt.colorbar(label="Image value")
plt.title("COSMOS IRAC1 (quicklook)")
plt.xlabel("Pixel X")
plt.ylabel("Pixel Y")
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

**COSMOS:**
If you use COSMOS ACS imaging data in published research,  please cite the dataset Digital Object Identifier (DOI): [10.26131/IRSA178](https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA178).

**Astropy:**
This work made use of [Astropy](http://www.astropy.org) a community-developed core Python package and an ecosystem of tools and resources for astronomy (Astropy Collaboration et al., 2013, Astropy Collaboration et al., 2018, Astropy Collaboration et al.,2022}.
