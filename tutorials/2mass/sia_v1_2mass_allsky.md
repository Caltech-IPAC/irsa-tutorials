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
---

# Searching for 2MASS Images with SIA v1

+++

## Learning Goals

By the end of this tutorial, you will:

* Learn how to access IRSA's 2MASS All-Sky Atlas images. via the Simple Image Access (SIA) service.
* Use Python to identify which of IRSA's 2MASS All-Sky Atlas images cover a specified coordinate.
* Visualize one of the identified images using Forefly.
* Create and display a cutout of the downloaded image.

+++

## Introduction

The Two Micron All Sky Survey (2MASS) project uniformly scanned the entire sky in three near-infrared bands to detect and characterize point sources brighter than about 1 mJy in each band, with signal-to-noise ratio (SNR) greater than 10. More information about 2MASS can be found at:

https://irsa.ipac.caltech.edu/Missions/2mass.html

The [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu) at Caltech is the archive for 2MASS images and catalogs. The 2MASS images that are the subject of this tutorial are made accessible via the [International Virtual Observatory Alliance (IVOA)](https://ivoa.net) [Simple Image Access (SIA)](https://wiki.ivoa.net/internal/IVOA/SiaInterface/SIA-V2-Analysis.pdf) protocol. IRSA's 2MASS SIA service is registered in the NASA Astronomical Virtual Observatory (NAVO) [Directory](https://vao.stsci.edu). 

```{note}
IRSA supports both SIA v1 and SIA v2 protocols. The version used depends on the specific dataset. This IRSA [website](https://irsa.ipac.caltech.edu/ibe/sia.html) provides information on which version each service uses and how to access them. Further information on how to access IRSA data with different techniques is available [here](https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html). This tutorial uses SIA v1 for 2MASS allsky images.
This SIA v1 service is based on an older set of SIA protocols and is limited to the WISE/NEOWISE, 2MASS, and PTF datasets. It allows for only position-based searches to a single table. The IRSA SIA v1 search service has been superseded by the SIA v2 service for all other datasets.```

+++

## Imports
- `astropy.table` for reading SIA CSV responses
- `re` for robust band-name matching in metadata
- `astropy.coordinates` for defining coordinates
- `astropy.nddata` for creating an image cutout
- `astropy.wcs` for interpreting the World Coordinate System header keywords of a fits file
- `matplotlib.pyplot` for plotting
- `astropy.utils.data` for downloading files
- `astropy.io` to manipulate FITS files

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
%pip install matplotlib astropy jupyter_firefly_extensions
```

```{code-cell} ipython3
from astropy.table import Table
import re
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
import astropy.units as u
from firefly_client import FireflyClient
```

## 1. Define the target

+++

Define coordinates of a bright star

```{code-cell} ipython3
ra = 314.30417
dec = 77.595559
pos = SkyCoord(ra=ra, dec=dec, unit='deg')
```

## 2. Discover 2MASS All-Sky Atlas images

+++

IRSA provides Simple Image Access (SIA) services for various datasets. A list of available datasets and their access URLs can be found [here](https://irsa.ipac.caltech.edu/ibe/sia_v1.html).
This tutorial uses SIA v1 for 2MASS allsky images.

The URL for IRSA SIA v1 queries takes the form of:

https://irsa/ipac.caltech.edu/ibe/sia/{mission}/{data-set}/{table-name}?{query-constraints}

From the linked table we can see that:
- mission is "twomass"
- data-set is "allsky"
- table-name is "allsky"

To use the below code to examine images from other v1 datasets, substitute their values from the table into the url below.

```{code-cell} ipython3
# first we need to know the name of the dataset on the IRSA system
# Also choose a size of image to return (note this is in degrees)
base = "https://irsa.ipac.caltech.edu/ibe/sia/twomass/allsky/allsky"
url = f"{base}?POS={ra},{dec}&SIZE=0.01,0.01&RESPONSEFORMAT=CSV"
```

## 3. Search for images
Which images in the IRSA 2MASS dataset include our target of interest?

```{code-cell} ipython3
#get a table of all images within 1 arcsecond of our target position
im_table = Table.read(url, format="ascii.csv")
```

```{code-cell} ipython3
# Inspect the table that is returned
im_table
```

```{code-cell} ipython3
# Look at a list of the column names included in this table
im_table.colnames
```

```{code-cell} ipython3
# Let's look at the values in one of the columns
im_table['sia_url']
```

## 4.Locate and visualize an image of interest

We start by filtering the image results for the H band images.  
Then look at the header of one of the resulting H band images of our target star.
Finally, we create an interactive FITS display of the H image by using [Firefly](https://caltech-ipac.github.io/firefly_client/index.html), an open-source interactive visualization tool for astronomical data.
To understand how to open the Firefly viewer in a new tab from your Python notebook, refer to [this documentation](https://caltech-ipac.github.io/firefly_client/usage/initializing-vanilla.html) on how to initialize FireflyClient.

```{code-cell} ipython3
# You can put the URL from the column "access_url" into a browser to download the file. 
# Or you can work with it in Python, as shown below.

# Filter to H-band by parsing the column `sia_title`
#    Titles look like: "2MASS All-Sky Data Release H-Band Atlas Image: ..."
titles = im_table["sia_title"].astype(str)

# Robust match: look for " H-Band " or "H Band" variants (case-insensitive)
is_h = [bool(re.search(r"\bH[-\s]?Band\b", s, flags=re.IGNORECASE)) for s in titles]

im_table_h = im_table[is_h]
print(f"H-band rows: {len(im_table_h)}")
```

```{code-cell} ipython3
# Lets look at the access_url of the first one:
row = im_table_h[0]
image_url = row["sia_url"]
print("URL:",image_url)

#Use Astropy to examine the header of the URL from the previous step.
hdulist = fits.open(image_url)
hdulist.info()
```

```{code-cell} ipython3
# Uncomment when opening a Firefly viewer in a tab within Jupyter Lab with jupyter_firefly_extensions installed
fc = FireflyClient.make_lab_client()

# Uncomment when opening Firefly viewer in contexts other than the above 
#fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")

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

# make 1' x 1' cutout
cutout = Cutout2D(data, position=pos, size=1 * u.arcmin, wcs=wcs)

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

**Authors:** IRSA Data Science Team, including Troy Raen, Brigitta Sipőcz, Jessica Krick, Andreas Faisst, Jaladh Singhal, Vandana Desai, Dave Shupe

**Updated:** 19 February 2026

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.

**Runtime:** As of the date above, this notebook takes about 20 seconds to run to completion on a machine with 8GB RAM and 4 CPU.
This runtime is dependent on archive servers which means runtime will vary for users.

+++

## Citations

**Astropy:**
To see the Bibtex references for this, uncomment the below cell

**2MASS:**
This publication makes use of data products from the Two Micron All Sky Survey, which is a joint project of the University of Massachusetts and the Infrared Processing and Analysis Center/California Institute of Technology, funded by the National Aeronautics and Space Administration and the National Science Foundation."

```{code-cell} ipython3
#import astropy

#astropy.__citation__
```
