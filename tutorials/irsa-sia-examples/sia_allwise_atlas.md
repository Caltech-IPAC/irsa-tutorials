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

# Searching for AllWISE Atlas Images

+++

## Learning Goals

By the end of this tutorial, you will:

* Learn how to access IRSA's WISE AllWISE Atlas (L3a) coadded images via the Simple Image Access (SIA) service.
* identify which of IRSA's AllWISE Atlas images cover a specified coordinate.
* Visualize one of the identified images using Forefly.
* Create and display a cutout of the downloaded image.

+++

## Introduction

The AllWISE program builds upon the work of the successful Wide-field Infrared Survey Explorer mission [(WISE; Wright et al. 2010)](http://adsabs.harvard.edu/abs/2010AJ....140.1868W) by combining data from the WISE cryogenic and NEOWISE [(Mainzer et al. 2011 ApJ, 731, 53)](http://adsabs.harvard.edu/abs/2011ApJ...731...53M) post-cryogenic survey phases to form the a comprehensive view of the full mid-infrared sky. The AllWISE Images Atlas is comprised of 18,240 4-band calibrated 1.56°x1.56° FITS images, depth-of-coverage and noise maps, and image metadata produced by coadding nearly 7.9 million Single-exposure images from all survey phases. For more information about the WISE mission, see:

https://irsa.ipac.caltech.edu/Missions/wise.html

The [NASA/IPAC Infrared Science Archive (IRSA)](https://irsa.ipac.caltech.edu) at Caltech is the archive for AllWISE images and catalogs. The AllWISE Atlas images that are the subject of this tutorial are made accessible via the [International Virtual Observatory Alliance (IVOA)](https://ivoa.net) [Simple Image Access (SIA)](https://wiki.ivoa.net/internal/IVOA/SiaInterface/SIA-V2-Analysis.pdf) protocol. 


```{note}
IRSA supports both SIA v1 and SIA v2 protocols. The version used depends on the specific dataset. This IRSA [website](https://irsa.ipac.caltech.edu/ibe/sia.html) provides information on which version each service uses and how to access them. Further information on how to access IRSA data with different techniques is available [here](https://irsa.ipac.caltech.edu/docs/program_interface/api_images.html)
```

+++

## Imports
- `astropy.coordinates` for defining coordinates
- `astropy.nddata` for creating an image cutout
- `astropy.wcs` for interpreting the World Coordinate System header keywords of a fits file
- `astropy.units` for attaching units to numbers passed to the SIA service
- `matplotlib.pyplot` for plotting
- `astropy.io` to manipulate FITS files
- `firefly_client` for visuzlizing images
- `astroquery/ipac.irsa` for data access

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
%pip install matplotlib astropy jupyter_firefly_extensions
```

```{code-cell} ipython3
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from firefly_client import FireflyClient
from astroquery.ipac.irsa import Irsa
```

## 1. Define the target

+++

Define coordinates of a bright star

```{code-cell} ipython3
ra = 314.30417
dec = 77.595559
pos = SkyCoord(ra=ra, dec=dec, unit='deg')
```

## 2. Discover AllWISE Atlas images

+++

IRSA provides Simple Image Access (SIA) services for various datasets. A list of available datasets and their access URLs can be found at:

https://irsa.ipac.caltech.edu/ibe/sia.html

This tutorial uses SIA v2 for AllWISE Atlas images.

```{code-cell} ipython3
names = Irsa.list_collections(filter="allwise")
names

# We see from the resulting table that the dataset collection we are interested in is called "wise_allwise"
```

## 3. Search for images 

```{code-cell} ipython3
#get a table of all images within 1 arcsecond of our target position

dataset_name = names['collection'][0]
im_table = Irsa.query_sia(pos=(pos, 1 * u.arcsec), collection=dataset_name)
```

Inspect the table that is returned

```{code-cell} ipython3
im_table
```

```{code-cell} ipython3
im_table.colnames
```

```{code-cell} ipython3
# Let's look at the values in one of the columns
im_table['energy_bandpassname']
```

## 4.Locate and visualize an image of interest

+++

Let's search the image results for the W3 band image.

```{code-cell} ipython3
# You can put the URL from the column "access_url" into a browser to download the file. 
# Or you can work with it in Python, as shown below.
w3_mask = im_table['energy_bandpassname'] == 'W3'
w3_table = im_table[w3_mask]

# Lets look at the access_url of the first one:
image_url = w3_table['access_url'][0]
image_url
```

```{code-cell} ipython3
#Use Astropy to examine the header of the URL from the previous step.

hdulist = fits.open(image_url)
hdulist.info()
```

Download the image and open it in Astropy

```{code-cell} ipython3
## 7. Visualize a SPHEREx Spectral Image MEF using the Firefly Python Client.


# We will use the open-source astronomy data visualization software Firefly. 
# Firefly has a Python client.

#Open a Firefly viewer in a tab within jupyterlab.

fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")
fc = FireflyClient.make_lab_client()

# Visualize an image by sending its URL to the viewer.

fc.show_fits_image(file_input=image_url,
             plot_id="image",
             Title="Image"
             )

#Try use the interactive tools in the viewer to explore the data.

```

## 5. Extract a cutout and plot it

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

+++

## About this notebook

**Authors:** IRSA Data Science Team, including Troy Raen, Brigitta Sipőcz, Jessica Krick, Andreas Faisst, Jaladh Singhal, Vandana Desai, Dave Shupe

**Updated:** 16 February 2026

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.

**Runtime:** As of the date above, this notebook takes about 20 seconds to run to completion on a machine with 8GB RAM and 4 CPU.
This runtime is dependent on archive servers which means runtime will vary for users.

+++

## Citations

+++

If you use `astropy` for published research, please cite the authors. Follow these links for more information about citing `astropy`:

* [Citing `astropy`](https://www.astropy.org/acknowledging.html)

+++

Please include the following in any published material that makes use of the WISE data products:

"This publication makes use of data products from the Wide-field Infrared Survey Explorer, which is a joint project of the University of California, Los Angeles, and the Jet Propulsion Laboratory/California Institute of Technology, funded by the National Aeronautics and Space Administration."

Please also cite the dataset Digital Object Identifier (DOI): [10.26131/IRSA153](https://www.ipac.caltech.edu/doi/irsa/10.26131/IRSA153)
