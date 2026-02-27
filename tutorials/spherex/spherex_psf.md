---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.3
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Understanding and Extracting the PSF Extension in a SPHEREx Cutout

+++

## 1. Learning Goals

* Determine how pixels in a SPHEREx cutout map to the pixels in the parent SPHEREx spectral image.
* Understand the structure of the PSF extension in a SPHEREx cutout (which is the same as the PSF extension in the parent spectral image)
* Learn which plane in a SPHEREx cutout PSF extension cube most accurately describes the coordinates you are interested in.

+++

## 2. SPHEREx Overview

SPHEREx is a NASA Astrophysics Medium Explorer mission that launched in March 2025.
During its planned two-year mission, SPHEREx will obtain 0.75-5 micron spectroscopy over the entire sky, with deeper data in the SPHEREx Deep Fields.
SPHEREx data will be used to:

* **constrain the physics of inflation** by measuring its imprints on the three-dimensional large-scale distribution of matter,
* **trace the history of galactic light production** through a deep multi-band measurement of large-scale clustering,
* **investigate the abundance and composition of water and biogenic ices** in the early phases of star and planetary disk formation.

The community will also mine SPHEREx data and combine it with synergistic data sets to address a variety of additional topics in astrophysics.

More information is available in the [SPHEREx Explanatory Supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR.pdf).

+++

## 3. Imports

The following packages must be installed to run this notebook.

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install astropy numpy pyvo
```

```{code-cell} ipython3
import http.client
import re
import time
import urllib.error

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pyvo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

# The time it takes to read SPHEREx files can exceed
# astropy's default timeout limit. Increase it.
from astropy.utils.data import conf
conf.remote_timeout = 120
```

## 4. Get SPHEREx Cutout

We first obtain a SPHEREx cutout for a given coordinate of interest from IRSA archive.
For this we define a coordinate and a size of the cutout.
Both should be defined using `astropy` units.
The goal is to obtain the cutout and then extract the PSF corresponding to the coordinates of interest.

```{tip}
To learn more about how to access SPHEREx spectral images and how to download cutouts, we refer to the [SPHEREx Intro Tutorial](#spherex-intro) and the [SPHEREx Cutouts Tutorial](#spherex-cutouts).
```

```{code-cell} ipython3
ra = 305.59875000000005 * u.degree
dec = 41.14888888888889 * u.degree
size = 0.01 * u.degree
```

Once we defined the coordinates of interest and the size of the cutout, we run a TAP query to gather all SPHEREx spectral images that cover the coordinates.

```{code-cell} ipython3
# Define the service endpoint for IRSA's Table Access Protocol (TAP)
# so that we can query SPHEREx metadata tables.
service = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

# Define a query that will search the appropriate SPHEREx metadata tables
# for spectral images that cover the chosen coordinates of interest.
# Return the cutout data access URL and the time of observation.
# Sort by observation time.
query = f"""
SELECT
    'https://irsa.ipac.caltech.edu/' || a.uri || '?center={ra.value},{dec.value}d&size={size.value}' AS uri,
    p.time_bounds_lower
FROM spherex.artifact a
JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', {ra.value}, {dec.value}), p.poly)
ORDER BY p.time_bounds_lower
"""

# Execute the query and return as an astropy Table.
t1 = time.time()
results = service.search(query)
print("Time to do TAP query: {:2.2f} seconds.".format(time.time() - t1))
print("Number of images found: {}".format(len(results)))
```

:::{note}
SPHEREx data are also available via SIA which can provide a simpler interface for many queries, as demonstrated in {ref}`spherex-intro`.
An advantage of the method shown above is that it provides access to data immediately after ingestion (which occurs weekly) and is not subject to the same ~1 day delay as SIA.
:::

For this example, we focus on the first one of the retrieved SPHEREx spectral images.

```{code-cell} ipython3
spectral_image_url = results['uri'][0]
print(spectral_image_url)
```

## 5. Read in a SPHEREx Cutout

Next, we use standard astropy tools to open the fits image and to read the different headers and data.

```{tip}
As we do below, you can use `hdul.info()` to print the list of FITS layers of the downloaded cutout.
```

```{code-cell} ipython3
Max number of times to retry transient read errors.
max_retries = 3
for attempt in range(max_retries):
    try:
        with fits.open(spectral_image_url) as hdul:
            hdul.info()
            cutout_header = hdul['IMAGE'].header
            psf_header = hdul['PSF'].header
            cutout = hdul['IMAGE'].data
            psfcube = hdul['PSF'].data
        break
    except (urllib.error.HTTPError, http.client.IncompleteRead):
        if attempt == max_retries - 1:
            raise
        time.sleep(10 * (attempt + 1))
```

The downloaded SPHEREx image cutout contains 5 FITS layers, which are described in the [SPHEREx Explanatory Supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR.pdf).
We focus in this example on the extensions `IMAGE` and `PSF`.
We have already loaded their data as well as their header.

```{code-cell} ipython3
psfcube.shape
```

The shape of the `psfcube` is (121,101,101).
This corresponds to a grid of 11x11 PSFs across the image, each of them of the size 101x101 pixels.

```{note}
Remember that the PSFs are oversampled by a factor of 10.
This means that the actual size of the PSFs is about 10x10 SPHEREx pixels, which corresponds to about 60x60 arcseconds.
```

+++

Let's look at a small part of the PSF header to understand its format:

```{code-cell} ipython3
psf_header[22:40]
```

We confirm that the oversampling factor (`OVERSAMP`) is 10.
The PSFs are distributed in an even grid with 11x11 zones.
Each of the 121 PSFs is responsible for one of these zones.
The PSF header therefore includes the center position of these zones as well as the width of the zones.
These center coordinate are specified with `XCTR_i` and `YCTR_i`, respectively, where i = 1...121.
The widths are specified with `XWID_i` and `YWID_i`, respectively, where again i = 1...121.
The zones have approximately equal widths and are arranged in an even grid.
The size of the zones is sufficient to capture well the changes of the PSF size and structure with wavelength and spatial coordinates.

The goal of this tutorial now is to find the PSF corresponding to our input coordinates of interest.

+++

## 6. Determine the Pixel Location on the Parent SPHEREx Image

To identify the zone which covers the coordinates of interest, we first need to translate these coordinates to the pixel coordinates on the parent large SPHEREx image from which the cutout was created.

We do this by first determining the pixel (x,y) coordinates of our coordinates of interest on the cutout itself.

```{code-cell} ipython3
wcs = WCS(cutout_header)
xpix_cutout, ypix_cutout = wcs.world_to_pixel(SkyCoord(ra=ra, dec=dec))

print(f"Pixel values of coordinates of interest on cutout image: x = {xpix_cutout}, y = {ypix_cutout}")
```

Next, we use the `CRPIX1A` and `CRPIX1A` header keywords (which describe the center of the cutout on the parent SPHEREx image) to shift the (x,y) coordinates of input to the parent SPHEREx image.

```{code-cell} ipython3
crpix1a = cutout_header["CRPIX1A"]
crpix2a = cutout_header["CRPIX2A"]

xpix_orig = 1 + xpix_cutout - crpix1a
ypix_orig = 1 + ypix_cutout - crpix2a

print(f"Pixel values of coordinates of interest on parent SPHEREx image: x = {xpix_orig}, y = {ypix_orig}")
```

## 7. Determine the PSF Corresponding to Coordinates of Interest

Since we now know the (x,y) pixel values of the coordinates of interest on the parent SPHEREx image, we can identify the PSF zone.
In the following we first extract the zone pixel coordinates from the `XCTR_*` and `YCTR_*` keys in the PSF header.

```{code-cell} ipython3
xctr = {}
yctr = {}

for key, val in psf_header.items():
    # Look for keys like XCTR* or YCTR*
    xm = re.match(r'(XCTR*)', key)
    if xm:
        xplane = int(key.split("_")[1])
        xctr[xplane] = val
    ym = re.match(r'(YCTR*)', key)
    if ym:
        yplane = int(key.split("_")[1])
        yctr[xplane] = val
```

Check that we got all of them!

```{code-cell} ipython3
len(xctr) == len(yctr)
```

Make a nice table so we can easily search for the distance between zone center and coordinates of interest.

```{code-cell} ipython3
tab = Table(names=["zone_id" , "x" , "y"], dtype=[int, float, float])
for zone_id in xctr.keys():
    tab.add_row([zone_id , xctr[zone_id] , yctr[zone_id]])
```

Once we have created this dictionary with zone pixel coordinates, we can simply search for the closest zone center to the coordinates of interest.
For this we first add the distance between zone center coordinates and coordinates of interest to the table. (Note that the x,y coordinates of the PSF zone centers are in 1,1 convention, therefore we have to subtract 1 pixels.)

```{code-cell} ipython3
tab["distance"] = np.sqrt((tab["x"]-1 - xpix_orig)**2 + (tab["y"]-1 - ypix_orig)**2)
```

Then we can sort the table and pick the closest zone to coordinates of interest.

```{code-cell} ipython3
tab.sort("distance")

psf_cube_plane = tab[0]["zone_id"]
distance_min = tab[0]["distance"]

print(f"The PSF zone corresponding to coordinates of interest is {psf_cube_plane} with a distance of {distance_min} pixels")
```

## 8. Extract and Show the PSF

Now that we know which zone corresponds to coordinates of interest, we can extract it and plot it.

```{code-cell} ipython3
psf = psfcube[psf_cube_plane-1]

fig = plt.figure(figsize=(5, 5))
ax1 = fig.add_subplot(1, 1, 1)

ax1.imshow(psf)

plt.show()
```

## 9. Using the SPHEREx PSF in Forward Modeling (e.g., Tractor)

The PSF returned by this notebook is oversampled relative to the native SPHEREx detector pixel grid. 
This is intentional: the PSF is evaluated on a fine sub-pixel grid so that it can represent different intra-pixel source positions accurately.

Tools such as Tractor do not expect an oversampled PSF directly. 
Instead, they require a PSF that is pixel-integrated at the native detector resolution and evaluated at the correct sub-pixel phase of the source. 
If you pass the oversampled PSF directly into Tractor without resampling, the effective PSF width and normalization will be incorrect, which can lead to systematic differences relative to the SPHEREx Spectrophotometry Tool.

To use this PSF for forward modeling or fitting, you must:
1. Shift the oversampled PSF to the source’s sub-pixel position,
2. Downsample (integrate) it onto the native SPHEREx pixel grid, and
3. Normalize the resulting PSF before passing it to Tractor.

+++

## Acknowledgements

- [Caltech/IPAC-IRSA](https://irsa.ipac.caltech.edu/)

## About this notebook

**Authors:** IRSA Data Science Team, including Vandana Desai, Andreas Faisst, Brigitta Sipőcz, Troy Raen

**Updated:** 24 October 2025

**Contact:** Contact [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.

**Runtime:** Approximately 30 seconds.
