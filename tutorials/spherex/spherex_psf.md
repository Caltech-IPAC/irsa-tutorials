---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
authors:
  - name: Vandana Desai
  - name: Jessica Krick
  - name: Andreas Faisst
  - name: Brigitta Sipőcz
  - name: Troy Raen
---

# Understanding and Extracting the PSF Extension in a SPHEREx Cutout

+++

## 1. Learning Goals

* Determine how pixels in a SPHEREx cutout map to the pixels in the parent SPHEREx spectral image.
* Understand the structure of the PSF extension in a SPHEREx cutout (which is the same as the PSF extension in the parent spectral image).
* Learn how to tell which version of the SPHEREx spectral image you are looking at, and how to interpret this information to obtain the correct PSF extension for the SPHEREx spectral images.
* Learn which plane in a SPHEREx cutout PSF extension cube most accurately describes the coordinates you are interested in.

+++

```{note}
This notebook is not intended for use of QR-1 data.
```

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
import copy
from packaging.version import Version

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

```{note}
SPHEREx data are also available via SIA which can provide a simpler interface for many queries, as demonstrated in {ref}`spherex-intro`.
An advantage of the method shown above is that it provides access to data immediately after ingestion (which occurs weekly) and is not subject to the same ~1 day delay as SIA.
```

For this example, we focus on the first one of the retrieved SPHEREx spectral images.

```{code-cell} ipython3
spectral_image_url = results['uri'][0]
print(spectral_image_url)
```

## 5. Read in a SPHEREx Cutout

Next, we use standard astropy tools to open the fits image and to read the different headers and data.
Transient read errors occur sometimes, so we'll catch those and retry a few times.

```{tip}
As we do below, you can use `hdul.info()` to print the list of FITS layers of the downloaded cutout.
```

```{code-cell} ipython3
# Max number of times to retry transient read errors.
max_retries = 3
for attempt in range(max_retries):
    try:
        # Read the data.
        with fits.open(spectral_image_url) as hdul:
            image_hdul = copy.deepcopy(hdul)
            cutout_header = hdul['IMAGE'].header
            psf_header = hdul['PSF'].header
            cutout = hdul['IMAGE'].data
            psfcube = hdul['PSF'].data
        break
    except (TimeoutError, urllib.error.HTTPError, http.client.IncompleteRead):
        if attempt == max_retries - 1:
            raise
        time.sleep(10 * (attempt + 1))
```

Let's examine the HDU list info.

```{code-cell} ipython3
hdul.info()
```

The downloaded SPHEREx image cutout contains 5 FITS layers, which are described in the [SPHEREx Explanatory Supplement](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR.pdf).
We focus in this example on the extensions `IMAGE` and `PSF`.
We have already loaded their data as well as their header.

```{code-cell} ipython3
psfcube.shape
```

The shape of the `psfcube` is output above.

```{note}
In the QR-2 data, the shape is (121,101,101), which corresponds to a grid of 11x11 PSF zones across the image.
The number of PSF zones may change in later versions of data products.
```

Each PSF has a size of 101x101 pixels.

```{note}
Remember that the PSFs are oversampled by a factor of 10.
This means that the actual size of the PSFs is about 10x10 SPHEREx pixels, which corresponds to about 60x60 arcseconds.
```

+++

Let's look at a small part of the PSF header to understand its format.

```{code-cell} ipython3
psf_header[0:30]
```

We confirm that the oversampling factor (`OVERSAMP`) is 10.
The PSFs are distributed in an even grid with NxM zones (in QR-2 data products it is N=M=11).
Each of the NxM PSFs is responsible for one of these zones.
The PSF header therefore includes the center position of these zones as well as the width of the zones.
These center coordinate are specified with `XCTR_i` and `YCTR_i`, respectively, where i = 1...(NxM).
The widths are specified with `XWID_i` and `YWID_i`, respectively, where again i = 1...(NxM).
The zones have approximately equal widths and are arranged in an even grid.
The size of the zones is sufficient to capture well the changes of the PSF size and structure with wavelength and spatial coordinates.

The goal of this tutorial now is to find the PSF corresponding to our input coordinates of interest.

+++

```{warning}
In the SPHEREx spectral image versions prior or equal to 6.5.5, there was a mismatch between the spatial layout of the PSF zones and the indexing of the PSF zones in the image header. This has now been fixed in versions 6.5.6 and beyond.

For more information about these changes, see the following webpage: [PSF Erratum](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/psfhdrerr.html)

**Users using the old versions will need to implement an extra step to update the image header. A function to update the header is given [in Section 5.1 below](#update-psf-header-function).**
```

Let's first check here if a header update is necessary. We can do that by printing the `VERSION` keyword in the header.

For comparing versions, we can use the Python-internal `Version()` function from the `packaging.version` package. Images that have already been reprocessed can have version names such as `6.5.4+psffix1` (which are superior to `6.5.4`, for example), and we can use `Version().local` to check for those.

```{code-cell} ipython3
this_version = Version(image_hdul['PRIMARY'].header["VERSION"])
contains_psffix1 = this_version.local is not None and "psffix1" in this_version.local
print(f"Current version is {this_version}")

if this_version <= Version("6.5.5") and not contains_psffix1:
    print("PSF header needs to be updated! -> Go to Section 5.1 :(")
else:
    print("PSF header is already up-to-date! -> Proceed to Section 6 :)")
```

If the version of the SPHEREx spectral image is less or equal than `6.5.5` and hasn't already been reprocessed, we will have to update the header. This is explained in Section 5.1. If the version is later than `6.5.5` or includes `"psffix1"`, the header is already updated and the PSF issue is fixed. In this case, proceed to Section 6 directly.

+++

(update-psf-header-function)=
### 5.1 Updating PSF Header (for SPHEREx Spectral Image versions $\leq$ 6.5.5)

+++

The function that can be used to update the header is shown below. The function
* first checks if a header update is necessary
* changes the PSF zone indexing and
* changes the version of the header such that it is consistent with the new released images

Note that this function can work as standalone function to process many images.

```{code-cell} ipython3
def update_psf_header(old_hdul):
    """
    Fix an old PSF FITS file header by rewriting only the per-plane header metadata
    so that plane k corresponds to x-fast ordering:
        k0 = iy * bins_x + ix

    The cube data are left untouched.

    Parameters
    ----------
    old_hdul : fits.HDUList
        Old SPHEREx Spectral Image HDUL

    Return
    ----------
    new_hdul : fits.HDUList
        New SPHEREx Spectral Image HDUL with updated PSF zone data in header and updated version number
    """

    VERSION_FIXED = Version("6.5.6")
    PSF_FIX_TAG = "psffix1"

    def psf_fix_applied(hdul) -> bool:
        """
        Return True if the PSF fix has been applied.

        Rules:
        - If the VERSION header is missing in the primary HDU, the fix is not applied.
        - If VERSION >= VERSION_FIXED, the fix is included in the software release.
        - Otherwise the local version tag (+...) must contain PSF_FIX_TAG.
        """
        header = hdul[0].header

        if "VERSION" not in header:
            return False

        v = Version(header["VERSION"])

        if v >= VERSION_FIXED:
            return True

        return v.local is not None and PSF_FIX_TAG in v.local

    if psf_fix_applied(old_hdul):
        return old_hdul

    ## Define some auxiliary functions -------
    def parse_ixiy_from_comment(comment):
        _zone_pat = re.compile(r"\((\d+)\s*,\s*(\d+)\)")
        m = _zone_pat.search(str(comment))
        if not m:
            raise ValueError(f"Could not parse zone indices from comment: {comment!r}")
        return int(m.group(1)), int(m.group(2))

    def infer_grid_shape_from_header_comments(hdr, nzone):
        max_ix = -1
        max_iy = -1

        for k1 in range(1, nzone + 1):
            key = f"XCTR_{k1}"
            if key not in hdr:
                raise KeyError(f"Missing required key: {key}")
            ix, iy = parse_ixiy_from_comment(hdr.comments[key])
            max_ix = max(max_ix, ix)
            max_iy = max(max_iy, iy)

        bins_x = max_ix + 1
        bins_y = max_iy + 1

        if bins_x * bins_y != nzone:
            raise ValueError(
                f"Inconsistent grid inferred from comments: "
                f"bins_x={bins_x}, bins_y={bins_y}, nzone={nzone}"
            )

        return bins_x, bins_y

    def collect_axis_values_by_zone(hdr, nzone):
        """
        Read the old header and collect unique x/y centers and widths by zone index
        labels found in the comments.

        This uses the old header only to recover the per-axis values for each ix, iy.
        It does NOT use the old plane ordering as truth.
        """
        x_center_by_ix = {}
        y_center_by_iy = {}
        x_width_by_ix = {}
        y_width_by_iy = {}

        for k1 in range(1, nzone + 1):
            ix, iy = parse_ixiy_from_comment(hdr.comments[f"XCTR_{k1}"])

            xck = f"XCTR_{k1}"
            yck = f"YCTR_{k1}"
            xwk = f"XWID_{k1}"
            ywk = f"YWID_{k1}"

            if xck in hdr:
                val = hdr[xck]
                if ix in x_center_by_ix and not np.isclose(x_center_by_ix[ix], val):
                    raise ValueError(
                        f"Inconsistent XCTR for ix={ix}: "
                        f"{x_center_by_ix[ix]} vs {val}"
                    )
                x_center_by_ix[ix] = val

            if yck in hdr:
                val = hdr[yck]
                if iy in y_center_by_iy and not np.isclose(y_center_by_iy[iy], val):
                    raise ValueError(
                        f"Inconsistent YCTR for iy={iy}: "
                        f"{y_center_by_iy[iy]} vs {val}"
                    )
                y_center_by_iy[iy] = val

            if xwk in hdr:
                val = hdr[xwk]
                if ix in x_width_by_ix and not np.isclose(x_width_by_ix[ix], val):
                    raise ValueError(
                        f"Inconsistent XWID for ix={ix}: "
                        f"{x_width_by_ix[ix]} vs {val}"
                    )
                x_width_by_ix[ix] = val

            if ywk in hdr:
                val = hdr[ywk]
                if iy in y_width_by_iy and not np.isclose(y_width_by_iy[iy], val):
                    raise ValueError(
                        f"Inconsistent YWID for iy={iy}: "
                        f"{y_width_by_iy[iy]} vs {val}"
                    )
                y_width_by_iy[iy] = val

        return x_center_by_ix, y_center_by_iy, x_width_by_ix, y_width_by_iy
    ## End defining some auxiliary functions --------

    ## Get Header
    extname = "PSF"
    hdu = old_hdul[extname]
    cube = np.asarray(hdu.data)
    hdr_in = hdu.header.copy()

    if cube.ndim != 3:
        raise ValueError(f"Expected 3D PSF cube, got shape {cube.shape}")

    nzone = cube.shape[0]
    bins_x, bins_y = infer_grid_shape_from_header_comments(hdr_in, nzone)

    print(f"Detected bins_x={bins_x}, bins_y={bins_y}, nzone={nzone}")

    x_center_by_ix, y_center_by_iy, x_width_by_ix, y_width_by_iy = collect_axis_values_by_zone(
        hdr_in, nzone
    )

    # Validate that all needed axis values were recovered
    missing = []
    for ix in range(bins_x):
        if ix not in x_center_by_ix:
            missing.append(f"x_center[{ix}]")
        if ix not in x_width_by_ix:
            missing.append(f"x_width[{ix}]")
    for iy in range(bins_y):
        if iy not in y_center_by_iy:
            missing.append(f"y_center[{iy}]")
        if iy not in y_width_by_iy:
            missing.append(f"y_width[{iy}]")

    if missing:
        raise ValueError(f"Missing axis metadata recovered from old header: {missing}")

    hdr_out = hdr_in.copy()

    # Rewrite only the per-plane metadata so plane k matches x-fast ordering.
    # plane k0 should correspond to:
    #   ix = k0 % bins_x
    #   iy = k0 // bins_x
    for k0 in range(nzone):
        ix = k0 % bins_x
        iy = k0 // bins_x
        k1 = k0 + 1

        hdr_out[f"XCTR_{k1}"] = (
            x_center_by_ix[ix],
            f"Center of x zone ({ix}, {iy})"
        )
        hdr_out[f"YCTR_{k1}"] = (
            y_center_by_iy[iy],
            f"Center of y zone ({ix}, {iy})"
        )
        hdr_out[f"XWID_{k1}"] = (
            x_width_by_ix[ix],
            f"Width of x zone ({ix}, {iy})"
        )
        hdr_out[f"YWID_{k1}"] = (
            y_width_by_iy[iy],
            f"Width of y zone ({ix}, {iy})"
        )

    # Optional but useful provenance note
    hdr_out["HISTORY"] = "Rewrote PSF per-plane zone metadata to x-fast ordering."
    hdr_out["HISTORY"] = "Cube plane data left unchanged."



    new_hdu = fits.ImageHDU(data=cube, header=hdr_out, name=hdu.name)

    ext_index = old_hdul.index_of(extname)
    new_hdul = fits.HDUList()
    for i, old in enumerate(old_hdul):
        if i == ext_index:
            new_hdul.append(new_hdu)
        else:
            new_hdul.append(old.copy())

    ## TO DO: UPDATE VERSION
    new_hdul['PRIMARY'].header["VERSION"] = new_hdul['PRIMARY'].header["VERSION"] + "+psffix1" # SET NEW VERSION HERE

    return(new_hdul)
```

We now run this function to create a new HDU list that we will use later.

```{code-cell} ipython3
new_image_hdul = update_psf_header(old_hdul=image_hdul)
```

Let's check if the version keywords was updated:

```{code-cell} ipython3
print(f"Old version: {image_hdul['PRIMARY'].header['VERSION']}")
print(f"Updated version: {new_image_hdul['PRIMARY'].header['VERSION']}")
```

Let's compare the new and old PSF headers to see the difference.

```{code-cell} ipython3
image_hdul['PSF'].header[22:40]
```

```{code-cell} ipython3
new_image_hdul['PSF'].header[22:40]
```

Now we have to update the variables we have set above.

```{code-cell} ipython3
cutout_header = new_image_hdul['IMAGE'].header
psf_header = new_image_hdul['PSF'].header
cutout = new_image_hdul['IMAGE'].data
psfcube = new_image_hdul['PSF'].data
```

With this fix, we are now ready to proceed!

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
        yctr[yplane] = val
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

**Updated:** 10 March 2026

**Contact:** Contact [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or problems.

**Runtime:** Approximately 30 seconds.
