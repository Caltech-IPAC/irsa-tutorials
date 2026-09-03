---
authors:
- name: Jaladh Singhal
- name: Jessica Krick
- name: Vandana Desai
- name: IRSA Team
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.3
kernelspec:
  display_name: irsa-tutorials
  language: python
  name: python3
---

# Using Firefly to Explore OpenUniverse2024 Simulated Roman and Rubin Images

+++

## Learning Goals

By the end of this tutorial, you will:

- Learn how to access cloud-hosted Roman and Rubin simulated images.

- Learn how to find the Roman and Rubin images covering a sky position with the IRSA Simple Image Access (SIA) service.

- Learn how to launch an interactive Firefly instance inside JupyterLab.

- Learn how to use the Firefly Jupyterlab extension to visualize cloud-hosted simulated images, overplot ds9 regions, overplot catalogs in Parquet format, and create 3 color images.

+++

## Introduction

The purpose of this tutorial is to become familiar with the simulated Roman and Rubin data published through OpenUniverse2024, and to become familiar with the Firefly JupyterLab Extension for visualizing astronomical data products.

OpenUniverse2024 is a project to simulate spatially overlapping imaging surveys to be carried out by the Nancy Grace Roman Telescope and the Vera C. Rubin Observatory. The simulations were carried out on Argonne's Theta cluster and consist of:

- The LSST ELAIS-S1 Deep Drilling Field (DDF)
- The Roman Time-Domain Survey (TDS) shifted to overlap the ELAIS region and LSST DDF
- Overlapping LSST Wide-Fast-Deep (WFD) survey (with rolling cadence)
- Overlapping Roman Wide-Area Survey (WAS) in the same region
- A deep-field calibration region of the Roman WAS in the same region

This tutorial works with the full simulation, which covers the entire survey footprint. More information about the dataset can be found at [IRSA's holding of this dataset](https://irsa.ipac.caltech.edu/data/theory/openuniverse2024/overview.html), and the [OpenUniverse2024 paper](https://doi.org/10.1093/mnras/staf1833) describes how the simulation was produced.

Firefly is an open-source web-based UI library for astronomical data archive access and visualization developed at Caltech and used by multiple space- and ground-based astrophysics archives. More information on Firefly can be found [here](https://github.com/Caltech-IPAC/firefly/blob/dev/README.md).

In addition to being used to make web applications, Firefly can be used from Python. More information on Firefly Python client can be found [here](https://caltech-ipac.github.io/firefly_client/usage/index.html).

The Firefly JupyterLab Extension makes it particularly easy to use Firefly to efficiently visualize cloud-hosted astronomical data using JupyterLab instances running locally or on cloud. More information on Firefly JupyterLab Extension can be found [here](https://github.com/Caltech-IPAC/jupyter_firefly_extensions/blob/master/README.md).

If you are new to OpenUniverse2024, the [Quickstart](openuniverse2024_quickstart) tutorial introduces the directory layout, the parquet catalogs, and the SIA image search that this notebook builds on.

+++

## Imports

- astropy.io.fits for accessing FITS files
- numpy for numerical computing
- matplotlib.pyplot for creating static visualizations of FITS images
- matplotlib.patches for annotating visualizations of FITS images
- astropy.wcs for dealing with astronomical world coordinate systems
- astropy.units for dealing with astronomical units
- astropy.coordinates.SkyCoord for dealing with astronomical coordinates
- astroquery.ipac.irsa.Irsa for asking IRSA which images cover a position
- hpgeom for converting a sky position into a HEALPix region index
- json for reading the cloud location returned by the image search
- firefly_client.FireflyClient for using the Firefly python client
- astropy.nddata.Cutout2D for making image cutouts
- reproject.reproject_interp for rotating cutouts so both telescopes' panels point the same way
- io.BytesIO for writing a fits file to an in-memory stream
- warnings for silencing the harmless FITS WCS warnings

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy astropy s3fs matplotlib firefly_client astroquery hpgeom reproject
```

```{code-cell} ipython3
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa
import hpgeom
import json
from firefly_client import FireflyClient
from astropy.nddata import Cutout2D
from reproject import reproject_interp
from io import BytesIO
import warnings

# Filter out the FITSFixedWarning, which is consequenceless and gets thrown every time you deal with a WCS
# in a Roman openuniverse simulated image using astropy.
warnings.simplefilter('ignore', category=wcs.FITSFixedWarning)
```

## 1. Learn where the OpenUniverse2024 data are hosted in the cloud

The OpenUniverse2024 data are hosted in the cloud via Amazon Web Services (AWS). To access these data, you need to create a client to read data from Amazon's Simple Storage Service (s3) buckets, and you need to know some information about those buckets. OpenUniverse2024 contains simulations of the Roman Wide-Area Survey (WAS) and the Roman Time Domain Survey (TDS). In this tutorial, we will focus on the WAS.

We ask IRSA's Simple Image Access (SIA) service to locate both the Roman and Rubin images by coordinate.

```{code-cell} ipython3
BUCKET_NAME = "nasa-irsa-simulations"
OU_PREFIX = "openuniverse2024"
TRUTH_FILES_PATH = f"{OU_PREFIX}/roman/full/roman_rubin_cats_v1.1.2_faint"
```

```{code-cell} ipython3
# Point the astroquery IRSA client at the simulated-data services, which are
# separate from the ones serving IRSA's observed data.
Irsa.sia_url = "https://irsa.ipac.caltech.edu/simulated/SIA"
Irsa.tap_url = "https://irsa.ipac.caltech.edu/simulated/TAP"

OU_ROMAN_SIA_COLLECTION = "simulated_roman_openuniverse2024"
OU_RUBIN_SIA_COLLECTION = "simulated_rubin_openuniverse2024"

# A small radius is all we need, since we only want images containing a given point.
SEARCH_RADIUS = 1 * u.arcsec
```

## 2. Roman Images

The Nancy Grace Roman Space Telescope will carry out a wide-area survey (WAS) in the near infrared. OpenUniverse2024 simulates that survey as individual exposures, each covering one detector's worth of sky. Bands include F184, H158, J129, K213, W146, Y106.

The survey stores those exposures by pointing and detector rather than by sky position, so in this section we define a function that asks the image search which exposure covers a position, and then reads it.

+++

### Define a function that retrieves a Roman simulated image given a sky position and filter.

The Roman collection holds both the Wide-Area Survey and the Time Domain Survey, so we keep only the WAS exposures and then only the band we want. Many exposures cover any given position, so we sort by observation time and take the earliest.

Having found the file, we prefix its path with `s3://` and open it with astropy. The science image sits in the first extension, and we take the WCS from the same header as an `astropy.wcs.WCS` object. The function returns a dictionary of both, along with the header itself so we can look at it later, and the path so that later cells can reuse it instead of searching a second time.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def get_roman_image(coord, filter):
    """Read the Roman WAS exposure covering a position in a given filter.

    Parameters
    ----------
    coord : `~astropy.coordinates.SkyCoord`
        Sky position to search around.
    filter : str
        Roman bandpass name, e.g. ``'H158'``.

    Returns
    -------
    dict
        Keys ``'data'`` (`~numpy.ndarray`), ``'wcs'`` (`~astropy.wcs.WCS`),
        ``'header'`` (`~astropy.io.fits.Header`) and ``'fpath'`` (str, the S3
        path of the exposure).

    Raises
    ------
    ValueError
        If no WAS exposure in that filter covers ``coord``.
    """
    results = Irsa.query_sia(pos=(coord, SEARCH_RADIUS),
                             collection=OU_ROMAN_SIA_COLLECTION)

    in_was = results[['WAS_simple_model' in str(r['obs_id']) for r in results]]
    in_band = in_was[[str(r['energy_bandpassname']).strip() == filter
                      for r in in_was]]
    if len(in_band) == 0:
        raise ValueError(f"No Roman WAS {filter}-band images cover "
                         f"{coord.to_string('hmsdms')}")
    
    in_band.sort('t_min')
    selected_image = in_band[0]

    cloud_info = json.loads(selected_image['cloud_access'])['aws']
    image_s3_fpath = f"{cloud_info['bucket_name']}/{cloud_info['key']}"

    with fits.open(f"s3://{image_s3_fpath}", fsspec_kwargs={"anon": True}) as hdul:
        # retrieve science data, which sits in the first extension
        image_data = hdul[1].data

        # make wcs using header
        image_wcs = wcs.WCS(hdul[1].header)

        return {'data': image_data, 'wcs': image_wcs, 'header': hdul[1].header,
                'fpath': image_s3_fpath}
```

### Inspect a simulated Roman image

+++

Choose a filter and a position within the survey footprint

```{code-cell} ipython3
coord = SkyCoord(ra=9.6205000, dec=-44.0641694, unit="deg")
filter_roman = 'H158' #F184, H158, J129, K213, W146, and Y106 are available
```

Retrieve the data and header information from the simulated Roman image corresponding to the chosen position and filter.

```{code-cell} ipython3
image_roman = get_roman_image(coord, filter_roman)
```

### Understand the size of a simulated Roman image.

```{code-cell} ipython3
# Number of pixels (Y, X)
image_roman['data'].shape
```

```{code-cell} ipython3
# Pixel size (scale Y, scale X) [degrees/pixel]
image_roman['wcs'].proj_plane_pixel_scales()
```

```{code-cell} ipython3
# Image size (FOV Y, FOV X)
[(num * size).to('arcsec') for num, size in zip(
    image_roman['data'].shape, image_roman['wcs'].proj_plane_pixel_scales())]
```

The field of view of a Roman exposure is ~450 arcsec, at a pixel scale of ~0.11 arcsec.

+++

### Look at everything else the image file contains

So far we have read the science image, but each exposure carries more than that. Listing the extensions shows the full picture.

```{code-cell} ipython3
with fits.open(f"s3://{image_roman['fpath']}",
               fsspec_kwargs={"anon": True}) as hdul:
    hdul.info()
```

Each file holds a primary header carrying the observation metadata but no pixels, followed by three 4088x4088 planes:

SCI = the simulated science image, in electrons per second

ERR = the estimated uncertainty on each science pixel

DQ = the data quality mask, flagging pixels that should not be trusted

+++

### Use the WCS from the Roman simulated header to convert the specified coordinate into a pixel position.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def coord_to_xy(w, coord):
    """Convert a sky position into (x, y) pixel indices.

    Parameters
    ----------
    w : `~astropy.wcs.WCS`
        WCS of the image the indices refer to.
    coord : `~astropy.coordinates.SkyCoord`
        Sky position to convert.

    Returns
    -------
    tuple of int
        Pixel indices as ``(x, y)``.
    """
    return w.world_to_array_index(coord)[::-1] #reverse since 0th axis is y, 1st axis is x
```

```{code-cell} ipython3
coord_arr_idx = coord_to_xy(image_roman['wcs'], coord)
coord_arr_idx
```

### Use matplotlib imshow to create a static visualization of the Roman simulated image and overplot the selected position.

Each telescope points independently, so neither survey's images arrive with north pointing up: this Roman exposure is turned about 170 degrees from north, and the Rubin visit in the next section about 160 degrees.
Every figure in this notebook is therefore resampled onto a north-up grid before plotting, so that all of them, and both telescopes, can be read the same way.
The target grid is a plain gnomonic (TAN) one, so the SIP distortion carried by the exposure headers is resampled away along with the rotation.
The rotation leaves the frame tilted inside a slightly larger box, with empty wedges in the corners.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def stretch_color(data, clipPercent):
    """Percentile clip limits for displaying an image.

    Parameters
    ----------
    data : `~numpy.ndarray`
        Image values; NaNs are ignored.
    clipPercent : float
        Percentage clipped off each end of the distribution.

    Returns
    -------
    tuple of float
        Lower and upper limits, for passing to ``clim``.
    """
    # nan-aware: rotating to north up leaves empty corners
    return np.nanpercentile(data, (0 + clipPercent, 100 - clipPercent))


def north_up(data, image_wcs, size, binning=1):
    """Resample onto a north-up grid, keeping the original pixel scale.

    Parameters
    ----------
    data : `~numpy.ndarray`
        Image to resample.
    image_wcs : `~astropy.wcs.WCS`
        WCS of ``data``.
    size : `~astropy.units.Quantity`
        Angular width of the output grid.
    binning : int, optional
        Sample this many input pixels per output pixel. Values greater than 1
        are enough for displaying a whole exposure and much quicker than
        resampling every pixel.

    Returns
    -------
    `~numpy.ndarray`
        The resampled image.
    `~astropy.wcs.WCS`
        WCS of the resampled image.
    """
    scale = abs(image_wcs.proj_plane_pixel_scales()[0].to(u.deg).value) * binning
    npix = int(round(size.to(u.deg).value / scale))

    ny, nx = data.shape
    center = image_wcs.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)

    target = wcs.WCS(naxis=2)
    target.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    target.wcs.crpix = [(npix + 1) / 2, (npix + 1) / 2]
    target.wcs.crval = [center.ra.deg, center.dec.deg]
    target.wcs.cdelt = [-scale, scale]

    return reproject_interp((data, image_wcs), target, shape_out=(npix, npix))[0], target


def full_frame_north_up(image, binning=4):
    """North-up view of a whole exposure, for display only.

    Parameters
    ----------
    image : dict
        Image dictionary with ``'data'`` and ``'wcs'`` keys, as returned by
        `get_roman_image` or `get_rubin_image`.
    binning : int, optional
        Sample this many input pixels per output pixel.

    Returns
    -------
    `~numpy.ndarray`
        The resampled image.
    `~astropy.wcs.WCS`
        WCS of the resampled image.
    """
    ny, nx = image['data'].shape
    fov = nx * image['wcs'].proj_plane_pixel_scales()[0].to(u.arcsec)
    return north_up(image['data'], image['wcs'], 1.2 * fov, binning=binning)
```

```{code-cell} ipython3
roman_display, roman_display_wcs = full_frame_north_up(image_roman)

plt.imshow(roman_display, origin='lower', 
           clim=stretch_color(roman_display, 1)
           )

plt.plot(*coord_to_xy(roman_display_wcs, coord), '+', mfc='none', mec='r',
         markersize=15, markeredgewidth=2)
plt.axis('off')
```

## 3. Rubin Images

OpenUniverse2024 includes simulated Rubin images in the following filters: u, g, r, i, z, y. As with Roman, these are stored by visit and detector rather than by sky position, so we again let the SIA service tell us which images cover our position. In this section, we define functions that retrieve a Rubin image for a chosen position and filter, returning data in the same structure as the function we defined above for Roman.

+++

### Retrieve Rubin images

The image search returns one row per image, with the cloud location of each tucked inside a `cloud_access` JSON string. Many visits cover any given position, so we sort by observation time and take the earliest.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def get_rubin_image_fpaths(coord, filters):
    """Find the Rubin images covering a position in each of several filters.

    Parameters
    ----------
    coord : `~astropy.coordinates.SkyCoord`
        Sky position to search around.
    filters : list of str
        Rubin bandpass names, e.g. ``['g', 'r', 'i']``.

    Returns
    -------
    dict
        Filter name to S3 path of the image covering ``coord``.

    Raises
    ------
    ValueError
        If any requested filter has no image covering ``coord``.
    """
    # One search covers every band, so ask once and sort the results out by band
    # afterwards rather than searching again for each one.
    results = Irsa.query_sia(pos=(coord, SEARCH_RADIUS),
                             collection=OU_RUBIN_SIA_COLLECTION)

    fpaths = {}
    for filter in filters:
        # Band names come back like "r_57", so compare only the part before the underscore.
        in_band = results[[str(r['energy_bandpassname']).split('_')[0] == filter
                           for r in results]]
        if len(in_band) == 0:
            raise ValueError(f"No Rubin {filter}-band images cover "
                             f"{coord.to_string('hmsdms')}")

        in_band.sort('t_min')
        cloud_info = json.loads(in_band['cloud_access'][0])['aws']
        fpaths[filter] = f"{cloud_info['bucket_name']}/{cloud_info['key']}"

    return fpaths
```

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def get_rubin_image(coord, filter):
    """Read the Rubin image covering a position in a given filter.

    Parameters
    ----------
    coord : `~astropy.coordinates.SkyCoord`
        Sky position to search around.
    filter : str
        Rubin bandpass name, e.g. ``'r'``.

    Returns
    -------
    dict
        Keys ``'data'`` (`~numpy.ndarray`), ``'wcs'`` (`~astropy.wcs.WCS`) and
        ``'fpath'`` (str, the S3 path of the image).
    """
    image_s3_fpath = get_rubin_image_fpaths(coord, [filter])[filter]

    with fits.open(f"s3://{image_s3_fpath}", fsspec_kwargs={"anon": True}) as hdul:
        # retrieve science data, which sits in the first extension
        image_data = hdul[1].data

        # make wcs using header
        image_wcs = wcs.WCS(hdul[1].header)

        # hand back the path too, so later cells can reuse it instead of
        # asking the image search for it a second time
        return {'data': image_data, 'wcs': image_wcs, 'fpath': image_s3_fpath}
```

### Inspect a simulated Rubin image

+++

Choose a filter and retrieve the data and header information from the simulated Rubin image covering our position.

```{code-cell} ipython3
filter_rubin = 'r'
image_rubin = get_rubin_image(coord, filter_rubin)
```

### Understand the size of a simulated Rubin image.

```{code-cell} ipython3
# Number of pixels (Y, X)
image_rubin['data'].shape
```

```{code-cell} ipython3
# Pixel size (scale Y, scale X) [degrees/pixel]
image_rubin['wcs'].proj_plane_pixel_scales()
```

```{code-cell} ipython3
# Image size (FOV Y, FOV X)
[(num * size).to('arcsec') for num, size in zip(
    image_rubin['data'].shape, image_rubin['wcs'].proj_plane_pixel_scales())]
```

The field of view of a Rubin detector is ~800 arcsec, at a pixel scale of ~0.2 arcsec.
A single Rubin detector therefore covers nearly twice the sky of one Roman exposure, at roughly twice the pixel scale.

+++

### Use matplotlib imshow to create a static visualization of the Rubin simulated image and overplot the selected position.

The two full-frame figures do not cover the same patch of sky. Each is centered on its own exposure and sized to its own field of view, so the Rubin panel is nearly twice as wide on the sky as the Roman one. What they share is the marked position, which falls inside both.

```{code-cell} ipython3
rubin_display, rubin_display_wcs = full_frame_north_up(image_rubin)

plt.imshow(rubin_display, origin='lower', 
           clim=stretch_color(rubin_display, 1)
           )

plt.plot(*coord_to_xy(rubin_display_wcs, coord), '+', mfc='none', mec='r',
         markersize=15, markeredgewidth=2)
plt.axis('off')
```

## 4. Compare simulated Roman and Rubin cutouts

+++

### Choose cutout size

```{code-cell} ipython3
cutout_size = 50*u.arcsec
```

### Create the cutouts

```{code-cell} ipython3
# Cut oversized: the rotation below needs data in the corners to draw from.
cutout_roman = Cutout2D(image_roman['data'], coord, size=1.5*cutout_size, wcs=image_roman['wcs'])
cutout_rubin = Cutout2D(image_rubin['data'], coord, size=1.5*cutout_size, wcs=image_rubin['wcs'])
```

### Rotate both cutouts so that north points up

The full frames above were rotated for display; here we do the same to the cutouts, using the same helper. The two telescopes are turned about 10 degrees apart from each other, so without this the panels would not line up. Each keeps its own pixel scale, so Roman still shows the finer sampling.

```{code-cell} ipython3
# Rotating a same-size grid would leave empty corners, so each cutout above was cut
# oversized and is resampled down to the size we actually want.
roman_data, roman_wcs = north_up(cutout_roman.data, cutout_roman.wcs, cutout_size)
rubin_data, rubin_wcs = north_up(cutout_rubin.data, cutout_rubin.wcs, cutout_size)
```

### Plot static side-by-side comparisons of the cutouts using matplotlib

```{code-cell} ipython3
fig, axs = plt.subplots(1, 2, figsize=(12, 6))


axs[0].imshow(roman_data, origin='lower', 
              clim=stretch_color(roman_data, .5)
              )
axs[0].set_title(f"ROMAN in filter {filter_roman}")

axs[1].imshow(rubin_data, origin='lower', 
              clim=stretch_color(rubin_data, .5)
              )
axs[1].set_title(f"RUBIN in filter {filter_rubin}")

fig.suptitle(f"Cutouts at ({coord.ra}, {coord.dec}) with {cutout_size} size", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.97])
```

Both panels cover the same patch of sky, at the same size and with north up in each. The difference is stark: the simulated Roman exposure resolves a field full of galaxies, while the same patch in the simulated Rubin image is close to blank. A single Rubin visit is short, so only the brightest sources here rise above its noise, and from the ground its point spread function is broad enough to blend close neighbors together. Roman's finer pixels and sharper point spread function keep those sources separate. Rubin reaches its full depth by stacking many visits, which is not what we are looking at.

## 5. Use Firefly to interactively explore the images

Static cutouts only help once you know where to look. **Choosing a position worth comparing means exploring the image interactively**, and for that we hand it to Firefly.

+++

### Launch and initialize Firefly
There are two ways to initialize a Firefly client from Python, depending on whether you're running the notebook in JupyterLab or not. Assuming you have `jupyter-firefly-extensions` set up in your environment as explained [here](https://github.com/Caltech-IPAC/jupyter_firefly_extensions/blob/master/README.md), you can use `make_lab_client()` in JupyterLab, which will open the Firefly viewer in a new tab within the Lab. Otherwise, you can use `make_client()` in a Jupyter Notebook (or even a Python shell), which will open the Firefly viewer in a new web browser tab.

You also need a Firefly server to communicate with your Firefly Python client. In this notebook, we use a public Firefly server: https://irsacloud.ipac.caltech.edu/firefly. However, you can also run a local Firefly server via a [Firefly Docker image](https://hub.docker.com/r/ipac/firefly) and access it at `http://localhost:8080/firefly`. The URL of the Firefly server is read by both `make_client()` and `make_lab_client()` through the environment variable `FIREFLY_URL`. However, `make_client()` also allows you to pass the URL directly as the `url` parameter.

```{code-cell} ipython3
# Uncomment when using within Jupyter Lab with jupyter_firefly_extensions installed
# fc = FireflyClient.make_lab_client()

# Uncomment for contexts other than above 
fc = FireflyClient.make_client(url="https://irsacloud.ipac.caltech.edu/firefly")

fc.reinit_viewer() # to clean the state, if this cell ran earlier
```

### Define a function that returns the URL for a given S3 filepath
Since the OpenUniverse2024 data are in a public S3 bucket, any file can also be reached over HTTPS. Firefly loads images and catalogs by URL rather than by S3 path, so we use this helper throughout the rest of the notebook.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def https_url(s3_fpath):
    """Convert an S3 path into a public HTTPS URL.

    Parameters
    ----------
    s3_fpath : str
        Path of the form ``'bucket/key'``.

    Returns
    -------
    str
        HTTPS URL for the same file.
    """
    if s3_fpath is None:
        return None
    
    s3_fpath_without_bucket = s3_fpath.split('/', 1)[1]
    return f"https://{BUCKET_NAME}.s3.amazonaws.com/{s3_fpath_without_bucket}"
```

Let's generate URL for the Rubin image we plotted above. Clicking on the returned URL will allow you to download this image locally.

```{code-cell} ipython3
image_s3_fpath_rubin = image_rubin['fpath']
https_url(image_s3_fpath_rubin)
```

### Send the simulated Rubin & Roman images to Firefly using show_fits.

For displaying the FITS image of the Rubin & Roman visits in Firefly, we use [`show_fits_image`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_fits_image):

```{code-cell} ipython3
image_ff_id_rubin = 'rubin-image-filter-r'
fc.show_fits_image(file_input=https_url(image_s3_fpath_rubin),
                   plot_id=image_ff_id_rubin,
                   Title=f"Rubin {filter_rubin}"
             )
```

```{code-cell} ipython3
image_ff_id_roman = 'roman-image-filter-h158'
fc.show_fits_image(file_input=https_url(image_roman['fpath']),
                   plot_id=image_ff_id_roman,
                   Title=f"Roman {filter_roman}"
             )
```

Now that both images are on display, let's align & lock them by WCS using [`align_images`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.align_images) method so that panning or zooming one moves the other in sync.

```{code-cell} ipython3
fc.align_images(lock_match=True)
```

### Use ds9 region syntax to outline the Roman exposure on the interactive display

The Firefly client includes several methods related to controlling ds9 region overlays. To 
overlay a region layer on the loaded FITS images, we can use [`overlay_region_layer`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.overlay_region_layer).

Region data is defined in ds9 region syntax that can be found [here](https://ds9.si.edu/doc/ref/region.html).

The Rubin visit on display covers more sky than the Roman exposure does, so it is worth seeing where the Roman data actually falls. The corners of the Roman exposure come straight from its WCS, and we draw them as a polygon.

```{code-cell} ipython3
roman_footprint = image_roman['wcs'].calc_footprint()
roman_footprint
```

```{code-cell} ipython3
# outline the Roman exposure as a polygon
corner_coords = ' '.join(f'{ra}d {dec}d' for ra, dec in roman_footprint)
roman_regions = [f'icrs;polygon {corner_coords} # text={{Roman {filter_roman}}} color=yellow']

roman_regions_id = 'roman_regions'
# Name the plot explicitly. Without plot_id the layer goes to whichever plot is
# active, which is not necessarily this one if cells are run out of order.
fc.overlay_region_layer(region_data=roman_regions,
                        title='Roman exposure', 
                        region_layer_id=roman_regions_id,
                        plot_id=image_ff_id_rubin)
```

## 6. Use Firefly to overlay sources on the images
Let's inspect the properties of sources within the Roman footprint we outlined above. For this we will use the input truth files present in S3 bucket.

OpenUniverse2024 includes the input truth files that were used to create the simulated images. These files are in Parquet and HDF5 format, and include information about the properties of galaxies, stars, and transients.

The catalogs are split by HEALPix sky region (nside=32, RING ordering) and the region index appears in each filename, so we can convert our sky position into that index with [`hpgeom`](https://hpgeom.readthedocs.io/en/latest/) rather than hunting through the directory listing. This is the same approach the [Quickstart](openuniverse2024_quickstart) tutorial takes.

Each region comes as two separate catalogs: a point source catalog holding the stars, and a galaxy catalog. We load both below and overlay each on the image.

```{code-cell} ipython3
region = hpgeom.angle_to_pixel(32, coord.ra.deg, coord.dec.deg, lonlat=True, nest=False)
print(f"HEALPix region for our position: {region}")
```

```{code-cell} ipython3
# Catalog table of star properties (in parquet format)
pointsource_cat_path = f"{BUCKET_NAME}/{TRUTH_FILES_PATH}/pointsource_{region}.parquet"
pointsource_cat_path
```

```{code-cell} ipython3
# Catalog table of galaxy properties (in parquet format)
galaxy_cat_path = f"{BUCKET_NAME}/{TRUTH_FILES_PATH}/galaxy_{region}.parquet"
galaxy_cat_path
```

### Use Firefly's show_table to overlay the catalogs on the interactive image

Each catalog spans a whole HEALPix region, which is much larger than the image we are looking at, so we filter the table down to the sky area the Roman exposure actually covers. (Note: you can remove filters through the table UI if you wish to see the entire data)

```{code-cell} ipython3
# Bound the catalogs by the footprint of the Roman image we sent to Firefly.
ra_bounds = (roman_footprint[:, 0].min(), roman_footprint[:, 0].max())
dec_bounds = (roman_footprint[:, 1].min(), roman_footprint[:, 1].max())
```

You can visualize catalogs interactively with Firefly using [`show_table`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_table). This capability can take many parameters. Here we will simply send our catalog to Firefly so that we can (a) see an interactive table; (b) see this table plotted over the image that we've already sent; and (c) use the GUI to quickly create exploratory plots. See if you can use the GUI to quickly determine approximately how many galaxies cover the Roman exposure and what the redshift distribution of these galaxies is.

```{code-cell} ipython3
cat_filters = [
    f'("ra" >= {ra_bounds[0]} AND "ra" <= {ra_bounds[1]})', 
    f'("dec" >= {dec_bounds[0]} AND "dec" <= {dec_bounds[1]})'
]
cat_filters
```

```{code-cell} ipython3
fc.show_table(file_input=https_url(pointsource_cat_path),
              title='Stars Catalog',
              tbl_id='stars_cat',
              filters=" AND ".join(cat_filters))
```

```{code-cell} ipython3
gal_cat_tbl_id = 'galaxy_cat'

# may take ~1.25min, because galaxy catalog is a big file
fc.show_table(file_input=https_url(galaxy_cat_path),
              title='Galaxy Catalog',
              tbl_id=gal_cat_tbl_id,
              filters=" AND ".join(cat_filters))
```

For each row in the table you can notice a marker in the image. Selecting a row or marker changes the corresponding marker or row, respectively. You can click on "Details" tab in the UI to show properties of each source selected in image/table.

### Use Firefly's apply_table_filters to show sources of interest

Let's say the apparent physical size of a galaxy's disk component is what interests us. So we filter the galaxy catalog table we sent to Firefly down to the galaxies with a disk half-light radius bigger than 1 arcsec. You can use whatever properties you are interested in to filter the table.
Notice how the table display, the image overlay, and the active chart all change as the sample shrinks.
You can remove this filter or add new ones through the GUI.

For filtering, we will use [`apply_table_filters`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.apply_table_filters) method on the galaxy table we loaded above.

```{code-cell} ipython3
fc.apply_table_filters(tbl_id=gal_cat_tbl_id,
                       filters=" AND ".join(cat_filters+['"diskHalfLightRadiusArcsec" > 1']))
```

You can play with the filters directly from the UI as well. Try removing adding more filters in the tables and see how markers change.

+++

## 7. Identify a region of interest to compare cutouts

Sections 5 and 6 gave us two complementary ways to explore the data interactively: panning and zooming around the displayed images, and filtering the catalog table by whatever properties catch our interest. In this section we put both to use to land on a single region worth a closer look, copy its coordinates into the notebook, mark it with a region overlay, and then see how Roman and Rubin compare in that region.

+++

### Use Firefly's pan/zoom/align capabilities to locate a region of interest

You can view the coordinates of your mouse pointer at the bottom left of the display window. To copy the coordinates for a specific coordinate:

- Toggle the "Click Lock" to "on" in the bottom right of the image display.
- Click on the position of interest, and notice that the coordinate display is now frozen.
- Click on "EQ-J2000", the coordinate label in the bottom left of the image display. The pop-up that opens lets you choose the coordinate system, and the format used when copying — pick the Python format so the value can be pasted straight into a `SkyCoord`.
- Close the pop-up and click the small clipboard icon next to the coordinate readout.

+++

### Copy the coordinates from the coordinate display to the Python notebook

Paste the position you copied from the display here.
We picked a large disk-shaped galaxy that the Roman image resolves clearly, but that looks like just a fuzzy blob in the Rubin image. Its disk half-light radius of ~2 arcsec is well above the threshold we filtered on in Section 6, which is what caught our eye while panning around.

```{code-cell} ipython3
coords_of_interest = SkyCoord('0h38m36.18s -44d00m12.3s', frame='icrs') # pasted from the UI
coords_of_interest
```

Let's also define the size of the region we are interested in, which will be used to create cutouts from the images:

```{code-cell} ipython3
cutout_size = 48*u.arcsec
```

### Use ds9 region syntax to overplot the selected position

For this we use the id of the region layer we defined above, and add more region data using [`add_region_data`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.add_region_data).

```{code-cell} ipython3
point_region = (f'icrs;circle {coords_of_interest.ra.value}d {coords_of_interest.dec.value}d {cutout_size.to(u.deg).value} '
                f'# text={{Region of interest}} color=cyan width=2')
fc.add_region_data(region_data=point_region, region_layer_id=roman_regions_id,
                   plot_id=image_ff_id_rubin)
```

### Plot static side-by-side comparisons of the cutouts using matplotlib

Now that interactive exploration has narrowed things down to a specific region of interest, rather than the arbitrary region we started with in Section 2, let's compare Roman and Rubin cutouts of it the same way we did in Section 4.

```{code-cell} ipython3
# As in section 4, cut 1.5x oversized so the rotation has data for the corners,
# then resample down onto the cutout_size grid that gets displayed.
cutout_roman = Cutout2D(image_roman['data'], coords_of_interest, size=1.5*cutout_size, wcs=image_roman['wcs'])
cutout_rubin = Cutout2D(image_rubin['data'], coords_of_interest, size=1.5*cutout_size, wcs=image_rubin['wcs'])

roman_data, roman_wcs = north_up(cutout_roman.data, cutout_roman.wcs, cutout_size)
rubin_data, rubin_wcs = north_up(cutout_rubin.data, cutout_rubin.wcs, cutout_size)
```

```{code-cell} ipython3
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

axs[0].imshow(roman_data, origin='lower', 
              clim=stretch_color(roman_data, .5)
              )
axs[0].set_title(f"ROMAN in filter {filter_roman}")

axs[1].imshow(rubin_data, origin='lower', 
              clim=stretch_color(rubin_data, .5)
              )
axs[1].set_title(f"RUBIN in filter {filter_rubin}")

fig.suptitle(f"Cutouts at ({coords_of_interest.ra:6f}, {coords_of_interest.dec:6f}) with {cutout_size} size", fontsize=14);
plt.tight_layout(rect=[0, 0, 1, 0.97])
```

## 8. Use Firefly to visualize your region in 3-color
Firefly can also combine multiple bands into a false-color image. Since a single Rubin visit doesn't detect faint sources at all, we make this 3-color image from the Roman exposures instead. We use the following bands:

```{code-cell} ipython3
ROMAN_RGB_FILTERS = ['H158', 'J129', 'Y106']
```

### Retrieve closely-overlapping Roman images in the 3 bands
`image_roman` is the H158 exposure covering `coords_of_interest` from Section 7. Roman's WAS filters aren't tiled identically, so rather than searching each RGB filter around `coords_of_interest` again, we search around the center of this exposure itself, which keeps the 3 filters closely aligned.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
tags: [hide-cell]
---
def get_roman_fpaths_at_center(image, filters):
    """Find the Roman WAS exposures in other filters closest to an image's center.

    Parameters
    ----------
    image : dict
        Image dictionary as returned by `get_roman_image`, whose center
        defines the search position.
    filters : list of str
        Roman bandpass names to look up at that position, e.g. ``['H158', 'J129', 'Y106']``.

    Returns
    -------
    dict
        Filter name to S3 path of the closest exposure, or ``None`` for a
        filter with no WAS exposure covering the center.
    """
    ny, nx = image['data'].shape
    center = image['wcs'].pixel_to_world((nx - 1) / 2.0, (ny - 1) / 2.0)

    # One search covers every band, so ask once rather than searching again for each filter
    results = Irsa.query_sia(pos=(center, SEARCH_RADIUS),
                             collection=OU_ROMAN_SIA_COLLECTION)
    in_was = results[['WAS_simple_model' in str(r['obs_id']) for r in results]]

    # Sort by distance to the center of the image
    image_centers = SkyCoord(in_was['s_ra'], in_was['s_dec'], unit='deg')
    in_was['dist_to_center'] = center.separation(image_centers).deg
    in_was.sort('dist_to_center')

    # Sorted by distance, so the first row seen for a filter is the closest exposure for that filter
    fpaths = {}
    for row in in_was:
        band = str(row['energy_bandpassname']).strip()
        if band in filters and band not in fpaths:
            cloud_info = json.loads(row['cloud_access'])['aws']
            fpaths[band] = f"{cloud_info['bucket_name']}/{cloud_info['key']}"

    return {filter: fpaths.get(filter) for filter in filters}
```

Let's print the S3 paths of the Roman exposures in the 3 bands that cover our region of interest.

```{code-cell} ipython3
roman_fpaths = get_roman_fpaths_at_center(image_roman, ROMAN_RGB_FILTERS)
roman_fpaths
```

Verify that the H158 exposure we retrieved in Section 2 is the same as the one we found here.

```{code-cell} ipython3
image_roman['fpath'] == roman_fpaths['H158']
```

### Use Firefly's show_fits_3color to create a 3 color image
For this we will use [`show_fits_3color`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_fits_3color) method:

```{code-cell} ipython3
three_color_params = [{
    'URL': https_url(roman_fpaths[filter_roman]),
    'Title': "Roman 3 color"
} for filter_roman in ROMAN_RGB_FILTERS]

image_ff_id_roman_3color = 'roman-3color'
fc.show_fits_3color(three_color_params=three_color_params,
                    plot_id=image_ff_id_roman_3color)
```

### Use Firefly's pan/zoom/align methods to highlight the region of interest
Now, let's pan & zoom to the region of interest. Also align & lock all images being displayed by WCS. For these operations we use [`set_pan`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.set_pan), [`set_zoom`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.set_zoom) and [`align_images`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.align_images).

```{code-cell} ipython3
fc.set_pan(plot_id=image_ff_id_roman_3color, x=coords_of_interest.ra.deg, y=coords_of_interest.dec.deg, coord='j2000')
fc.set_zoom(plot_id=image_ff_id_roman_3color, factor=0.5)
fc.align_images(lock_match=True)
```

### Use Firefly's set_stretch method to change the stretch of the image display via Python

The image has a lot of noise that obscures our sources of interest. You can use the Firefly GUI to change the stretch of the image display. We identify that linear stretch from -1 to 30 sigma highlights the colors of our sources better. You can also use the Firefly client's [`set_stretch`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.set_stretch) to do this via Python. This is helpful for reproducibility and for scaling up to many images.

```{code-cell} ipython3
fc.set_stretch(plot_id=image_ff_id_roman_3color, stype='sigma', algorithm='linear', 
               band='ALL', lower_value=-1, upper_value=30)
```

***

## About This Notebook

**Updated:** 2026-09-02

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
