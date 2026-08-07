---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.3
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
authors:
  - name: Jaladh Singhal
  - name: Vandana Desai
  - name: IRSA Team
---

# Using Firefly to Explore OpenUniverse2024 Simulated Roman and Rubin Images

+++

## Learning Goals

By the end of this tutorial, you will:

- Learn how to access cloud-hosted Roman and Rubin simulated images.

- Learn how to locate a Roman coadd block from a sky position, and how to find Rubin images covering that position with the IRSA Simple Image Access (SIA) service.

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

This tutorial works with the full simulation, which covers the entire survey footprint. More information about the dataset can be found at [IRSA's holding of this dataset](https://irsa.ipac.caltech.edu/data/theory/openuniverse2024/overview.html), and the [OpenUniverse2024 paper](https://arxiv.org/abs/2501.05632) describes how the full simulation differs from the smaller preview subset that is also available.

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
- reproject.reproject_interp to convert Roman coadds from STG to TAN projection
- io.BytesIO for writing a fits file to an in-memory stream

```{code-cell} ipython3
# Uncomment the next line to install dependencies if needed.
# !pip install numpy astropy matplotlib firefly_client reproject astroquery hpgeom
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
```

## 1. Learn where the OpenUniverse2024 data are hosted in the cloud

The OpenUniverse2024 data are hosted in the cloud via Amazon Web Services (AWS). To access these data, you need to create a client to read data from Amazon's Simple Storage Service (s3) buckets, and you need to know some information about those buckets. OpenUniverse2024 contains simulations of the Roman Wide-Area Survey (WAS) and the Roman Time Domain Survey (TDS). In this tutorial, we will focus on the WAS.

The two telescopes need two different approaches. The Roman coadds sit on a fixed grid of sky positions, so once we know the grid we can work out a file path ourselves. The Rubin images are individual visits scattered across the sky, so instead we ask IRSA's Simple Image Access (SIA) service which ones cover the position we care about. We set up both routes here.

```{code-cell} ipython3
BUCKET_NAME = "nasa-irsa-simulations"
OU_PREFIX = "openuniverse2024"
ROMAN_COADD_PATH = f"{OU_PREFIX}/roman/full/RomanWAS/images/coadd"
TRUTH_FILES_PATH = f"{OU_PREFIX}/roman/full/roman_rubin_cats_v1.1.2_faint"
```

```{code-cell} ipython3
# Point the astroquery IRSA client at the simulated-data services, which are
# separate from the ones serving IRSA's observed data.
Irsa.sia_url = "https://irsa.ipac.caltech.edu/simulated/SIA"
Irsa.tap_url = "https://irsa.ipac.caltech.edu/simulated/TAP"

OU_RUBIN_SIA_COLLECTION = "simulated_rubin_openuniverse2024"

# A small radius is all we need, since we only want images containing a given point.
SEARCH_RADIUS = 1 * u.arcsec
```

## 2. Roman Coadds

The Nancy Grace Roman Space Telescope will carry out a wide-area survey (WAS) in the near infrared. OpenUniverse2024 includes coadded mosaics of simulated WAS data, created with the IMCOM algorithm (Rowe et al. 2011). Bands include F184, H158, J129, K213, Y106. In this section, we define some functions that make it convenient to retrieve a given cloud-hosted simulated Roman coadd based on position and filter.

+++

### Describe the grid of Roman simulated "blocks"

The simulated Roman coadds are arranged in blocks, as described in Hirata et al. 2024. Rather than listing the position of every block, we can describe the entire grid with a handful of numbers, because all 1296 blocks share one projection centered on the survey and differ only in where that center falls within each block. These values are recorded in the header of every coadd file, and we check them against a real header once we have opened one.

```{code-cell} ipython3
SURVEY_CRVAL = (9.55, -44.1)         # deg, the projection center shared by every block
PIXEL_SCALE = 1.0850694444444e-05    # deg/pixel
BLOCK_NPIX = 2688                    # pixels along each side of a block
BLOCK_STEP = 2560                    # pixel offset between the centers of adjacent blocks
N_BLOCKS = 36                        # the survey is a 36 x 36 grid of blocks
REF_BLOCK = 18                       # the block whose CRPIX sits at the survey center
REF_CRPIX = 64.5                     # that block's CRPIX, in the 1-based FITS convention
```

Each block is wider than the spacing between blocks, so neighbors overlap slightly rather than butting up against each other.

```{code-cell} ipython3
block_size = (BLOCK_NPIX * PIXEL_SCALE * u.deg).to(u.arcsec)
block_spacing = (BLOCK_STEP * PIXEL_SCALE * u.deg).to(u.arcsec)

print(f"Each block is {block_size:.1f} across, laid down every {block_spacing:.1f}")
```

### Define a function that returns the WCS of any block

Because the blocks differ only in their CRPIX, we can write down the WCS of any block without opening its file.

```{code-cell} ipython3
def block_wcs(col, row):
    """
    Build the WCS of the Roman WAS coadd block in a given column and row.

    Parameters
    ----------
    col, row : int
        Block indices, each running from 0 to 35.

    Returns
    -------
    astropy.wcs.WCS
        The WCS of that block.
    """
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---STG", "DEC--STG"]
    w.wcs.crval = list(SURVEY_CRVAL)
    w.wcs.cdelt = [-PIXEL_SCALE, PIXEL_SCALE]
    w.wcs.crpix = [REF_CRPIX - (col - REF_BLOCK) * BLOCK_STEP,
                   REF_CRPIX - (row - REF_BLOCK) * BLOCK_STEP]
    w.array_shape = (BLOCK_NPIX, BLOCK_NPIX)
    return w
```

### Define a function that returns the block containing a given sky position

To find the block holding a position, we express that position in the pixel grid of the reference block, then count how many block widths it lands away from that block's center.

```{code-cell} ipython3
def find_block(coord):
    """
    Find the Roman WAS coadd block that contains a sky position.

    Parameters
    ----------
    coord : astropy.coordinates.SkyCoord
        The position of interest.

    Returns
    -------
    tuple of int
        The (column, row) of the block containing that position.
    """
    x, y = block_wcs(REF_BLOCK, REF_BLOCK).world_to_pixel(coord)
    block_center = (BLOCK_NPIX - 1) / 2

    col = REF_BLOCK + int(np.round((x - block_center) / BLOCK_STEP))
    row = REF_BLOCK + int(np.round((y - block_center) / BLOCK_STEP))

    if not (0 <= col < N_BLOCKS and 0 <= row < N_BLOCKS):
        raise ValueError(f"{coord.to_string('hmsdms')} is not covered by the "
                         "OpenUniverse2024 simulated Roman coadds")

    return col, row
```

### Define a function that retrieves a Roman simulated coadd given a sky position and filter.

+++

Each of the cloud-hosted simulated Roman coadds can be accessed via a S3 filepath. This function returns the access path for the simulated Roman coadd that includes a specified position on the sky and desired filter.

```{code-cell} ipython3
def get_roman_coadd_fpath(coord, filter):
    col, row = find_block(coord)

    # Construct the coadd filename from the chosen filter, column, and row.
    coadd_fname_root = f"prod_{filter[0]}_{col:02d}_{row:02d}_map.fits.gz"
    coadd_fpath = f"{BUCKET_NAME}/{ROMAN_COADD_PATH}/{filter}/{coadd_fname_root}"
    return coadd_fpath
```

Now we prefix that path with `s3://` and open it with astropy.

We use [`.section`](https://docs.astropy.org/en/stable/io/fits/usage/cloud.html#subsetting-fits-files-hosted-in-amazon-s3-cloud-storage) to pull out just the science plane as a 2D `numpy.array`, so the full 15-plane cube is never assembled in memory, and we take the WCS from the fits header as an `astropy.wcs.WCS` object. The function below returns a dictionary of both, along with the header itself so we can look at it later.

A note on speed: `.section` can also cut down how much of a cloud-hosted file travels across the network, but only when the file can be read out of order. These coadds are gzipped, and a compressed stream has to be read from the beginning, so astropy works through the whole file to reach the plane we want. Expect each of these reads to take twenty to thirty seconds.

```{code-cell} ipython3
def get_roman_coadd(coord, filter):
    # retrive fits file of block/tile from the coadd mosiac
    coadd_s3_fpath = get_roman_coadd_fpath(coord, filter)
    coadd_s3_uri = f"s3://{coadd_s3_fpath}"

    with fits.open(coadd_s3_uri, fsspec_kwargs={"anon": True}) as hdul:
        # retrieve science data from coadd fits
        coadd_data = hdul[0].section[0, 0, :, :]  # has (1, 15, 2688, 2688) shape, with 0th layer in the cube as science image

        # make wcs using header
        coadd_wcs = wcs.WCS(hdul[0].header, naxis=2)

        return {'data': coadd_data, 'wcs': coadd_wcs, 'header': hdul[0].header}
```

### Inspect a simulated Roman Coadd

+++

Choose a filter and a position within the survey footprint

```{code-cell} ipython3
coord = SkyCoord(ra=9.6055383, dec=-44.1895542, unit="deg")
filter_roman = 'H158' #F184, H158, J129, K213, and Y106 are available
```

```{code-cell} ipython3
print("This position falls in block (column, row):", find_block(coord))
```

Retrieve the data and header information from the simulated Roman coadd corresponding to the chosen position and filter.

```{code-cell} ipython3
coadd_roman = get_roman_coadd(coord, filter_roman)
```

With a real header in hand, we can confirm that the grid description above reproduces it exactly.

```{code-cell} ipython3
predicted = block_wcs(*find_block(coord))

print("CRPIX from the file:", coadd_roman['header']['CRPIX1'], coadd_roman['header']['CRPIX2'])
print("CRPIX from our grid:", *predicted.wcs.crpix)
```

### Understand the size of a simulated Roman coadd.

```{code-cell} ipython3
# Number of pixels (Y, X)
coadd_roman['data'].shape
```

```{code-cell} ipython3
# Pixel size (scale Y, scale X) [degrees/pixel]
coadd_roman['wcs'].proj_plane_pixel_scales()
```

```{code-cell} ipython3
# Coadd size (FOV Y, FOV X)
[(num * size).to('arcsec') for num, size in zip(
    coadd_roman['data'].shape, coadd_roman['wcs'].proj_plane_pixel_scales())]
```

The field of view of Roman coadd is ~100 arcsec.

+++

### Look at everything else the coadd file contains

So far we have read a single plane out of the primary HDU, but each coadd file carries a good deal more than the science image. Listing the extensions shows the full picture. Finding each extension means walking past the one before it, so this cell has to work through the entire compressed file and takes a few seconds.

```{code-cell} ipython3
with fits.open(f"s3://{get_roman_coadd_fpath(coord, filter_roman)}",
               fsspec_kwargs={"anon": True}) as hdul:
    hdul.info()
```

The Primary HDU for the coadded image is a cube with 15 layers, i.e., its shape is 1x15x2688x2688. The layers are as follows:

0 = simulated "Science" image (Roman+Rubin simulation, units of e/(0.11 arcsec)^2/exposure)

1 = lab noise: based on dark frames from the April 2023 test, masked at 3 e/p/s. Units: e/(0.11 arcsec)^2/s

2 = GalSim stars, on HEALPix resolution 14 grid, normalized to total flux of 1

3 = noisy stars, on HEALPix resolution 14 grid, normalized to total flux of 2.4e5 e with self-Poisson noise, including 86 e^2/input pixel background variance

4 = stars, on HEALPix resolution 14 grid, total flux 1, but on in only one of the passes (to test transient response)

5 = stars, on HEALPix resolution 14 grid, with total flux that varies by 5% from center to edge of the focal plane (to test what happens when the filter bandpass varies; 5% is highly exaggerated)

6 = GalSim extended objects, on HEALPix resolution 14 grid, right now exponential profiles. The scale radius is log-distributed between 0.125 and 0.500 arcsec, and the ellipticity (g1,g2) is uniformly distributed in the disc of radius 0.5, i.e., g1^2+g2^2<0.5^2.

7,8,9 = same objects as layer 6, but with applied shear of 0.02. The shear orientations are spaced by 60° in tangent vector space, so that in the (g1,g2)-space they are spaced by 120° and can be used for finite differences. Specifically, the directions are: layer 7 -> in East-West direction (shear PA = 270°). (g1,g2) = (0.02,0) layer 8 -> in NNW-SSE direction (shear PA = 330°). (g1,g2) = (-0.02/2,0.02√3/2) layer 9 -> in NNE-SSW direction (shear PA = 30°). (g1,g2) = (-0.02/2,-0.02√3/2)

10 = coadded 1/f noise map, normalized to variance per ln f of 1

11,12,13,14 = coadded white noise maps, different seeds, normalized to variance of 1 in each input pixel

The remaining HDUs contain additional information:

CONFIG = the configuration file

INDATA = the input images used, as a binary table. The columns are: obsid (int32) -> observation ID sca (int16) -> SCA (1 through 18, inclusive) ra (float64) -> right ascension of pointing center in degrees dec (float64) -> declination of pointing center in degrees pa (float64) -> position angle of pointing in degrees valid (logical) -> input science data file is valid (should be True)

INWEIGHT = the mean input weights for how much each 1.25x1.25 arcsec postage stamp depends on each input exposure. The shape is 1 x Nin x 84 x 84, where Nin is the number of input images listed in INDATA. Note that each postage stamp is 32 output pixels, so 84x32=2688.

If summed on axis 1, this will normally be something close to 1. Deviations of ~10% are common, due to plate scale variations and the normalization issues introduced by diffraction spikes.

INWTFLAT = a reshape of INWEIGHT suitable for display as a single image in a viewer such as DS9. The contributions from the Nin input exposures are rearranged into a 1 x 84 x (N_in*84) array.

FIDELITY, SIGMA, INWTSUM, EFFCOVER = maps of U_alpha/C, S_alpha, sum_i T_{alpha i}, and the effective coverage as rescaled int16 maps. See Rowe et al. (2011) for details on the definitions of these quantities. The comment in the 'UNIT' keyword indicates how to rescale these.

+++

### Use the WCS from the Roman simulated header to convert the specified coordinate into a pixel position.

```{code-cell} ipython3
def coord_to_xy(w, coord):
    return w.world_to_array_index(coord)[::-1] #reverse since 0th axis is y, 1st axis is x

coord_arr_idx = coord_to_xy(coadd_roman['wcs'], coord)
coord_arr_idx
```

### Use matplotlib imshow to create a static visualization of the Roman simulated coadd and overplot the selected position.

```{code-cell} ipython3
def stretch_color(data, clipPercent):
    return np.percentile(data, (0 + clipPercent, 100 - clipPercent))

plt.imshow(coadd_roman['data'], origin='lower', 
           clim=stretch_color(coadd_roman['data'], 1)
           )

plt.plot(*coord_arr_idx, 'r+', markersize=15)
```

## 3. Rubin Images

OpenUniverse2024 includes simulated Rubin images in the following filters: u, g, r, i, z, y. These are individual visits rather than a fixed grid of mosaics, so instead of building a path we let the SIA service tell us which images cover our position. In this section, we define functions that retrieve a Rubin image for a chosen position and filter, returning data in the same structure as the functions we defined above for Roman.

+++

### Retrieve Rubin images

The image search returns one row per image, with the cloud location of each tucked inside a `cloud_access` JSON string. Many visits cover any given position, so we sort by observation time and take the earliest; that keeps the notebook reproducible from one run to the next.

```{code-cell} ipython3
def get_rubin_image_fpath(coord, filter):
    results = Irsa.query_sia(pos=(coord, SEARCH_RADIUS),
                             collection=OU_RUBIN_SIA_COLLECTION)

    # Band names come back like "r_57", so compare only the part before the underscore.
    in_band = results[[str(r['energy_bandpassname']).split('_')[0] == filter
                       for r in results]]
    if len(in_band) == 0:
        raise ValueError(f"No Rubin {filter}-band images cover "
                         f"{coord.to_string('hmsdms')}")

    in_band.sort('t_min')
    cloud_info = json.loads(in_band['cloud_access'][0])['aws']
    return f"{cloud_info['bucket_name']}/{cloud_info['key']}"
```

```{code-cell} ipython3
def get_rubin_image(coord, filter):
    image_s3_fpath = get_rubin_image_fpath(coord, filter)

    with fits.open(f"s3://{image_s3_fpath}", fsspec_kwargs={"anon": True}) as hdul:
        # retrieve science data, which sits in the first extension
        image_data = hdul[1].data

        # make wcs using header
        image_wcs = wcs.WCS(hdul[1].header)

        return {'data': image_data, 'wcs': image_wcs}
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

A single Rubin detector covers a far larger patch of sky than one Roman coadd block, at a coarser pixel scale.

+++

### Use matplotlib imshow to create a static visualization of the Rubin simulated image and overplot the selected position.

```{code-cell} ipython3
plt.imshow(image_rubin['data'], origin='lower', 
           clim=stretch_color(image_rubin['data'], 1)
           )

plt.plot(*coord_to_xy(image_rubin['wcs'], coord), 'r+', markersize=15)
```

### Define a function that returns the URL for a given S3 filepath
Since the OpenUniverse2024 data is available through a public S3 bucket, we can access a given S3 file using HTTPS URL as follows:

```{code-cell} ipython3
def https_url(s3_fpath):
    s3_fpath_without_bucket = s3_fpath.split('/', 1)[1]
    return f"https://{BUCKET_NAME}.s3.amazonaws.com/{s3_fpath_without_bucket}"
```

Let's generate URL for the Rubin image we plotted above. Clicking on the returned URL will allow you to download this image locally.

```{code-cell} ipython3
image_s3_fpath_rubin = get_rubin_image_fpath(coord, filter_rubin)
https_url(image_s3_fpath_rubin)
```

## 4. Compare simulated Roman and Rubin cutouts for a selected position

+++

### Choose cutout size

```{code-cell} ipython3
cutout_size = 50*u.arcsec
```

### Create the cutouts

```{code-cell} ipython3
cutout_roman = Cutout2D(coadd_roman['data'], coord, size=cutout_size, wcs=coadd_roman['wcs'])
cutout_rubin = Cutout2D(image_rubin['data'], coord, size=cutout_size, wcs=image_rubin['wcs'])
```

### Use matplotlib imshow to plot static side-by-side comparisons of the cutouts

```{code-cell} ipython3
fig, axs = plt.subplots(1, 2, figsize=(12, 6))


axs[0].imshow(cutout_roman.data, origin='lower', 
              clim=stretch_color(cutout_roman.data, .5)
              )
axs[0].plot(*coord_to_xy(cutout_roman.wcs, coord), 'r+', markersize=15)
axs[0].set_title(f"ROMAN in filter {filter_roman}")

axs[1].imshow(cutout_rubin.data, origin='lower', 
              clim=stretch_color(cutout_rubin.data, .5)
              )
axs[1].plot(*coord_to_xy(cutout_rubin.wcs, coord), 'r+', markersize=15)
axs[1].set_title(f"RUBIN in filter {filter_rubin}")

fig.suptitle(f"Cutouts at ({coord.ra}, {coord.dec}) with {cutout_size} size", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.97])
```

## 5. Use Firefly to interactively identify a blended source

Clearly, the simulated Roman coadd has higher spatial resolution than the simulated Rubin image. Let's try to locate blended objects to compare in the simulated Rubin and Roman images. We will use Firefly's interactive visualization to make this task easier.

+++

### Launch and initialize Firefly
There are two ways to initialize a Firefly client from Python, depending on whether you're running the notebook in JupyterLab or not. Assuming you have `jupyter-firefly-extensions` set up in your environment as explained [here](https://github.com/Caltech-IPAC/jupyter_firefly_extensions/blob/master/README.md), you can use `make_lab_client()` in JupyterLab, which will open the Firefly viewer in a new tab within the Lab. Otherwise, you can use `make_client()` in a Jupyter Notebook (or even a Python shell), which will open the Firefly viewer in a new web browser tab.

You also need a Firefly server to communicate with your Firefly Python client. In this notebook, we use a public Firefly server: the IRSA Viewer (https://irsa.ipac.caltech.edu/irsaviewer). However, you can also run a local Firefly server via a [Firefly Docker image](https://hub.docker.com/r/ipac/firefly) and access it at `http://localhost:8080/firefly`. The URL of the Firefly server is read by both `make_client()` and `make_lab_client()` through the environment variable `FIREFLY_URL`. However, `make_client()` also allows you to pass the URL directly as the `url` parameter.

```{code-cell} ipython3
# Uncomment when using within Jupyter Lab with jupyter_firefly_extensions installed
# fc = FireflyClient.make_lab_client()

# Uncomment for contexts other than above 
fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")

fc.reinit_viewer() # to clean the state, if this cell ran earlier
```

### Send the simulated Rubin image to Firefly using show_fits.

For displaying the FITS image of the Rubin visit in Firefly, we use [`show_fits`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_fits):

```{code-cell} ipython3
image_ff_id_rubin = 'rubin-image-filter-r'
fc.show_fits(url=https_url(image_s3_fpath_rubin),
             plot_id=image_ff_id_rubin,
             Title="Rubin Image"
             )
```

### Use ds9 region syntax to overplot the simulated Roman image blocks on the interactive display

The Firefly client includes several methods related to controlling ds9 region overlays. To 
overlay a region layer on the loaded FITS images, we can use [`overlay_region_layer`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.overlay_region_layer).

Region data is defined in ds9 region syntax that can be found [here](https://ds9.si.edu/doc/ref/region.html).

Drawing all 1296 blocks would clutter the display, so we outline only the ones in the neighborhood of our target. Their centers come straight from the WCS of each block.

```{code-cell} ipython3
def block_center(col, row):
    """Sky position of the center of a given coadd block."""
    return block_wcs(col, row).pixel_to_world((BLOCK_NPIX - 1) / 2,
                                              (BLOCK_NPIX - 1) / 2)

col, row = find_block(coord)
nearby_blocks = [(c, r)
                 for c in range(max(col - 3, 0), min(col + 4, N_BLOCKS))
                 for r in range(max(row - 3, 0), min(row + 4, N_BLOCKS))]
```

```{code-cell} ipython3
# mark the roman coadd blocks as boxes
roman_regions = [
    f'icrs;box {center.ra.deg}d {center.dec.deg}d {block_size.value}" {block_size.value}" 0d'
    for center in (block_center(c, r) for c, r in nearby_blocks)
]

roman_regions_id = 'roman_regions'
fc.overlay_region_layer(region_data=roman_regions,
                        title='Roman Mosiac', 
                        region_layer_id=roman_regions_id)
```

### Use Firefly's pan and zoom capabilities to locate a region of interest (a blended source)

You can view the coordinates of your mouse pointer at the bottom left of the display window. To copy the coordinates for a specific coordinate:

- Toggle the "Click Lock" to "on" in the bottom right of the image display.
- Click on the position of interest, and notice that the coordinate display is now frozen.
- Click on "EQ-J2000", the coordinate label in the bottom left of the image display. In the dialog that opens, change copy options to "[Python] Astropy SkyCoord" so that we can directly work with them in python.
- Close the dialog and click on the copy icon next to the coordinate values display.

+++

### Copy the coordinates from the coordinate display to the Python notebook

We have provided an example. You can change this based on your interests.

```{code-cell} ipython3
coords_of_interest = SkyCoord('0h38m25.35s -44d00m10.1s', frame='icrs') # located and copied through UI
coords_of_interest
```

We can now use this [astropy `SkyCoord` object](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord) to compare our images.

+++

### Use ds9 region syntax to overplot the selected position

For this we use the id of the region layer we defined above, and add more region data using [`add_region_data`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.add_region_data).

```{code-cell} ipython3
point_region = f'icrs;point {coords_of_interest.ra.value}d {coords_of_interest.dec.value}d # point=cross 15 text={{Blended source}}'
fc.add_region_data(region_data=point_region, region_layer_id=roman_regions_id)
```

## 6. Plot cutouts of the identified blended source

```{code-cell} ipython3
coadd_roman = get_roman_coadd(coords_of_interest, filter_roman)
image_rubin = get_rubin_image(coords_of_interest, filter_rubin)
```

```{code-cell} ipython3
cutout_size = 20*u.arcsec
```

```{code-cell} ipython3
cutout_roman = Cutout2D(coadd_roman['data'], coords_of_interest, size=cutout_size, wcs=coadd_roman['wcs'])
cutout_rubin = Cutout2D(image_rubin['data'], coords_of_interest, size=cutout_size, wcs=image_rubin['wcs'])
```

```{code-cell} ipython3
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

axs[0].imshow(cutout_roman.data, origin='lower', 
              clim=stretch_color(cutout_roman.data, .5)
              )
axs[0].set_title(f"ROMAN in filter {filter_roman}")

# Let's also encircle the blended source we identified 
axs[0].add_patch(patches.Circle(coord_to_xy(cutout_roman.wcs, coords_of_interest), 
                                radius=50, color='r', fill=False, linewidth=2))
# and a bonus blended source that is close to it
other_coords = SkyCoord(coords_of_interest.ra-9.2*u.arcsec, coords_of_interest.dec+5*u.arcsec)
axs[0].add_patch(patches.Circle(coord_to_xy(cutout_roman.wcs, other_coords),
                                radius=36, color='cyan', fill=False, linewidth=2))


axs[1].imshow(cutout_rubin.data, origin='lower', 
              clim=stretch_color(cutout_rubin.data, .5)
              )
axs[1].set_title(f"RUBIN in filter {filter_rubin}")

# Let's also encircle the source we identified 
axs[1].add_patch(patches.Circle(coord_to_xy(cutout_rubin.wcs, coords_of_interest),
                                radius=10, color='r', fill=False, linewidth=2))
axs[1].add_patch(patches.Circle(coord_to_xy(cutout_rubin.wcs, other_coords),
                                radius=8, color='cyan', fill=False, linewidth=2))


fig.suptitle(f"Cutouts at ({coords_of_interest.ra:6f}, {coords_of_interest.dec:6f}) with {cutout_size} size", fontsize=14);
plt.tight_layout(rect=[0, 0, 1, 0.97])
# plt.savefig("plot.pdf", bbox_inches='tight', pad_inches=0.2)
```

## 7. Use Firefly to visualize the OpenUniverse2024 catalogs
Let's inspect the properties of sources in the Rubin image. For this we will use the input truth files present in S3 bucket.

OpenUniverse2024 includes the input truth files that were used to create the simulated images. These files are in Parquet and HDF5 format, and include information about the properties of galaxies, stars, and transients.

The catalogs are split by HEALPix sky region (nside=32, RING ordering) and the region index appears in each filename, so we can convert our sky position into that index with [`hpgeom`](https://hpgeom.readthedocs.io/en/latest/) rather than hunting through the directory listing. This is the same approach the [Quickstart](openuniverse2024_quickstart) tutorial takes.

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

Each catalog spans a whole HEALPix region, which is much larger than the image we are looking at, so we filter the table down to the sky area the Rubin image actually covers. (Note: you can remove filters through the table UI if you wish to see the entire data)

```{code-cell} ipython3
# Bound the catalogs by the footprint of the Rubin image we sent to Firefly.
footprint = image_rubin['wcs'].calc_footprint()
ra_bounds = (footprint[:, 0].min(), footprint[:, 0].max())
dec_bounds = (footprint[:, 1].min(), footprint[:, 1].max())
```

+++

You can visualize catalogs interactively with Firefly using [`show_table`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_table). This capability can take many parameters. Here we will simply send our catalog to Firefly so that we can (a) see an interactive table; (b) see this table plotted over the image that we've already sent; and (c) use the GUI to quickly create exploratory plots. See if you can use the GUI to quickly determine approximately how many galaxies cover the Rubin image and what the redshift distribution of these galaxies is.

```{code-cell} ipython3
cat_filters = [
    f'("ra" >= {ra_bounds[0]} AND "ra" <= {ra_bounds[1]})', 
    f'("dec" >= {dec_bounds[0]} AND "dec" <= {dec_bounds[1]})'
]
cat_filters
```

```{code-cell} ipython3
fc.show_table(url=https_url(pointsource_cat_path),
              title='Stars Catalog',
              tbl_id='stars_cat',
              filters=" AND ".join(cat_filters))
```

```{code-cell} ipython3
gal_cat_tbl_id = 'galaxy_cat'

# may take ~1.25min, because galaxy catalog is a big file
fc.show_table(url=https_url(galaxy_cat_path),
              title='Galaxy Catalog',
              tbl_id=gal_cat_tbl_id,
              filters=" AND ".join(cat_filters))
```

For each row in the table you can notice a marker in the image. Selecting a row or marker changes the corresponding marker or row, respectively. You can click on "Details" tab in the UI to show properties of each source selected in image/table.

### Use Firefly's apply_table_filters to show only high-redshift galaxies

High redshift galaxies are the most interesting, so let's filter the table we sent to Firefly to only include z>3 galaxies. Notice how the table display and image overlay change. Notice how the chart becomes a scatterplot from a heatmap because the sources reduce. You can remove this filter or add new ones through the GUI.

For filtering, we will use [`apply_table_filters`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.apply_table_filters) method on the galaxy table we loaded above.

```{code-cell} ipython3
fc.apply_table_filters(tbl_id=gal_cat_tbl_id,
                       filters=" AND ".join(cat_filters+['"redshift" > 3']))
```

You can play with the filters directly from the UI as well. Try removing adding more filters in the tables and see how markers change.

+++

### Use Firefly's show_fits_3color to create a 3 color image of the simulated Rubin images

```{code-cell} ipython3
# [R, G, B]
ROMAN_RGB_FILTERS = ['H158', 'J129', 'Y106']
RUBIN_RGB_FILTERS = ['r', 'g', 'u']
```

We already have a Rubin image with the catalog overlaid, let's make a 3 color image to see colors of marked objects more clearly. For this we will use [`show_fits_3color`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_fits_3color) method:

```{code-cell} ipython3
image_ff_id_rubin_3color = 'rubin-image-3color'
threeC = [
    dict(url=https_url(get_rubin_image_fpath(coord, filter_name)),
         Title="Rubin 3 color")
    for filter_name in RUBIN_RGB_FILTERS
]

fc.show_fits_3color(three_color_params=threeC,
                    plot_id=image_ff_id_rubin_3color)
```

### Use Firefly's interactivity to identify a region of interest 

For example, we found a region of the sky that seems to have a high number of high redshift sources and we copy it from the image display:

```{code-cell} ipython3
# located and copied through UI
high_z_gal_coords = SkyCoord('0h38m00.77s -44d12m10.2s', frame='icrs')
high_z_gal_coords
```

```{code-cell} ipython3
# let's also mark it in our region layer, so that it's easy to pinpoint later
point_region = f'icrs;point {high_z_gal_coords.ra.value}d {high_z_gal_coords.dec.value}d # point=cross 15 text={{z>3 mock galaxies}}'
fc.add_region_data(region_data=point_region, region_layer_id=roman_regions_id)
```

## 8. Plot 3-color Roman coadd containing your region of interest
Let's inspect WCS of Roman coadd first

```{code-cell} ipython3
coadd_roman['wcs']
```

### Prepare Roman coadds for displaying in Firefly

Roman coadds have STG projection which cannot be read by Firefly yet. Firefly can display any FITS image but it needs to read the WCS for overlaying catalogs and other interactive features. So unlike the Rubin 3 color image where we directly passed URLs to Firefly, we will read Roman coadd files in Python, reproject them from STG to TAN, and write them back to FITS to pass them to Firefly.

Let's first define functions to do so:

```{code-cell} ipython3
def reproject_to_TAN(coadd_roman):
    # Define a new WCS with TAN projection (in CTYPE key)
    output_wcs = coadd_roman['wcs'].deepcopy()
    output_wcs.wcs.ctype = [ctype.replace('STG', 'TAN') for ctype in coadd_roman['wcs'].wcs.ctype]
    
    # Use reproject to convert a given data and wcs, to a desired wcs and shape
    reprojected_data, _ = reproject_interp(
        (coadd_roman['data'], coadd_roman['wcs']),
        output_projection=output_wcs,
        shape_out=coadd_roman['data'].shape
    )

    return {'data': reprojected_data, 'wcs': output_wcs}
```

```{code-cell} ipython3
def get_fits_stream(coadd_roman):
    # Create a FITS PrimaryHDU object with the coadd data
    hdu = fits.PrimaryHDU(data=coadd_roman['data'], header=coadd_roman['wcs'].to_header())

    # Write the HDU to the in-memory stream (to save I/O time)
    fits_stream = BytesIO()
    hdu.writeto(fits_stream, overwrite=True)
    fits_stream.seek(0) # to bring reading pointer to the beginning of file

    return fits_stream
```

Then we perform all 3 operations we mentioned above for the RGB filters of Roman. This reads another coadd for each band, so the cell takes a couple of minutes:

```{code-cell} ipython3
coadds_rgb = []
coadds_rgb_reprojected = []
coadds_rgb_fits_stream = []

for filter_name in ROMAN_RGB_FILTERS:
    print(f'\nFILTER: {filter_name}')
    print('Retrieving Roman coadd...')
    coadd_roman = get_roman_coadd(high_z_gal_coords, filter_name)
    coadds_rgb.append(coadd_roman)

    print('Reprojecting to TAN...')
    coadd_roman_reprojected = reproject_to_TAN(coadd_roman)
    coadds_rgb_reprojected.append(coadd_roman_reprojected)

    print('Writing back to fits stream...')
    coadd_roman_fits_stream = get_fits_stream(coadd_roman_reprojected)
    coadds_rgb_fits_stream.append(coadd_roman_fits_stream)
```

### Use Firefly's show_fits_3color to create a 3 color image of Roman coadds

Now we upload each fits stream (in-memory fits file) to firefly using [`upload_fits_data()`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.upload_fits_data) and prepare color params to pass to the [`show_fits_3color()`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.show_fits_3color).

```{code-cell} ipython3
three_color_params = [
    {
        'file': fc.upload_fits_data(fits_stream),
        'Title': "Roman Coadd 3 color"
    } for fits_stream in coadds_rgb_fits_stream]
```

```{code-cell} ipython3
coadd_ff_id_roman_3color = 'roman-coadd-3color-high_z_gal'
fc.show_fits_3color(three_color_params=three_color_params,
                    plot_id=coadd_ff_id_roman_3color)
```

We can see 3 color image of Roman coadd containing the high-redshift galaxy sources. Try panning and zomming out, you can notice it spans over one block compared to the Rubin image which is much larger.

### Use Firefly's pan/zoom/align methods to locate high redshift sources
Now, let's pan & zoom to the region where we located high-redshift galaxy sources. Also align & lock all images being displayed by WCS. For these operations we use these 3 methods (click on them to see their documentation):
- [`set_pan`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.set_pan)
- [`set_zoom`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.set_zoom)
- [`align_images`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.align_images)

```{code-cell} ipython3
fc.set_pan(plot_id=coadd_ff_id_roman_3color, x=high_z_gal_coords.ra.deg, y=high_z_gal_coords.dec.deg, coord='j2000')
fc.set_zoom(plot_id=coadd_ff_id_roman_3color, factor=1)
fc.align_images(lock_match=True)
```

### Use Firefly's set_stretch method to change the stretch of the image display via Python

The image has a lot of noise that obscures our high redshift sources of interest. You can use the Firefly GUI to change the stretch of the image display. We identify that squared stretch from -2 to 10 sigma highlights the colors of our sources better. You can also use the Firefly client's [`set_stretch`](https://caltech-ipac.github.io/firefly_client/api/firefly_client.FireflyClient.html#firefly_client.FireflyClient.set_stretch) to do this via Python. This is helpful for reproducibility and for scaling up to many images.

```{code-cell} ipython3
fc.set_stretch(plot_id=coadd_ff_id_roman_3color, stype='sigma', algorithm='squared', 
               band='ALL', lower_value=-2, upper_value=10)
```

***

## About This Notebook

**Updated:** 2026-08-05

**Contact:** [the IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions or reporting problems.
