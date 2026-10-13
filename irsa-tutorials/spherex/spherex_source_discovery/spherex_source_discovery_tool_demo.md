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
  - name: Gabriela Torrini
  - name: Andreas Faisst
  - name: SPHEREx Science Data Center (SSDC)
  - name: Troy Raen
  - name: Brigitta Sipőcz
  - name: Jaladh Singhal
---

# SPHEREx Source Discovery Tool IRSA Demo

The SPHEREx Source Discovery Tool is the Python package `spherex_source_discovery_tool` (included in this directory), which is used to discover and extract sources from SPHEREx Spectral Images and visualize their spectra.
This notebook demonstrates how to use it.

## 1. Learning Goals

* Search and download SPHEREx images

* Implement a simple image subtraction to find red sources in the SPHEREx images

* Visualize images and catalog data in _Firefly_ as well as interactive _bokeh_ plots.

* Select sources interactively and extract their SPHEREx spectra

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

## 3. Requirements

+++

### 3.1 On local machines

All necessary packages and tools (e.g., <code>SExtractor</code>) are listed in the `conda-spherex_sdt.yml` file.

To create a conda environment with the dependencies on your local machine use:
  ```
  conda env create --file conda-spherex_sdt.yml
  ```

Then, make the environment available in the list of kernels for your JupyterLab:
  ```
  conda activate spherex_sdt
  python -m ipykernel install --user --name=spherex_sdt
  ```

To use the environment in your Jupyter notebooks, either start JupyterLab in that environment by typing

```
conda activate spherex_sdt
jupyter-lab
```

or select the environment `spherex_sdt` in your Jupyter Notebook using the dropdown on the upper left.

### 3.2 On Fornax

Installing the conda environment on the [NASA Fornax Science Console](https://science.nasa.gov/astrophysics/programs/physics-of-the-cosmos/community/the-fornax-initiative/) needs slightly different steps. These can be reviewed in the documentation [create a new environment](https://docs.fornax.sciencecloud.nasa.gov/compute-environments/#create-new-env).

In order to install this specific conda environment on Fornax, the file name of the `yml` file specifically needs to be in the format `conda-*.yml`. The `yml` file distributed here is already provided in that format (`conda-spherex_sdt.yml`). Once this is set, open a new terminal (click on the large "+" button right under the "File" menu tab) and type to following command _inside_ the same directory where the `yml` file is located:

```
setup-conda-env --user
```

Note that we use the `--user` option here, which will keep the environment available for subsequent Fornax sessions.

To use the environment in your Jupyter notebooks on Fornax, directly select the environment `spherex_sdt` in your Jupyter Notebook using the dropdown on the upper left.

+++

## 4. Imports

```{code-cell} ipython3
# Standard library imports
import concurrent.futures
import copy
import json
import os
import io
import pickle
import time
from concurrent.futures import as_completed

# Related 3rd-party imports
import numpy as np
import pyvo
import sep
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.units import Unit
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa
from bokeh.io import output_notebook
from bokeh.resources import INLINE
from firefly_client import FireflyClient
from reproject import reproject_exact, reproject_interp

# Local application imports
from spherex_source_discovery_tool.aperture_photometry import build_phot_table, grab_star
from spherex_source_discovery_tool.bokeh_viz import (
    plot_apertures, plot_spectrum, plot_subtracted_trio, plot_overlap_trio)
from spherex_source_discovery_tool.firefly_viz import preview_query
from spherex_source_discovery_tool.sdt_utils import (
    format_extracted, get_filename, get_exp_id, get_obs_id, results_summary, get_lambda_range)
from spherex_source_discovery_tool.source_extraction import run_sextractor, get_sextractor_file
```

```{code-cell} ipython3
# Suppress logging temporarily to prevent astropy
# from repeatedly printing out warning notices related to alternate WCSs
import logging
logging.getLogger('astropy').setLevel(logging.ERROR)

# The time it takes to read SPHEREx files can exceed
# astropy's default timeout limit. Increase it.
from astropy.utils.data import conf
conf.remote_timeout = 120

# For interactive bokeh visualizations
output_notebook(resources=INLINE)
```

## 5. Initialize Ephemeral Query Parameters

```{attention}
The query parameters below are valid as of **February 2026**, but will need to be updated as the mission proceeds.
Specifically, the ephemeral query parameters _will_ change for future SPHEREx data releases.
```

```{code-cell} ipython3
# General
data_release = "qr2"

# Spectral image queries
sia_collection = f"spherex_{data_release}"

# Solid angle pixel map queries
dataset_type="solid_angle_pixel_map"
sapm_collection = "cal-sapm-v2-2025-164"
sapm_s3 = [
    f"nasa-irsa-spherex/{data_release}/{dataset_type}/{sapm_collection}/{k}/"\
    f"{dataset_type}_D{k}_spx_{sapm_collection}.fits" for k in range(1, 7)
]
```

## 6. Search for SPHEREx Spectral Images

+++

We search first search for all the available SPHEREx images around a given position on the sky.

Here, we choose the random position at R.A. = 150.7817877 degrees and Decl. = 3.3502204 degrees, which is close to the COSMOS field. We also choose a search radius of 50 arc-seconds.

```{code-cell} ipython3
coord = SkyCoord(150.7817877, 3.3502204, unit='deg')
search_radius = 50 * u.arcsec
```

We then search for the SPHEREx spectral images that overlap with the position in the search radius defined above.
For this we use the [IRSA module in astroquery](https://astroquery.readthedocs.io/en/latest/ipac/irsa/irsa.html) and the Simple Image Access (SIA) API.

```{code-cell} ipython3
results = Irsa.query_sia(pos=(coord, search_radius), collection=sia_collection)
results_summary(results)
```

Note that we have used the collection `spherex_qr2` in this case. All the available collections are documented at [SPHEREx Data Access: Application Program Interfaces (APIs)](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access#application-program-interfaces-apis).
+++

```{tip}
The IRSA SIA collections can be listed using using the ``list_collections`` method, we can filter on the ones containing "spherex" in the collection name:

    Irsa.list_collections(filter='spherex')
```

+++

+++

Next, we set up the _Firefly_ client and use it to visualize one image with the search coordinates (and search radius) indicated.

There are two ways to run _Firefly_ from the JupyterLab environment:
* Open it in another browser tab
* Open it in another Jupyter Notebook tab

The second option requires setting some environment variables and installing the _Firefly_ JupyterLab extension. See [here](https://github.com/Caltech-IPAC/jupyter_firefly_extensions/blob/master/README.md) for more information on how to do that.
More information can be found in the [_Firefly_ documentation](https://caltech-ipac.github.io/firefly_client/usage/).


```{tip}
If you open _Firefly_ in a new Jupyter Notebook tab, you can drag the _Firefly_ tab to the right to enable split view. That way you can see the _Firefly_ window next to this notebook, which makes it more convenient for visualization.
```

Below, we give the two options to open _Firefly_.

```{code-cell} ipython3
# Set up Firefly client to open in a new browser tab
#fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")

# Set up Firefly client to open in new JupyterLab tab
fc = FireflyClient.make_lab_client()

# Initiate Firefly
fc.reinit_viewer()

# Show image
preview_query(fc, results['access_url'][0], coord, search_radius)
```

Now, we select two images, one each in Detector 2 (D2) and Detector (D4). For simplicity, we just select the first image in Detector 2 and the first in Detector 4.
Detector 2 covers 1.10 - 1.62 $\mu$m and Detector 4 covers 2.42 - 3.82 $\mu$m.

We can use two different methods to obtain the images:
- using a boolean mask selecting `energy_bandpassname` as "SPHEREx_D2"
- using the `astropy.table.TableGroups`

We show both methods below, the first one for the detector 2 spectral image and the second one for the detector 4 spectral image.

```{code-cell} ipython3
# Using a boolean mask
d2_result = results[results['energy_bandpassname'] == 'SPHEREx-D2'][0]  # first D2 image

# Using `astropy.table.TableGroups`
results_by_detector = results.group_by('energy_bandpassname')
print(results_by_detector.groups.keys)
d4_result = results_by_detector.groups[3][0]  # first D4 image
```

## 7. Load the Data

Next, we download the two SPHEREx spectral images. We download them directly from the cloud.

For this, we first define a handy function to obtain the S3 cloud path.

```{code-cell} ipython3
# Function to access data from cloud
def get_s3_fpath(cloud_access):
    cloud_info = json.loads(cloud_access) # converts str to dict
    bucket_name = cloud_info['aws']['bucket_name']
    key = cloud_info['aws']['key']
    return f'{bucket_name}/{key}'
```

Next we download the images.

```{code-cell} ipython3
# HDUs for D2 image
s3_fpath_d2 = get_s3_fpath(d2_result['cloud_access'])
with fits.open(f's3://{s3_fpath_d2}', fsspec_kwargs={"anon": True}) as hdul2:
    img_hdu2 = hdul2["IMAGE"].copy()
    img_data2 = hdul2["IMAGE"].data
    img_hdr2 = hdul2["IMAGE"].header
    flags_data2 = hdul2["FLAGS"].data

# HDUs for D4 image
s3_fpath_d4 = get_s3_fpath(d4_result['cloud_access'])
with fits.open(f's3://{s3_fpath_d4}', fsspec_kwargs={"anon": True}) as hdul4:
    img_hdu4 = hdul4["IMAGE"].copy()
    img_data4 = hdul4["IMAGE"].data
    img_hdr4 = hdul4["IMAGE"].header
    flags_data4 = hdul4["FLAGS"].data
```

```{note}
If you want to download the SPHEREx spectral images from IPAC directly (not from the cloud), simply use instead, for example,


    with fits.open(d2_result['access_url']) as hdul2:
        ...


```

+++

## 8. Remove Local Background, Create Masks, and Reproject

Before we can subtract the images, we have to do some preprocessing. This includes:
- Converting from surface brightness to flux density units
- Removing the background
- Creating masks to mask out bad pixels
- Reproject the images on a common pixel grid (necessary for pixel-by-pixel subtraction of images later)

First, we convert from surface brightness to flux units using the solid angle pixel map for each detector. We download these calibrations from the cloud using the S3 paths hardcoded below. Then, we use `astropy.units` to convert the image data from MJy/sr to $\mu$Jy/arcsec$^2$. This quantity is multiplied by the solid angle pixel map in arcsec$^2$ to get image data in units of $\mu$Jy.

```{tip}
You can find more calibration products for the SPHEREx spectral images on the S3 bucket:
[https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/index.html](https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/index.html)

```

```{code-cell} ipython3
s3_fpath_sapm2 = sapm_s3[1] # list idx 1 is value for D2 path
with fits.open(f's3://{s3_fpath_sapm2}', fsspec_kwargs={"anon": True}) as hdul_sapm2:
    sapm_data2 = hdul_sapm2["IMAGE"].data
    sapm_hdr2 = hdul_sapm2["IMAGE"].header

s3_fpath_sapm4 = sapm_s3[3] # list idx 3 is value for D4 path
with fits.open(f"s3://{s3_fpath_sapm4}", fsspec_kwargs={"anon": True}) as hdul_sapm4:
    sapm_data4 = hdul_sapm4["IMAGE"].data
    sapm_hdr4 = hdul_sapm4["IMAGE"].header
```

After downloading the solid angle pixel map calibration products for each detector, we convert the images from brightness (MJy/sr) to spectral flux density ($\mu$Jy) and update the corresponding units in the headers of the imates.

```{code-cell} ipython3
print(f"Spectral image BUNIT: {img_hdr2['BUNIT']}")
print(f"Solid angle pixel map BUNIT: {sapm_hdr2['BUNIT']}")

img_data2 = ((img_data2 * Unit(img_hdr2["BUNIT"])).to(u.uJy / u.arcsec**2) * (sapm_data2 * Unit(sapm_hdr2["BUNIT"]))).value
img_hdu2.header["BUNIT"] = "uJy"
img_hdr2["BUNIT"] = "uJy"

img_data4 = ((img_data4 * Unit(img_hdr4["BUNIT"])).to(u.uJy / u.arcsec**2) * (sapm_data4 * Unit(sapm_hdr4["BUNIT"]))).value
img_hdu4.header["BUNIT"] = "uJy"
img_hdr4["BUNIT"] = "uJy"
```

Next, we remove the background from both images. We estimate the background using the `sep` Python package, which allows fitting a non-parametric background to an image in user-defined zones.

```{warning}
`sep` only support native byte order arrays. We therefore have to apply

    img.astype(img.dtype.newbyteorder("="))

before feeding it to `sep`.
```

```{code-cell} ipython3
bkg2 = sep.Background(img_data2.astype(img_data2.dtype.newbyteorder("=")), bw=64, bh=64, fw=11, fh=11)
img_data2 = img_data2 - bkg2

bkg4 = sep.Background(img_data4.astype(img_data4.dtype.newbyteorder("=")), bw=64, bh=64, fw=11, fh=11)
img_data4 = img_data4 - bkg4
```

Next we create masks (which also need to be reprojected). The merged masks will be applied later to the difference image.
Note that flags `2097152` (which corresponds to $2^{21}$, so bit 21), is assigned to pixel which are not contaminated by image errors, transients, cosmic rays, etc.
Note that we ignore the `OVERFLOW` flag, which is coincident with the `SUR_ERROR` flag. See the [SPHEREx bitmask flag description](https://irsa.ipac.caltech.edu/onlinehelp/spherex/spherex/sp.html#explorebitflags) for more information about the different flag bits.
We therefore choose these as good pixels and flag the other ones.
For simplicity, we set the bad pixels to `np.nan`.

```{code-cell} ipython3
mask2 = np.zeros(flags_data2.shape) * np.nan
mask2[(flags_data2 == 2097152) | (flags_data2 == 2) | (flags_data2 == 0)] = 1

mask4 = np.zeros(flags_data4.shape) * np.nan
mask4[(flags_data4 == 2097152) | (flags_data2 == 2) | (flags_data4 == 0)] = 1
```

Finally, we reproject the images and masks to the same pixel grid. We reproject D2 on D4. For this we use `reproject_exact()` (for flux conservation) and `reproject_interp()` (for mask) from the _reproject_ Python package.

```{warning}
- `reproject_exact` conserves the total flux, but it takes slightly longer than non-exact methods. The cell below should take between 2-3 minutes to execute.
- Due to the bilinear method applied in `reproject_interp`, single-pixel NaNs get conservatively extended into multi-pixel NaN regions.
```

```{code-cell} ipython3
img_data2_reproj, footprint = reproject_exact(input_data=(img_data2, img_hdr2), output_projection=img_hdr4)
footprint[footprint == 0] = np.nan  # Set zero values to NaN
mask2_reproj, _ = reproject_interp(input_data=(mask2, img_hdr2), output_projection=img_hdr4)
```

We use _bokeh_ to visualize the two detector images as well as the reprojected D4 image. See the file `spherex_source_discovery_tool/bokeh_viz.py` for more information on the used function `plot_overlap_trio()` to learn how to use _bokeh_.

```{warning}
The initial bokeh plot has a ~20 second rendering penalty due to static asset loading. Subsequent bokeh plots should take less time to display.
```

```{code-cell} ipython3
# Plot data
plot1_info = {
    "img_data": img_data4,
    "title": "D4 Image",
    "units": img_hdu4.header["BUNIT"]
}
plot2_info = {
    "img_data": img_data2_reproj * footprint,
    "title": "D2 Reprojected on D4",
    "units": img_hdu2.header["BUNIT"]
}
plot3_info = {
    "img_data": img_data4 * footprint,
    "title": "D4 Image with D2 Footprint Mask",
    "units": img_hdu4.header["BUNIT"]
}

plot_overlap_trio(plot1_info, plot2_info, plot3_info, clip=True)
```

## 9. Generate Cutouts and Subtracted Color Image

+++

For this example, we only use parts of the image. For this, we first generate cutouts of the D4 image as well as the reprojected D2 image.

+++

Next, we subtract the images. For keeping track of the wavelength and file names, we create two handy dictionaries.

Finally, we use _bokeh_ to show interactive representation of the two images as well as the subtracted image _inside_ this Jupyter Notebook.

```{code-cell} ipython3
# Set up data directory and cutout paths
data_dir = './data/'
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

cut_path2 = os.path.join(data_dir , f"cutout_{get_filename(d2_result['access_url'])}")
cut_path4 = os.path.join(data_dir , f"cutout_{get_filename(d4_result['access_url'])}")

# Generate cutouts from input images
# NOTE: We use D4's WCS for both cutouts since we reprojected
cut2 = Cutout2D(img_data2_reproj, position=coord, size=(150, 150), wcs=WCS(img_hdu4.header))
cut2_hdu = copy.deepcopy(img_hdu4)  # make copy to avoid modifying HDU of full-size image
cut2_hdu.data = cut2.data
cut2_hdu.header.update(cut2.wcs.to_header())
cut2_hdu.writeto(cut_path2, overwrite=True)
cut2_mask = Cutout2D(mask2_reproj, position=coord, size=(150, 150), wcs=WCS(img_hdu4.header))

cut4 = Cutout2D(img_data4*footprint, position=coord, size=(150, 150), wcs=WCS(img_hdu4.header))
cut4_hdu = copy.deepcopy(img_hdu4)  # make copy to avoid modifying HDU of full-size image
cut4_hdu.data = cut4.data
cut4_hdu.header.update(cut4.wcs.to_header())
cut4_hdu.writeto(cut_path4, overwrite=True)
cut4_mask = Cutout2D(mask4, position=coord, size=(150, 150), wcs=WCS(img_hdu4.header))
```

```{code-cell} ipython3
d2_range = get_lambda_range(d2_result)
d2_info = {
    "img_data": cut2.data,
    "title": rf"{d2_range[0]:.2f}-{d2_range[1]:.2f} um (D2)",
    "units": img_hdr2["BUNIT"]
}

d4_range = get_lambda_range(d4_result)
d4_info = {
    "img_data": cut4.data,
    "title": rf"{d4_range[0]:.2f}-{d4_range[1]:.2f} um (D4)",
    "units": img_hdr4["BUNIT"]
}

subtracted_info = {
    "img_data": np.subtract(cut2.data, cut4.data),
    "title": "Subtracted (D2-D4)",
    "units": img_hdr4["BUNIT"]
}

subtracted_masked_info = {
    "img_data": np.subtract(cut2.data, cut4.data) * cut2_mask.data * cut4_mask.data,
    "title": "Subtracted (D2-D4)",
    "units": img_hdr4["BUNIT"]
}

plot_subtracted_trio(d2_info, d4_info, subtracted_masked_info, clip=True)
```

## 10. Run <code>SExtractor</code> on the Subtracted Color Image

A use case could be to search for red sources such as galaxies with an Active Galactic Nucleus (AGN), which are manifested by red mid-infrared colors. Such red sources can be identified on the subtracted image.

In the following, we run [<code>SExtractor</code>](https://sextractor.readthedocs.io/en/latest/Introduction.html) on all the images (reprojected D2, D4, and subtracted image) to identify sources. Before that, however, we have to save the subtracted image to disk (note that the D2 and D4 cutouts are already saved to disk).

```{code-cell} ipython3
# Save subtracted image to fits file for source extraction
subtracted_cut_path = os.path.join(data_dir , f"subtracted_{get_exp_id(img_hdu2)}_{get_exp_id(img_hdu4)}.fits")
sub_hdu = copy.deepcopy(cut2_hdu)  # make copy to avoid modifying original
sub_hdu.data = subtracted_info["img_data"]
sub_hdu.writeto(subtracted_cut_path, overwrite=True)
```

Now, we run <code>SExtractor</code>. This software needs a few parameter files such as a configuration and parameter files as well as some others (for example a convolution filter). These files are provided in the `spherex_source_discovery_tool` directory. For <code>SExtractor</code> to find them, we have to give it absolute path names to these files. We use `os.path.realpath()` in the following to translate relative paths to absolute paths.
Finally, we run <code>SExtractor</code> on the last lines.

```{code-cell} ipython3
sxt_config = get_sextractor_file("default_sdt.sex", package="spherex_source_discovery_tool")
sxt_params = get_sextractor_file("default_sdt.param", package="spherex_source_discovery_tool")
sxt_nnw = get_sextractor_file("default.nnw", package="spherex_source_discovery_tool")
sxt_conv = get_sextractor_file("default.conv", package="spherex_source_discovery_tool")

subtracted_cat = os.path.realpath( os.path.join(data_dir , f"subtracted_{get_exp_id(img_hdu2)}_{get_exp_id(img_hdu4)}.cat")) # need to use absolute path here!
d2_cat = os.path.realpath( os.path.join(data_dir , f"{get_exp_id(img_hdu2)}.cat") ) # need to use absolute path here!
d4_cat = os.path.realpath( os.path.join(data_dir , f"{get_exp_id(img_hdu4)}.cat") ) # need to use absolute path here!

run_sextractor(os.path.realpath(subtracted_cut_path), sxt_config, sxt_params, subtracted_cat, sxt_nnw, sxt_conv)
run_sextractor(os.path.realpath(cut_path2), sxt_config, sxt_params, d2_cat, sxt_nnw, sxt_conv)
run_sextractor(os.path.realpath(cut_path4), sxt_config, sxt_params, d4_cat, sxt_nnw, sxt_conv)
```

## 11. Select Sources Interactively

In the following, we show two ways to visualize the extracted sources on the images and how to interactively select interesting sources. (Note that you can always select the source automatically in Python using, for example, `astropy` tables instead of doing it interactively.)

First, we read in the <code>SExtractor</code> catalog. Note that we want to change `ALPHA_J2000` and `DELTA_J2000` to the more universal `ra` and `dec` column names. This is especially important for Firefly to be able to read the table out of the box and mark the source in the table on the images.

```{code-cell} ipython3
extracted = Table.read(subtracted_cat , format="ascii")
extracted.rename_column("ALPHA_J2000", "ra")
extracted.rename_column("DELTA_J2000", "dec")
extracted.sort(keys="NUMBER")
```

### 11.1 Using _bokeh_ visualizations

First, we demonstrate how we can use _bokeh_ to create an interactive dashboard to select sources.

```{warning}
Currently, _bokeh_ tables cannot be displayed in Fornax. Therefore the example below will not work. The user is encouraged to use _Firefly_, which use is demonstrated in Section 10.2 below.
```

```{code-cell} ipython3
plot_apertures(subtracted_cut_path, source_tab=extracted, label="Extracted Sources", clip=True)
```

Alternatively, the sources can be easily selected in Python itself by defining a selection function.

```{code-cell} ipython3
# Create custom selection function
def select_sources(sxt_tab):
    """Selects the 25 brightest sources from the SExtractor-generated catalog.

    Parameters
    ----------
    sxt_tab: astropy.table.Table
        Table of extracted sources, including named IDs and converted fluxes

    Returns
    -------
    bright_tab: astropy.table.Table
        Table of selected sources
    """
    # Select bright sources (using automatic aperture flux for now)
    sxt_tab.sort("FLUX_AUTO", reverse=True)
    bright_tab = sxt_tab[0:25]

    return bright_tab

# Preview selected sources
selected = select_sources(extracted)
plot_apertures(subtracted_cut_path, source_tab=selected, label="25 Brightest Sources", clip=True)
```

### 11.2 Using _Firefly_ visualization

Alternatively to _bokeh_, we can use the _Firefly_ application to show the image including the data table. First we reinitiate _Firefly_, which will again open a new tab called "Firefly Viewer".
As mentioned above, we give again two options for initiating _Firefly_; in a new browser tab or a new JupterLab tab.

```{code-cell} ipython3
# Set up Firefly client in browser tab
#fc = FireflyClient.make_client(url="https://irsa.ipac.caltech.edu/irsaviewer")

# Set up Firefly client in new JupyterLab tab
fc = FireflyClient.make_lab_client()

fc.reinit_viewer()
```

Then, we load the images and finally the catalog to the application. _Firefly_ takes care of the rest!

```{tip}
We do not have to save the table to disk. Instead, we can use a "stream" and directly feed it to _Firefly_. We do this by using

    io.BytesIO()

For more tips and tricks on how to use _Firefly_ in Python, we encourage reviewing the [_Firefly_ manual webpage](https://caltech-ipac.github.io/firefly_client/usage/).
```

```{tip}
You can align and lock the images by their WCS by adding

    fc.align_images(lock_match=True)

```

```{code-cell} ipython3
# Load images
fc.show_fits_image(
        file_input=cut_path2,
        plot_id="D2",
        Title="Detector 2"
    )
fc.show_fits_image(
        file_input=cut_path4,
        plot_id="D4",
        Title="Detector 4"
    )
fc.show_fits_image(
        file_input=subtracted_cut_path,
        plot_id="subtracted",
        Title="Subtracted image"
    )
fc.align_images(lock_match=True)

# Load table
tbl_stream = io.BytesIO()
extracted.write(tbl_stream, format='votable')
fc.show_table(file_input=tbl_stream,
              tbl_id="SExtractor Table",
              title='Table')
```

## 12. Perform Aperture Photometry on Selected Sources

Finally, we compute the SPHEREx spectra for selected sources. This includes downloading all the available SPHEREx spectral image cutouts for the sources and then measuring their photometry to generate a spectrum. For simplicity, we here perform aperture photometry but include background subtraction.

### 12.1 Download SPHEREx spectral image cutouts

Here we select one source to proceed (multiple sources can be measured by looping over the code below). Note that we here directly give it the unique source ID that we defined in the catalog column `NUMBER`.

```{code-cell} ipython3
selected_ids = [65]
matched_stars = extracted[extracted['NUMBER'] == selected_ids[0]]["NUMBER", "ra", "dec"].to_pandas()
```

We then set up the TAP service as well as the cutout size to run a query for all SPHEREx spectral images that contain the sources.
```{note}
Note that the following procedure follows the one shown in the SPHEREx [Image Cutout tutorial notebook](https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/).
```

```{code-cell} ipython3
# Define the service endpoint for IRSA's Table Access Protocol (TAP)
# so that we can query SPHEREx metadata tables.
service = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

# Specify cutout size (approximately 60 pixels on a side)
size = 0.1 * u.degree

# Define a query that will search the appropriate SPHEREx metadata tables
# for spectral images that cover the chosen coordinates and match the
# specified bandpass. Return the cutout data access URL and the time of observation.
# Sort by observation time.
query = f"""
SELECT
    'https://irsa.ipac.caltech.edu/' || a.uri || '?center={matched_stars['ra'][0]},{matched_stars['dec'][0]}d&size={size.value}' AS uri,
    p.time_bounds_lower
FROM spherex.artifact a
JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', {matched_stars['ra'][0]}, {matched_stars['dec'][0]}), p.poly)
ORDER BY p.time_bounds_lower
"""

# Execute the query and return as an astropy Table.
t1 = time.time()
tap_results = service.search(query)
print("Time to do TAP query: {:2.2f} seconds.".format(time.time() - t1))
print("Number of images found: {}".format(len(tap_results)))
```

We then define a handy function that downloads the images.

```{code-cell} ipython3
def process_cutout(row, ra, dec, cache, retries=3, delay=5):
    '''
    Downloads the cutouts given in a row of the table including all SPHEREx images overlapping with a position.

    Parameters
    ----------
    row : `astropy.table.Row`
        Row of a table that will be changed in place by this function. The table
        is created by the SQL TAP query.
    ra, dec : `float`
        RA and Dec coordinates (same as used for the TAP query) in degrees.
    cache : `bool`
        If set to `True`, the output of cached and the cutout processing will run faster next time.
        Turn this feature off by setting `cache = False`.
    retries : `int`, optional
        Number of times to retry downloading the cutout in case of failure.
    delay : `int`, optional
        Delay in seconds between retries.
    '''
    for attempt in range(retries):
        try:
            with fits.open(row["uri"]) as hdulist:
                # There are seven HDUs:
                # 0 contains minimal metadata in the header and no data.
                # 1 through 6 are: IMAGE, FLAGS, VARIANCE, ZODI, PSF, WCS-WAVE
                header = hdulist[1].header

                # Fetch spatial and spectral WCS
                spatial_wcs = WCS(header)
                spectral_wcs = WCS(header, fobj=hdulist, key="W")
                spectral_wcs.sip = None  # disable SIP distortions
                row["spectral_wcs"] = spectral_wcs

                # Compute wavelength/bandpass at source position
                x_src, y_src = spatial_wcs.world_to_pixel(SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs"))
                wl_src, bp_src = spectral_wcs.pixel_to_world(x_src, y_src)
                row["central_wl_src"] = wl_src.to(u.micrometer).value

                # Collect the HDUs for this cutout and append the row's cutout_index to the EXTNAME.
                hdus = []
                for hdu in hdulist[1:]:  # skip the primary header
                    hdu.header["EXTNAME"] = f"{hdu.header['EXTNAME']}{row['cutout_index']}"
                    hdus.append(hdu.copy())  # Copy so the data is available after the file is closed
                row["hdus"] = hdus
        except TimeoutError:
            if attempt == retries - 1:
                raise
            time.sleep(delay)
```

We add some placeholders and other information to the table so we can keep track of things.

```{code-cell} ipython3
# Add columns to hold relevant cutout information
results_table = tap_results.to_table()
results_table["cutout_index"] = range(1, len(results_table) + 1)
results_table["central_wl_src"] = np.full(len(results_table), np.nan)
results_table["hdus"] = np.full(len(results_table), None)
results_table["spectral_wcs"] = np.full(len(results_table), None)  # since we are not saving primary HDU
```

Finally, we can run the download of the cutouts in parallel via the function created above. The cutouts are saved directly in the table.

```{warning}
The cell below can take a couple of minutes to run depending on the number of observations available at a given sky coordinate position.
```

```{code-cell} ipython3
# Process cutouts in parallel
t1 = time.time()
n_timeouts = 0
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(process_cutout, row, matched_stars['ra'][0], matched_stars['dec'][0], False) for row in results_table]

    # Raise any exceptions caught during cutout processing
    for future in concurrent.futures.as_completed(futures):
        try:
            future.result()
        except Exception as e:
            n_timeouts += 1
            print(f"An error occurred: {e}")

# Check if any cutouts failed. If so, you might need to re-run this cell with
# fewer workers, cache enabled, or increased timeout/retries.
n_failed = np.count_nonzero(results_table['hdus'] == None)
if n_timeouts > 0 or n_failed > 0:
    print(f"\n!!!!! WARNING !!!!!\n{n_timeouts} cutouts had timeouts and {n_failed} cutouts"
          f" failed to process.\n")

print("Time to create cutouts in serial mode: {:2.2f} minutes.".format((time.time() - t1) / 60))
```

### 12.2 Obtain spectrum

Having obtained the cutouts, we can now perform aperture photometry on the cutouts and measure the fluxes as a function of wavelength.
See `spherex_source_discovery_tool/aperture_photometry.py` for more information on the functions used below.

```{code-cell} ipython3
phot_path = os.path.realpath( "./data/ap_phot_sources.pkl" ) # this is where the photometry results are saved.
aperture_radii = [3.0, 4.0, 5.0]  # aperture radii in pixels
matched_stars = matched_stars.reset_index(drop=True)

# Perform aperture photometry
build_phot_table(results_table, matched_stars, aperture_radii, phot_path, sapm_s3)
```

Load the photometry file we just created.

```{code-cell} ipython3
if os.path.exists(phot_path):
    with open(phot_path, 'rb') as f:
        all_stars_band_phot = pickle.load(f)
```

Apply flagging to the extracted 1D-Spectrum.

```{code-cell} ipython3
# Define bad photometry flags
bad_flags = ['TRANSIENT', 'OVERFLOW', 'SUR_ERROR', 'PHANTOM', 'REFERENCE', 'NONFUNC', 'DICHROIC',
             'MISSING_DATA', 'HOT', 'COLD', 'FULLSAMPLE', 'PHANMISS', 'NONLINEAR', 'PERSIST', 'OUTLIER']
ap_rad = 3.0  # set your preferred aperture radius, from the `aperture_radii` list you specified above
targ = 0  # pick index from `all_stars_band_phot` table
wvl, bdw, abmags, abmag_errs, flx, var, n_bad = \
                grab_star(obs_log_entry=all_stars_band_phot['obs_log'][targ], ap_rad=ap_rad, bad_flags=bad_flags)
```

Finally visualize the spectrum using _bokeh_ in flux and magnitudes!

```{code-cell} ipython3
plot_spectrum(wvl, bdw, flx, var, abmags, abmag_errs, n_bad, type="flux",
              source_id=all_stars_band_phot['NUMBER'][targ], ap_size=3.0)
plot_spectrum(wvl, bdw, flx, var, abmags, abmag_errs, n_bad, type="magnitude",
              source_id=all_stars_band_phot['NUMBER'][targ], ap_size=3.0)
```

## Acknowledgements

- The Source Discovery Tool aperture photometry functions are adapted from code by Zafar Rustamkulov.
- [Caltech/IPAC-IRSA](https://irsa.ipac.caltech.edu/)

This work has made use of the NASA/IPAC Infrared Science Archive, which is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology. We acknowledge support from the SPHEREx project under a contract from the NASA/Goddard Space Flight Center to the California Institute of Technology. This research was partly carried out at the California Institute of Technology under a contract with the National Aeronautics and Space Administration (80GSFC18C0011).  The DOI for SPHEREx QR2 data is 10.26131/IRSA629 (https://doi.org/10.26131/IRSA629).

+++

## About this Notebook

**Contact:** [IRSA Helpdesk](https://irsa.ipac.caltech.edu/docs/help_desk.html) with questions
or problems.

**Updated:** 18 February 2026

**Runtime:** approximately 5 minutes
